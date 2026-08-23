"""SIFTS UniProt <-> PDB residue-number mapping (EBI PDBe).

## Why this exists

Everything upstream numbers residues by **UniProt** position: `T315I` means
threonine at UniProt position 315. Everything downstream — the coordinates, the
distance measurement, the mutation placement — indexes into a **PDB** file,
which numbers residues however the depositor chose. The two agree often enough
to be dangerous and not often enough to assume.

When they disagree, nothing fails. The pipeline measures from a real residue at
a real position that is simply the wrong residue, and reports a confident
distance for it. `structure.covers_residue()` and `vina_dock`'s residue-identity
check catch the subset where the wrong residue happens to be a different amino
acid, which is most of the time but not all: a T->T collision passes silently.

SIFTS is the authoritative, curated mapping maintained by EBI for exactly this
problem. This module reads it so the assumption becomes a check.

## What it does not do

It does not renumber structures or rewrite coordinates. It answers one question
-- "which author residue number in this PDB entry corresponds to UniProt
position N?" -- and says plainly when it cannot.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import httpx

from secondlook.http_retry import with_retry

DEFAULT_PDBE_API_BASE = "https://www.ebi.ac.uk/pdbe/api"


class SiftsLookupError(RuntimeError):
    """SIFTS could not be consulted. Distinct from "consulted, no mapping"."""


@dataclass(frozen=True)
class ResidueMapping:
    """Where one UniProt position lands in one PDB chain."""

    pdb_id: str
    uniprot_accession: str
    chain_id: str
    #: The number written in the PDB file's ATOM records — what a parser sees.
    author_residue_number: int
    uniprot_position: int

    @property
    def offset(self) -> int:
        """author - uniprot. Zero means the two numbering schemes coincide."""
        return self.author_residue_number - self.uniprot_position


@dataclass(frozen=True)
class MappingSegment:
    """One contiguous UniProt->PDB correspondence in a single chain.

    `author_start` is optional because PDBe does not always populate
    `author_residue_number` — EGFR 8A27 and ALK 4Z55 both return null for it
    while still giving a valid `unp_start`/`unp_end` range. Those entries can
    still answer "is this residue crystallised at all?", which is worth having,
    so the segment is kept rather than discarded. It simply cannot answer "what
    number is it called", and says so by returning None.
    """

    chain_id: str
    unp_start: int
    unp_end: int
    author_start: int | None = None

    def covers(self, uniprot_position: int) -> bool:
        return self.unp_start <= uniprot_position <= self.unp_end

    def to_author(self, uniprot_position: int) -> int | None:
        # Segments are contiguous by construction, so one offset is exact within
        # a segment. Across segments it is not, which is why each is kept
        # separately rather than collapsed into a single global offset.
        if self.author_start is None:
            return None
        return self.author_start + (uniprot_position - self.unp_start)


class SiftsMapper:
    """Fetches and applies SIFTS mappings. One HTTP call per PDB entry, cached."""

    def __init__(self, base_url: str | None = None, client: httpx.Client | None = None) -> None:
        self.base_url = (
            base_url or os.environ.get("PDBE_API_BASE") or DEFAULT_PDBE_API_BASE
        ).rstrip("/")
        self._client = client
        self._cache: dict[str, dict[str, list[MappingSegment]]] = {}

    def segments(self, pdb_id: str) -> dict[str, list[MappingSegment]]:
        """All UniProt->PDB segments for `pdb_id`, keyed by UniProt accession."""
        key = pdb_id.lower()
        if key in self._cache:
            return self._cache[key]

        payload = self._get_json(f"/mappings/uniprot/{key}")
        entry = payload.get(key) or {}
        by_accession: dict[str, list[MappingSegment]] = {}
        for accession, body in (entry.get("UniProt") or {}).items():
            segments: list[MappingSegment] = []
            for m in body.get("mappings") or []:
                start, end = m.get("start") or {}, m.get("end") or {}
                unp_start, unp_end = m.get("unp_start"), m.get("unp_end")
                chain = m.get("chain_id")
                # Without the UniProt range and a chain there is nothing to say.
                # Inventing a default here would produce a mapping that looks
                # authoritative and is not.
                if None in (unp_start, unp_end, chain):
                    continue

                author_start = start.get("author_residue_number")
                if author_start is None:
                    # Fall back to deriving the start from the segment end.
                    # BRAF 8C7X populates author_end but not author_start, and
                    # the segment is contiguous, so this is exact rather than a
                    # guess.
                    author_end = end.get("author_residue_number")
                    if author_end is not None:
                        author_start = int(author_end) - (int(unp_end) - int(unp_start))

                segments.append(
                    MappingSegment(
                        chain_id=str(chain),
                        unp_start=int(unp_start),
                        unp_end=int(unp_end),
                        author_start=None if author_start is None else int(author_start),
                    )
                )
            if segments:
                by_accession[accession.upper()] = segments

        self._cache[key] = by_accession
        return by_accession

    def map_position(
        self, pdb_id: str, uniprot_accession: str, uniprot_position: int, chain_id: str | None = None
    ) -> ResidueMapping | None:
        """Author residue number for a UniProt position, or None if unmapped.

        None is a real answer: the position may lie outside every crystallised
        segment (disordered, truncated, or simply absent from the construct).
        Callers must treat it as "SIFTS says this residue is not in this
        structure", not as an error to retry.
        """
        by_accession = self.segments(pdb_id)
        segments = by_accession.get(uniprot_accession.upper())
        if not segments:
            return None
        for seg in segments:
            if chain_id is not None and seg.chain_id != chain_id:
                continue
            if seg.covers(uniprot_position):
                author = seg.to_author(uniprot_position)
                if author is None:
                    continue  # covered, but this entry cannot supply numbering
                return ResidueMapping(
                    pdb_id=pdb_id.upper(),
                    uniprot_accession=uniprot_accession.upper(),
                    chain_id=seg.chain_id,
                    author_residue_number=author,
                    uniprot_position=uniprot_position,
                )
        return None

    def is_crystallised(
        self, pdb_id: str, uniprot_accession: str, uniprot_position: int
    ) -> bool | None:
        """Whether SIFTS places this UniProt position inside the entry at all.

        Answerable even when the entry supplies no author numbering, which is
        the common case. Separating it from `map_position` matters: "the residue
        is absent from this structure" and "the residue is present but I cannot
        tell you its number" are different facts, and only the first should stop
        a structure being used.
        """
        segments = self.segments(pdb_id).get(uniprot_accession.upper())
        if not segments:
            return None
        return any(seg.covers(uniprot_position) for seg in segments)

    def numbering_agrees(
        self, pdb_id: str, uniprot_accession: str, uniprot_position: int
    ) -> bool | None:
        """True when author numbering equals UniProt numbering at this position.

        The pipeline assumes this. Returns None when SIFTS has no mapping, so a
        caller can tell "assumption verified false" from "could not verify".
        """
        mapping = self.map_position(pdb_id, uniprot_accession, uniprot_position)
        return None if mapping is None else mapping.offset == 0

    def _get_json(self, path: str) -> dict:
        url = f"{self.base_url}{path}"

        def attempt() -> httpx.Response:
            if self._client is not None:
                response = self._client.get(url, timeout=30.0)
            else:
                response = httpx.get(url, timeout=30.0, follow_redirects=True)
            response.raise_for_status()
            return response

        try:
            return with_retry(attempt).json()
        except httpx.HTTPStatusError as exc:
            # 404 means PDBe holds no SIFTS mapping for this entry — a fact, not
            # a transport failure, so it is reported as "no mapping" rather than
            # raised.
            if exc.response.status_code == 404:
                return {}
            raise SiftsLookupError(f"SIFTS request failed for {url}") from exc
        except httpx.HTTPError as exc:
            raise SiftsLookupError(f"SIFTS request failed for {url}") from exc
        except ValueError as exc:
            raise SiftsLookupError(f"SIFTS returned malformed JSON for {url}") from exc
