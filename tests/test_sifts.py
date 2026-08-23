"""Tests for SIFTS UniProt<->PDB residue mapping.

Offline: every test injects a fake httpx client. The live checks are marked
`integration`.
"""

from __future__ import annotations

import httpx
import pytest

from secondlook.sifts import MappingSegment, SiftsLookupError, SiftsMapper


def client_returning(payload, status=200):
    def handler(request):
        return httpx.Response(status, json=payload)

    return httpx.Client(transport=httpx.MockTransport(handler))


def entry(mappings, pdb="5hu9", accession="P00519"):
    return {pdb: {"UniProt": {accession: {"identifier": "TEST", "mappings": mappings}}}}


def mapping(unp_start, unp_end, author_start=None, author_end=None, chain="A"):
    return {
        "chain_id": chain,
        "unp_start": unp_start,
        "unp_end": unp_end,
        "start": {"author_residue_number": author_start},
        "end": {"author_residue_number": author_end},
    }


# --- the mapping itself -------------------------------------------------------


def test_maps_a_uniprot_position_to_its_author_residue_number():
    m = SiftsMapper(client=client_returning(entry([mapping(229, 500, author_start=229)])))
    r = m.map_position("5hu9", "P00519", 315)
    assert r is not None
    assert r.author_residue_number == 315
    assert r.offset == 0
    assert r.chain_id == "A"


def test_a_nonzero_offset_is_applied_not_ignored():
    """The case the whole module exists for: PDB numbering that does not match
    UniProt. Measuring at the UniProt number here would hit residue 315, which
    is a real residue and the wrong one."""
    m = SiftsMapper(client=client_returning(entry([mapping(229, 500, author_start=29)])))
    r = m.map_position("5hu9", "P00519", 315)
    assert r.author_residue_number == 115
    assert r.offset == -200


def test_a_position_outside_every_segment_returns_none():
    """Not an error. The residue is genuinely absent from the construct —
    disordered, truncated, or never crystallised."""
    m = SiftsMapper(client=client_returning(entry([mapping(229, 500, author_start=229)])))
    assert m.map_position("5hu9", "P00519", 42) is None
    assert m.is_crystallised("5hu9", "P00519", 42) is False


def test_an_unknown_accession_returns_none_rather_than_guessing():
    m = SiftsMapper(client=client_returning(entry([mapping(229, 500, author_start=229)])))
    assert m.map_position("5hu9", "P99999", 315) is None
    assert m.is_crystallised("5hu9", "P99999", 315) is None


# --- entries that do not supply author numbering ------------------------------


def test_author_start_is_derived_from_author_end_when_missing():
    """BRAF 8C7X populates author_end but not author_start. The segment is
    contiguous, so deriving backwards is exact rather than a guess."""
    m = SiftsMapper(client=client_returning(entry([mapping(444, 721, author_end=721)])))
    r = m.map_position("5hu9", "P00519", 600)
    assert r is not None
    assert r.author_residue_number == 600
    assert r.offset == 0


def test_a_segment_with_no_author_numbering_still_reports_coverage():
    """EGFR 8A27 and ALK 4Z55 return null for both ends. "Present but I cannot
    tell you its number" and "absent" are different facts, and only the second
    should stop a structure being used."""
    m = SiftsMapper(client=client_returning(entry([mapping(700, 1022)])))
    assert m.map_position("5hu9", "P00519", 790) is None
    assert m.is_crystallised("5hu9", "P00519", 790) is True
    assert m.is_crystallised("5hu9", "P00519", 1500) is False


def test_numbering_agrees_separates_false_from_unverifiable():
    """A caller must be able to tell "the assumption is wrong" from "I could not
    check" — collapsing them into False would reject good structures."""
    agrees = SiftsMapper(client=client_returning(entry([mapping(229, 500, author_start=229)])))
    disagrees = SiftsMapper(client=client_returning(entry([mapping(229, 500, author_start=29)])))
    unknown = SiftsMapper(client=client_returning(entry([mapping(700, 1022)])))
    assert agrees.numbering_agrees("5hu9", "P00519", 315) is True
    assert disagrees.numbering_agrees("5hu9", "P00519", 315) is False
    assert unknown.numbering_agrees("5hu9", "P00519", 790) is None


def test_a_segment_missing_its_uniprot_range_is_discarded():
    bad = {"chain_id": "A", "unp_start": None, "unp_end": 500,
           "start": {"author_residue_number": 229}, "end": {}}
    m = SiftsMapper(client=client_returning(entry([bad])))
    assert m.map_position("5hu9", "P00519", 315) is None


# --- transport ----------------------------------------------------------------


def test_a_404_means_no_mapping_not_a_failure():
    """PDBe holds no SIFTS record for some entries. That is a fact about the
    entry, not a transport problem to retry or raise on."""
    m = SiftsMapper(client=client_returning({}, status=404))
    assert m.map_position("9zzz", "P00519", 315) is None


def test_a_server_error_raises_rather_than_reporting_no_mapping():
    """Silently converting an outage into "residue not present" would drop good
    structures and look like a data problem."""
    m = SiftsMapper(client=client_returning({}, status=503))
    with pytest.raises(SiftsLookupError):
        m.map_position("5hu9", "P00519", 315)


def test_one_http_call_per_entry():
    calls = []

    def handler(request):
        calls.append(str(request.url))
        return httpx.Response(200, json=entry([mapping(229, 500, author_start=229)]))

    m = SiftsMapper(client=httpx.Client(transport=httpx.MockTransport(handler)))
    m.map_position("5hu9", "P00519", 315)
    m.map_position("5hu9", "P00519", 316)
    m.is_crystallised("5hu9", "P00519", 400)
    assert len(calls) == 1


def test_segment_maths_is_exact_across_the_whole_range():
    seg = MappingSegment(chain_id="A", unp_start=229, unp_end=500, author_start=29)
    assert seg.to_author(229) == 29
    assert seg.to_author(500) == 300
    assert seg.covers(229) and seg.covers(500)
    assert not seg.covers(228) and not seg.covers(501)


# --- live ---------------------------------------------------------------------


@pytest.mark.integration
def test_live_sifts_confirms_the_numbering_assumption_on_our_structures():
    """The pipeline assumes PDB author numbering equals UniProt numbering. This
    is that assumption checked against the authoritative source, on the exact
    structures the gold-standard run uses."""
    m = SiftsMapper()
    cases = [
        ("5hu9", "P00519", 315), ("8c7x", "P15056", 600),
        ("8pqd", "P10721", 816), ("8s14", "P10721", 560),
    ]
    for pdb, acc, pos in cases:
        assert m.is_crystallised(pdb, acc, pos) is True, f"{pdb} {pos} not crystallised"
        assert m.numbering_agrees(pdb, acc, pos) is True, f"{pdb} {pos} numbering differs"


@pytest.mark.integration
def test_live_sifts_confirms_coverage_where_it_cannot_confirm_numbering():
    """EGFR 8A27 and ALK 4Z55 supply no author numbering. They must still
    confirm the residue is in the structure."""
    m = SiftsMapper()
    for pdb, acc, pos in [("8a27", "P00533", 790), ("4z55", "Q9UM73", 1202)]:
        assert m.is_crystallised(pdb, acc, pos) is True
        assert m.numbering_agrees(pdb, acc, pos) is None
