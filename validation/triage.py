"""Fast gate triage for the gold-standard cases — no docking.

A full gold-standard run costs ~2 hours, almost all of it Vina search. But most
questions about *why* a case did or did not produce a delta are settled before
docking starts: the covalent-mechanism gate, which structure was chosen, which
protomer, and whether the mutation is inside the 8 A contact range. Those are
arithmetic on coordinates and take seconds.

This reports, per case, the gate it hits and the outcome that implies. Use it to
check a change immediately; use run_gold_standard.py when you need the number.

    python validation/triage.py                # cached structures only, seconds
    python validation/triage.py --fetch        # fetch any PDB not cached yet
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import httpx  # noqa: E402

from secondlook.covalent import classify_covalent  # noqa: E402
from secondlook.validation import GOLD_STANDARD_CASES  # noqa: E402
from secondlook.vina_dock import (  # noqa: E402
    MUTATION_CONTACT_MAX_ANGSTROM,
    NoBindingSiteError,
    chains_with_ligand,
    chains_with_residue,
    extract_chain,
    residue_min_distance_to_ligand,
    select_docking_chain,
)

CACHE = REPO_ROOT / "validation" / "cache"
PDB_CACHE = Path("/tmp/athena-pdb-cache")

BOLD, DIM, RESET = "\033[1m", "\033[2m", "\033[0m"
GREEN, YELLOW, RED = "\033[32m", "\033[33m", "\033[31m"


def cached_structure_id(gene: str, mutation: str) -> str | None:
    path = CACHE / f"{gene}_{mutation}.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    return ((payload.get("structure") or {}).get("id")) or None


def pdb_text(pdb_id: str, fetch: bool) -> str | None:
    PDB_CACHE.mkdir(parents=True, exist_ok=True)
    path = PDB_CACHE / f"{pdb_id}.pdb"
    if path.exists():
        return path.read_text()
    if not fetch:
        return None
    text = httpx.get(
        f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60.0, follow_redirects=True
    ).text
    path.write_text(text)
    return text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fetch", action="store_true", help="download PDBs not already cached")
    args = ap.parse_args()

    print(f"{BOLD}Gate triage — gold-standard cases (no docking){RESET}")
    print(f"{DIM}contact range {MUTATION_CONTACT_MAX_ANGSTROM} A; "
          f"a mutation beyond it cannot be scored by docking{RESET}\n")
    header = (f"{'case':<26}{'drug':<13}{'PDB':<6}{'chain':<6}{'closest':>9}"
              f"  {'protomers':<22}{'gate'}")
    print(header)
    print("-" * len(header))

    would_dock = refused = unknown = 0
    for case in GOLD_STANDARD_CASES:
        label = f"{case.gene} {case.mutation}"
        position = int(case.mutation[1:-1])

        verdict = classify_covalent(case.drug)
        if verdict.is_covalent:
            print(f"{label:<26}{case.drug:<13}{'—':<6}{'—':<6}{'—':>9}  "
                  f"{'—':<22}{YELLOW}covalent — refused before structure{RESET}")
            refused += 1
            continue

        pdb_id = cached_structure_id(case.gene, case.mutation)
        if not pdb_id:
            print(f"{label:<26}{case.drug:<13}{'?':<6}{'?':<6}{'?':>9}  "
                  f"{'—':<22}{DIM}no cached run — run run_gold_standard.py{RESET}")
            unknown += 1
            continue

        text = pdb_text(pdb_id, args.fetch)
        if text is None:
            print(f"{label:<26}{case.drug:<13}{pdb_id:<6}{'?':<6}{'?':>9}  "
                  f"{'—':<22}{DIM}PDB not cached — re-run with --fetch{RESET}")
            unknown += 1
            continue

        common = sorted(set(chains_with_residue(text, position)) & set(chains_with_ligand(text)))
        per_chain = {
            ch: residue_min_distance_to_ligand(extract_chain(text, ch), position) for ch in common
        }
        shown = " ".join(f"{ch}:{d:.1f}" for ch, d in per_chain.items() if d is not None) or "—"
        try:
            chain = select_docking_chain(text, position)
        except NoBindingSiteError:
            print(f"{label:<26}{case.drug:<13}{pdb_id:<6}{'—':<6}{'—':>9}  "
                  f"{shown:<22}{RED}no chain holds both residue and ligand{RESET}")
            refused += 1
            continue

        distance = per_chain.get(chain)
        if distance is None:
            print(f"{label:<26}{case.drug:<13}{pdb_id:<6}{chain:<6}{'—':>9}  "
                  f"{shown:<22}{RED}no measurable distance{RESET}")
            unknown += 1
            continue

        inside = distance <= MUTATION_CONTACT_MAX_ANGSTROM
        gate = (f"{GREEN}within {MUTATION_CONTACT_MAX_ANGSTROM:.0f} A — would dock{RESET}"
                if inside else
                f"{YELLOW}outside {MUTATION_CONTACT_MAX_ANGSTROM:.0f} A — refused{RESET}")
        would_dock += inside
        refused += not inside
        disagree = ""
        if len({round(d, 1) for d in per_chain.values() if d is not None}) > 1:
            disagree = f"  {DIM}(protomers disagree){RESET}"
        print(f"{label:<26}{case.drug:<13}{pdb_id:<6}{chain:<6}{distance:>9.2f}  "
              f"{shown:<22}{gate}{disagree}")

    print()
    print(f"  would dock: {would_dock}   refused at a gate: {refused}   unknown: {unknown}")
    print(f"  {DIM}'would dock' is an upper bound on how many cases can produce a delta.{RESET}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
