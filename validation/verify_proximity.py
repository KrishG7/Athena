"""Independent verification of the proximity measurement.

Evaluates the four criteria pre-registered in `docs/proximity-validation-plan.md`.

**This module deliberately does not import `secondlook.proximity`.** It parses
PDB `ATOM`/`HETATM` records itself and does its own distance arithmetic. Sharing
code with the thing under test would re-run the same arithmetic and prove
nothing -- the point is an independent second opinion on the same question.

Run:
    python validation/verify_proximity.py --cache ../Athena_tier_2/validation/cache
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import httpx

#: Criterion C1: reported vs independently recomputed distance.
C1_TOLERANCE_ANGSTROM = 0.05
#: Criterion C4: cross-run reproducibility, at the precision earlier runs recorded.
C4_TOLERANCE_ANGSTROM = 0.01

#: Three-letter -> one-letter, for the C2 residue-identity check.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

#: Never treated as the bound drug: water, and the ions/buffer components that
#: co-crystallise incidentally. A distance measured to a sulfate is not a
#: distance to a ligand.
NOT_A_LIGAND = {
    "HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "ACT", "CL", "NA", "MG", "ZN",
    "CA", "K", "MN", "PEG", "DMS", "TRS", "IOD", "BR", "FMT", "NO3", "ACY",
}


@dataclass
class Atom:
    chain: str
    resseq: int
    resname: str
    element: str
    x: float
    y: float
    z: float


def parse_pdb(text: str) -> tuple[list[Atom], list[Atom]]:
    """Return (protein atoms, hetero atoms). Heavy atoms only -- hydrogens are
    frequently absent from crystal structures, so including them would make the
    measurement depend on whether a given entry happens to have them."""
    protein: list[Atom] = []
    hetero: list[Atom] = []
    for line in text.splitlines():
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
        if element == "H":
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):      # primary conformer only
            continue
        try:
            atom = Atom(
                chain=line[21],
                resseq=int(line[22:26]),
                resname=line[17:20].strip(),
                element=element,
                x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
            )
        except ValueError:
            continue
        (protein if rec == "ATOM  " else hetero).append(atom)
    return protein, hetero


def min_distance(a: list[Atom], b: list[Atom]) -> float | None:
    if not a or not b:
        return None
    best = math.inf
    for p in a:
        for q in b:
            d = math.dist((p.x, p.y, p.z), (q.x, q.y, q.z))
            if d < best:
                best = d
    return None if best is math.inf else best


def fetch_pdb(pdb_id: str, cache_dir: Path) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{pdb_id}.pdb"
    if path.exists():
        return path.read_text()
    text = httpx.get(
        f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60.0, follow_redirects=True
    ).text
    path.write_text(text)
    return text


@dataclass
class CaseCheck:
    gene: str
    mutation: str
    drug: str
    pdb_id: str
    position: int
    expected_wt: str
    reported_distance: float
    reported_ligand: str
    ligand_is_candidate: bool
    # filled in by verification
    recomputed: float | None = None
    chain: str | None = None
    observed_residue: str | None = None
    c1: bool = False
    c2: bool = False
    c3: bool = False
    c4: bool | None = None      # None until a baseline exists to compare against
    notes: str = ""

    @property
    def key(self) -> str:
        return f"{self.gene} {self.mutation}/{self.drug}"



def verify(case: CaseCheck, pdb_text: str) -> CaseCheck:
    protein, hetero = parse_pdb(pdb_text)

    # C3 first: is the named ligand actually in this structure?
    lig_atoms_by_chain: dict[str, list[Atom]] = {}
    for h in hetero:
        if h.resname == case.reported_ligand:
            lig_atoms_by_chain.setdefault(h.chain, []).append(h)
    if not lig_atoms_by_chain:
        case.notes = f"ligand {case.reported_ligand!r} not present in {case.pdb_id}"
        return case
    if case.reported_ligand in NOT_A_LIGAND:
        case.notes = f"{case.reported_ligand} is a crystallisation additive, not a drug"
        return case
    case.c3 = True

    # Find the chain where BOTH the residue and the ligand live, and measure
    # there. Scanning rather than trusting a recorded chain is deliberate: it
    # independently re-derives the choice the pipeline made.
    best: tuple[float, str, str] | None = None
    for chain, ligs in lig_atoms_by_chain.items():
        res = [a for a in protein if a.chain == chain and a.resseq == case.position]
        if not res:
            continue
        d = min_distance(res, ligs)
        if d is not None and (best is None or d < best[0]):
            best = (d, chain, res[0].resname)

    if best is None:
        case.notes = (
            f"no chain holds both residue {case.position} and ligand "
            f"{case.reported_ligand} (ligand chains: {sorted(lig_atoms_by_chain)})"
        )
        return case

    case.recomputed, case.chain, case.observed_residue = best
    case.c1 = abs(case.recomputed - case.reported_distance) <= C1_TOLERANCE_ANGSTROM
    case.c2 = THREE_TO_ONE.get(case.observed_residue or "") == case.expected_wt
    if not case.c2:
        case.notes = (
            f"expected {case.expected_wt} at {case.position}, "
            f"structure has {case.observed_residue}"
        )
    return case


def load_cases(cache_dir: Path) -> list[CaseCheck]:
    cases: list[CaseCheck] = []
    for f in sorted(cache_dir.glob("*.json")):
        payload = json.loads(f.read_text())
        struct = payload.get("structure") or {}
        for item in payload.get("results", []):
            prox = item.get("proximity") or {}
            dist = prox.get("distance_angstrom")
            if dist is None:
                continue
            mut = payload["mutation"]
            cases.append(
                CaseCheck(
                    gene=payload["gene"],
                    mutation=mut,
                    drug=item.get("drug", "?"),
                    pdb_id=struct.get("id", ""),
                    position=int(struct.get("annotated_position") or mut[1:-1]),
                    expected_wt=mut[0],
                    reported_distance=float(dist),
                    reported_ligand=prox.get("measured_to_ligand") or "",
                    ligand_is_candidate=bool(prox.get("ligand_is_candidate_drug")),
                )
            )
    return cases


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", default="../Athena_tier_2/validation/cache")
    ap.add_argument("--pdb-cache", default="/tmp/athena-pdb-cache")
    ap.add_argument(
        "--reference",
        default=str(Path(__file__).parent / "proximity-reference.json"),
        help="baseline of previously recomputed distances, for the C4 cross-run check",
    )
    args = ap.parse_args()

    cache_dir = Path(args.cache)
    if not cache_dir.is_dir():
        print(f"cache directory not found: {cache_dir}", file=sys.stderr)
        return 2

    cases = load_cases(cache_dir)
    if not cases:
        print("no cached cases with a distance; run validation/run_gold_standard.py first")
        return 2

    print(f"Verifying {len(cases)} measurements against independently parsed coordinates\n")
    header = f"{'case':<26} {'PDB':<6} {'lig':<5} {'reported':>9} {'recomputed':>11} {'delta':>7}  C1 C2 C3"
    print(header)
    print("-" * len(header))

    for c in cases:
        pdb_text = fetch_pdb(c.pdb_id, Path(args.pdb_cache))
        verify(c, pdb_text)
        rec = f"{c.recomputed:.4f}" if c.recomputed is not None else "—"
        delta = f"{abs(c.recomputed - c.reported_distance):.4f}" if c.recomputed is not None else "—"
        label = f"{c.gene} {c.mutation}/{c.drug}"
        print(
            f"{label:<26} {c.pdb_id:<6} {c.reported_ligand:<5} "
            f"{c.reported_distance:>9.4f} {rec:>11} {delta:>7}  "
            f"{'Y' if c.c1 else 'N'}  {'Y' if c.c2 else 'N'}  {'Y' if c.c3 else 'N'}"
        )
        if c.notes:
            print(f"    note: {c.notes}")

    print()
    # C4: does this run reproduce the distances an earlier run recorded? The
    # baseline is written on first use, so the very first run establishes it
    # rather than pretending to check against itself.
    ref_path = Path(args.reference)
    measured = {c.key: c.recomputed for c in cases if c.recomputed is not None}
    if ref_path.exists():
        reference = json.loads(ref_path.read_text())
        for c in cases:
            prior = reference.get(c.key)
            if prior is None or c.recomputed is None:
                c.c4 = None
            else:
                c.c4 = abs(c.recomputed - float(prior)) <= C4_TOLERANCE_ANGSTROM
                if not c.c4:
                    print(
                        f"    C4 drift: {c.key} was {float(prior):.4f}, now {c.recomputed:.4f}"
                    )
        c4_state = "checked"
    else:
        ref_path.write_text(json.dumps(measured, indent=2, sort_keys=True) + "\n")
        c4_state = "baseline"

    for crit, name, ok in (
        ("C1", f"geometric accuracy (<= {C1_TOLERANCE_ANGSTROM} A)", all(c.c1 for c in cases)),
        ("C2", "residue identity", all(c.c2 for c in cases)),
        ("C3", "ligand identity stated and present", all(c.c3 for c in cases)),
    ):
        n = sum(1 for c in cases if getattr(c, crit.lower()))
        print(f"{crit}  {name:<45} {n}/{len(cases)}  {'PASS' if ok else 'FAIL'}")

    c4_name = f"cross-run reproducibility (<= {C4_TOLERANCE_ANGSTROM} A)"
    if c4_state == "baseline":
        print(f"C4  {c4_name:<45} {'—':>3}  baseline recorded ({ref_path.name}); "
              "re-run to check")
    else:
        comparable = [c for c in cases if c.c4 is not None]
        n = sum(1 for c in comparable if c.c4)
        verdict = "PASS" if comparable and n == len(comparable) else "FAIL"
        if not comparable:
            verdict = "no overlap with baseline"
        print(f"C4  {c4_name:<45} {n}/{len(comparable)}  {verdict}")

    substituted = [c for c in cases if not c.ligand_is_candidate]
    print(
        f"\n{len(substituted)} of {len(cases)} measured to a ligand that is NOT the candidate "
        "drug. Per C3 that is acceptable only because the output says so explicitly "
        "(`ligand_is_candidate_drug: false`, and the ligand named in the description)."
    )
    passed = all(c.c1 and c.c2 and c.c3 for c in cases) and all(
        c.c4 for c in cases if c.c4 is not None
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
