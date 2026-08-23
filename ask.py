#!/usr/bin/env python3
"""Ask Athena about one gene + mutation + drug.

Tier 1 (documented evidence in the CIViC knowledge graph) is consulted first.
If it holds strong evidence, that IS the answer and Tier 2 does not run --
a published Level A/B finding beats a structural prediction every time.
Only when Tier 1 has nothing, or only weak evidence, does Tier 2 compute a
structural second opinion.

    python ask.py ABL1 T315I --drug imatinib
    python ask.py ABL1 T315I --drug imatinib --cancer-type "chronic myeloid leukemia"
    python ask.py BRAF V600E --drug vemurafenib --tier2-anyway
    python ask.py EGFR T790M --drug osimertinib --skip-tier1

Tier 1 needs FalkorDB (`docker compose up -d`). Tier 2 hits live services
(UniProt, RCSB, DGIdb, ...) and takes minutes per drug.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from redis.exceptions import RedisError  # noqa: E402  (falkordb's transport)
from secondlook.pipeline import run_tier2  # noqa: E402
from secondlook.tier1_adapter import Tier1RetrievalPolicy, Tier1UnavailableError  # noqa: E402

#: FalkorDB speaks the Redis protocol, so a database that is down surfaces as a
#: redis-py error, not as Tier1UnavailableError. Both are caught at the call
#: site -- an unreachable graph is a "start Docker" message, not a traceback.
TIER1_DOWN = (Tier1UnavailableError, RedisError, OSError)

#: "T315I", "V600E", "R273*" -- CIViC's own name for a missense variant.
SHORTHAND = re.compile(r"^[A-Z]\d+[A-Z*]$")

#: Protein-level HGVS, the only form Tier 1's `hgvs_p` column holds:
#: "p.Thr315Ile", or RefSeq-prefixed "NP_005148.2:p.Thr315Ile".
HGVS_P = re.compile(r"(^|:)p\.")

BOLD, DIM, RESET = "\033[1m", "\033[2m", "\033[0m"
GREEN, YELLOW, RED = "\033[32m", "\033[33m", "\033[31m"


ONE_TO_THREE = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}


def to_three_letter(mutation: str) -> str | None:
    """`T315I` -> `p.Thr315Ile`. None if either residue is not standard."""
    m = SHORTHAND.match(mutation)
    if not m:
        return None
    wt, pos, mut = mutation[0], mutation[1:-1], mutation[-1]
    wt3, mut3 = ONE_TO_THREE.get(wt), ONE_TO_THREE.get(mut)
    if wt3 is None or mut3 is None:
        return None
    return f"p.{wt3}{pos}{mut3}"


def variant_keys(mutation: str) -> dict:
    """The single best Tier 1 match key for this mutation string."""
    return tier1_lookups(mutation)[0]


def tier1_lookups(mutation: str) -> list[dict]:
    """Every Tier 1 match key worth trying, best first.

    The loaded graph uses three different conventions, verified against it on
    2026-08-23, so no single key finds everything:

      | variant                  | `name`        | `hgvs_p`                    |
      |--------------------------|---------------|-----------------------------|
      | ERCC2 K751Q              | `K751Q`       | (empty)                     |
      | FGFR1 N546K              | `N546K`       | `NP_075598.2:p.Asn546Lys`   |
      | ABL1 (Tier 2 stub)       | (empty)       | `p.Thr315Ile`               |

    `name` is populated for essentially every variant that carries evidence, so
    shorthand goes first. `hgvs_p` is tried second because exact equality also
    fails on the RefSeq prefix -- that lookup usually only succeeds via Mode 2's
    relaxed fallback, which is still better than reporting nothing.
    """
    if HGVS_P.search(mutation):
        return [{"mutation": mutation}]
    if not SHORTHAND.match(mutation):
        # "Fusion", "Amplification", "Loss", "Overexpression", ... CIViC's names
        # for everything that is not a point substitution. These live in `name`;
        # `hgvs_p` is empty for all of them, so sending them there matches
        # nothing and silently falls through to Mode 2's relaxed gene-level
        # query -- which downgrades a real exact hit to `weak_hit`.
        return [{"variant_name": mutation}]
    lookups = [{"variant_name": mutation}]
    three = to_three_letter(mutation)
    if three:
        lookups.append({"mutation": three})
    return lookups


def rule(title: str) -> None:
    print(f"\n{BOLD}{title}{RESET}\n{'-' * len(title)}")


def show_tier1(decision) -> None:
    colour = {"strong_hit": GREEN, "weak_hit": YELLOW, "no_hit": DIM}.get(decision.state, "")
    print(f"  state   {colour}{decision.state}{RESET}")
    print(f"  reason  {decision.reason}")
    if decision.is_placeholder:
        print(f"  {YELLOW}note    placeholder policy, not live Tier 1 retrieval{RESET}")
    if not decision.tier1_results:
        print(f"  {DIM}no evidence items{RESET}")
        return
    print()
    for item in decision.tier1_results:
        drug = item.drug or "—"
        print(f"  [{item.evidence_level}] {BOLD}{drug}{RESET}  ({item.source} {item.citation_id})")
        print(f"      {item.summary}")
        print(f"      {DIM}{item.citation_url}{RESET}")


def show_tier2(out) -> None:
    print(f"  status  {out.status}")
    if out.structure:
        pdb = getattr(out.structure, "id", None) or "—"
        print(f"  structure  {pdb}")
    if not out.results:
        print(f"  {DIM}no candidate drugs returned{RESET}")
        for f in out.failures:
            print(f"  {RED}failure{RESET}  {getattr(f, 'message', f)}")
        return
    for item in out.results:
        print()
        print(f"  {BOLD}{item.drug}{RESET}   {item.signal_type}")
        if item.signal_type == "binding_delta":
            print(f"      label   {BOLD}{item.label}{RESET}")
            print(f"      delta   {item.delta_score:+.3f} kcal/mol  (method: {item.method})")
        else:
            print(f"      {YELLOW}no binding prediction made{RESET}")
            if item.no_delta_reason_code:
                # The code matters as much as the prose: `docking_timeout` is a
                # setting to raise, `receptor_prep_failed` is a bug to fix, and
                # `outside_pocket` is a real scientific limit. See ISSUES.md #17.
                print(f"      code    {BOLD}{item.no_delta_reason_code}{RESET}")
            if item.no_delta_reason:
                print(f"      why     {item.no_delta_reason}")
        if item.binding_site_distance_angstrom is not None:
            print(f"      pocket  {item.binding_site_distance_angstrom:.2f} A from the bound ligand")
        band = (item.proximity or {}).get("band")
        if band:
            print(f"      band    {band}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("gene", help="HGNC symbol, e.g. ABL1")
    ap.add_argument("mutation", help="shorthand (T315I) or protein HGVS (p.Thr315Ile)")
    ap.add_argument("--drug", action="append", dest="drugs",
                    help="restrict to this drug; repeatable. Omit to score every candidate.")
    ap.add_argument("--cancer-type", default=None, help="disease context for Tier 1")
    ap.add_argument("--tier2-anyway", action="store_true",
                    help="run Tier 2 even if Tier 1 has strong evidence (the UI override)")
    ap.add_argument("--skip-tier1", action="store_true",
                    help="go straight to Tier 2 (use when FalkorDB is not running)")
    args = ap.parse_args()

    drugs = tuple(args.drugs) if args.drugs else None
    print(f"{BOLD}Athena — second opinion{RESET}")
    print(f"  gene      {args.gene}")
    print(f"  mutation  {args.mutation}")
    print(f"  drug      {', '.join(drugs) if drugs else 'all candidates'}")
    if args.cancer_type:
        print(f"  cancer    {args.cancer_type}")

    policy = None
    if not args.skip_tier1:
        rule("TIER 1 — documented evidence (CIViC knowledge graph)")
        policy = Tier1RetrievalPolicy()
        try:
            decision = None
            for keys in tier1_lookups(args.mutation):
                decision = policy.decide(
                    gene=args.gene,
                    cancer_type=args.cancer_type,
                    manual_override=args.tier2_anyway,
                    **keys,
                )
                label = ", ".join(f"{k}={v!r}" for k, v in keys.items())
                found = len(decision.tier1_results)
                print(f"  {DIM}tried {label} -> {decision.state} ({found} item(s)){RESET}")
                if found or args.tier2_anyway:
                    break
        except TIER1_DOWN as exc:
            print(f"  {RED}Tier 1 unavailable{RESET}  {exc}")
            print(f"  {DIM}start it with:  docker compose up -d{RESET}")
            print(f"  {DIM}or skip it with: --skip-tier1{RESET}")
            return 2
        show_tier1(decision)
        if not decision.should_run_tier2:
            print(f"\n  {GREEN}Tier 1 answered. Tier 2 not run.{RESET}")
            print(f"  {DIM}Override with --tier2-anyway.{RESET}")
            return 0
        print(f"\n  {YELLOW}Tier 1 is not sufficient — handing to Tier 2.{RESET}")

    rule("TIER 2 — structural prediction")
    print(f"  {DIM}hitting live services, this takes a few minutes per drug...{RESET}")
    out = run_tier2(
        args.gene,
        args.mutation,
        cancer_type=args.cancer_type,
        restrict_to_drugs=drugs,
        activation_policy=None,   # Tier 1 already ran above; don't ask twice
    )
    show_tier2(out)
    print()
    print(f"  {DIM}A 'proximity_only' result is not a failure — it is the honest{RESET}")
    print(f"  {DIM}output where no defensible binding delta could be computed.{RESET}")
    print(f"  {DIM}See ISSUES.md §1 before quoting any delta.{RESET}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
