# Proximity validation plan — pre-registered

**Written 2026-08-22, before any evaluation was run.** Criteria below are fixed
at the moment of writing and are not to be adjusted after seeing results, under
the same rule that governs `validation-plan.md`.

## Why this document exists

`validation-plan.md` pre-committed criteria for a **binding-change prediction**:
70% directional agreement across nine cases, with two hard-required positive
controls. That evaluation ran and **failed** — 0/9, both controls failed, zero
candidates scored by any method. That result stands and is not revisited here.

Per that plan's own fallback clause, the product falls back to a narrower claim:
binding-site proximity plus the AlphaMissense flag, with the binding delta and
label dropped from the UI. **That narrower claim has never been validated.** It
is what the pipeline currently ships, so it needs criteria of its own.

This document supplies them.

## What is and is not already known

Honesty about pre-registration, since the distances are not secret:

- **Already known:** the nine measured distances (3.5, 3.5, 10.1, 3.5, 16.6,
  20.0, 12.2, 3.7, 6.9 Å) and the fact that they reproduce exactly across runs.
- **Not known, and what this evaluates:** whether those numbers are *correct* —
  whether they match an independent recomputation, whether the residue measured
  from is the right residue, and whether the ligand measured to is the right
  ligand.

A number can be perfectly reproducible and still wrong. Reproducibility was
established; accuracy was not. That gap is the subject here.

## The claim being tested

> For a given gene, protein mutation and drug, the pipeline reports the distance
> between the mutated residue and a co-crystallised ligand in an experimental
> structure, correctly and with the ligand's identity stated.

Note what is **not** claimed: that the distance predicts resistance, sensitivity,
or the magnitude of any binding change. Case §4 of `../ISSUES.md` and the BRAF
V600E result both show it does not.

## Criteria

All four must pass. A single failure means the proximity claim is not validated
and the fallback is not shippable either.

### C1 — Geometric accuracy (hard requirement)

For every case producing a distance, the reported value must match an
independent recomputation from the raw PDB coordinates to within **0.05 Å**.

The recomputation must be written separately from `proximity.py` and must not
import it — otherwise it re-tests the same arithmetic and proves nothing. It
parses `ATOM`/`HETATM` records directly and takes the minimum heavy-atom
distance between the mutated residue and the ligand group.

### C2 — Residue identity (hard requirement)

For every case, the residue at the reported position in the chosen structure
must match the wild-type residue named in the mutation. `T790M` must measure
from a THR at 790. A mismatch means the structure's numbering does not agree
with UniProt numbering and the distance is measured from the wrong atom.

### C3 — Ligand identity is stated and correct (hard requirement)

For every case, the output must name the ligand actually measured, and that
ligand must be present in the structure under the reported identifier.

This exists because it has already failed once: an earlier build reported 3.5 Å
for both gefitinib and osimertinib from structure 8A27, whose bound ligand is
KY9 — neither drug. The distance was real; the attribution was not.

**A case that measures to a ligand which is not the candidate drug is not a
failure of C3**, provided the output says so plainly. Silence about the
substitution is the failure.

### C4 — Reproducibility

Two independent runs of the same case must produce distances identical to
**0.01 Å**, and must select the same structure and chain.

## Pass threshold

**4 of 4.** Unlike the binding-prediction plan there is no partial-credit
threshold, because these are correctness properties of a measurement rather than
accuracy of a prediction. A measurement is right or it is not.

## What passing does and does not license

**Passing licenses:** reporting the distance, its band
(`in_contact` / `pocket_adjacent` / `distant`), the structure it came from, and
the ligand it was measured against.

**Passing does not license:** any statement about whether the drug will work,
any resistance or sensitivity label, or any implication that a short distance
means impact and a long distance means none. BRAF V600E measures 12.2 Å from
vemurafenib and has a large measured binding effect (MdrDB, PLX-4720,
−2.58 kcal/mol). Distance and effect are different quantities.

**Failing means:** the proximity fallback is not shippable either, and the
honest product is AlphaMissense pathogenicity plus documented Tier 1 evidence,
with no structural output at all.


---

## Result of the first evaluation (2026-08-22)

Run with `validation/verify_proximity.py`, which re-derives every distance from
raw PDB coordinates and does not import `proximity.py`.

| Criterion | Result |
|---|---|
| C1 geometric accuracy | **FAIL** — 8 of 9 |
| C2 residue identity | PASS — 9 of 9 |
| C3 ligand identity stated and present | PASS — 9 of 9 |
| C4 reproducibility | PASS |

**Overall: FAIL** (threshold was 4 of 4).

### What C1 caught

BRAF V600E: the pipeline reported **12.2452 Å**, independent recomputation gave
**7.2077 Å** — 100x the 0.05 Å tolerance.

The cause was not arithmetic. 8C7X is a homodimer, and the two protomers give
12.2452 Å (chain A) and 7.2077 Å (chain B). Both same-chain, both correct
measurements of different copies. The pipeline measured in whichever chain the
*docking* chain-selector returned — a rule written for a task that needs exactly
one receptor — and reported it as the structure's answer. Nothing recorded that
the other copy disagreed by 5 Å, or that the disagreement crossed a band
boundary and changed `distant` into `pocket_adjacent`.

The other eight measurements were confirmed exact to four decimal places.

### Fix

`measure_proximity` now measures in every chain holding both the residue and a
ligand, reports the closest approach, retains all per-chain values in
`ProximitySignal.protomer_distances`, and states the disagreement in the
description when protomers fall in different bands. Seven tests cover it.

Re-verification after the fix is required before this criterion can be marked
passed. That re-run is recorded below.

### Why this vindicates the criterion

C1 exists because a number can be reproducible and still wrong. BRAF V600E
reproduced exactly across two full pipeline runs at two different pipeline
versions, and was measuring the wrong protomer throughout. Any check built on
re-running and comparing would have confirmed it every time.

That is the argument for independent re-derivation over repetition, and it is
worth remembering the next time a validation shortcut looks tempting.

---

## Result of the second evaluation (2026-08-23)

Re-run of `validation/verify_proximity.py` after the protomer fix, against PDB
files re-downloaded from RCSB rather than the ones cached during the first
evaluation.

| Criterion | Result |
|---|---|
| C1 geometric accuracy | **PASS** — 9 of 9 |
| C2 residue identity | PASS — 9 of 9 |
| C3 ligand identity stated and present | PASS — 9 of 9 |
| C4 cross-run reproducibility | **PASS** — 9 of 9 |

**Overall: PASS** (threshold was 4 of 4). Exit status 0.

BRAF V600E now reports **7.2077 Å**, matching independent recomputation exactly.
All nine distances agree with independently parsed coordinates to four decimal
places (delta 0.0000 Å in every row).

### C4 was not actually implemented at the first evaluation

`C4_TOLERANCE_ANGSTROM` was defined but never used; the first evaluation's
"C4 PASS" was recorded without a check behind it. The verifier now writes a
baseline of independently recomputed distances to
`validation/proximity-reference.json` on first run and compares against it on
every subsequent run, and C4 participates in the exit status.

The 9/9 above is therefore the first genuine C4 result: baseline written, PDB
cache deleted, files re-fetched from RCSB, all nine reproduced within 0.01 Å.

### Scope of this pass

This validates the **proximity measurement only** — that the reported distance is
the distance the coordinates support. It says nothing about the docking deltas,
which are a separate measurement with a separate problem; see `ISSUES.md` §15.
