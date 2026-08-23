# Vina vs mCSM-lig: what actually scored, and what that costs

`tier2-implementation-spec.md` deliverable 5. Previously unwritable because
nothing scored; the numbers below are measured from the run that produced
`validation/results.md`.

## The measured split

**Every scored result came from AutoDock Vina. mCSM-lig produced none.**

| Method | Scored | Share |
|---|---|---|
| docking (AutoDock Vina) | 2 | **100%** |
| mCSM-lig | 0 | 0% |

Two of nine drug-pairs scored at all; the other seven returned proximity only.

## Why this matters more than it looks

The two methods do not carry the same evidential weight, and the pipeline is
explicit about it:

| | mCSM-lig | AutoDock Vina |
|---|---|---|
| `confidence` | `moderate` | `low` |
| Published accuracy for this task | correlation up to rho = 0.67 | **none** |
| Trained on | Platinum (mutation-induced affinity change) | general protein-ligand docking |
| Answers | "how does this mutation change binding?" | "how well does this ligand dock here?" |

mCSM-lig was built for the question we are asking. Vina was not — it scores a
pose, and we infer a mutation effect by scoring twice and subtracting. The
subtraction is our construction, not something Vina was validated to support.

**So the one figure people will reach for — rho = 0.67 — describes a method that
contributed nothing to these results.** Quoting it alongside a Vina-derived
delta would borrow credibility across methods. `labeling.py` prevents this
structurally: `DOCKING_CALIBRATION` carries `confidence="low"` and its own
`accuracy_note` that does not mention 0.67, and both ride on every label.

## Why mCSM-lig scored nothing

It is a web service driven through a browser form, and it needs a PDB code, a
chain, and the **ligand's HET code** — the crystallised ligand, not the drug we
are asking about. Seven of our nine pairs measure to a ligand that is not the
candidate drug (`ligand_is_candidate_drug: false` throughout), so the question
mCSM-lig would answer is "how does this mutation change binding of the ligand
that happens to be in the crystal", which is not the question asked.

That is a real limitation of the pairing, not a bug to fix by submitting
anyway. Submitting anyway would return a confident number about the wrong
molecule.

## What this means for reporting

1. **Do not quote rho = 0.67 anywhere near these results.** It belongs to a
   method that produced none of them.
2. **The scored deltas are `confidence: low` and `calibration_status:
   provisional`**, and both fields travel into the graph on every signal.
3. **A 100% docking split is not a preference, it is a fallback.** The stronger
   method was unavailable for structural reasons, and the weaker one filled in.
   Anyone reading "2/2 docking" should read it as "the better tool could not be
   used here".

## What would change this

Getting mCSM-lig to contribute needs structures whose crystallised ligand *is*
the candidate drug. That is a structure-selection problem — prefer entries whose
HET code matches the drug's InChIKey connectivity layer — and
`ligand_identity.py` already computes exactly that match. It is currently used
to *label* whether the ligand is the candidate drug, not to *rank* candidate
structures. Wiring it into `rcsb.py`'s ranking is the concrete next step.
