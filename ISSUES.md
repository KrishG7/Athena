# Athena — known issues and how to fix them

**For the team.** Everything here is measured or read out of code. Nothing is
estimated; where something is unverified it says so.

Each section is: what's wrong → why → what to do. Sections are independent, so
you can take one without reading the rest.

**Read §1 before presenting this project to anyone.**

## Status at a glance

| Area | State |
|---|---|
| Both tiers install and import as one package | Working |
| 611 unit tests | Passing, 0 failures (also with the legacy source trees off `sys.path`, so a fresh clone is self-contained) |
| Tier 1 retrieval Modes 1 & 2 against live CIViC graph | Working |
| Tier 1 → Tier 2 activation handoff | Verified live; **variant-key routing downgrades strong hits** — §16 |
| Tier 2 → graph writes, full traversal | Verified live |
| WT-vs-mutant docking runs end to end | Verified (2Z4O D130N) |
| Binding-site proximity | Independently re-derived; all four pre-registered criteria **PASS** 9/9 — §13 |
| Docking protomer selection | Was using the wrong homodimer copy, wrongly refusing BRAF V600E — **fixed**, re-run needed — §18 |
| **Validated binding-change prediction** | **1/9 — but that run docked with a 120s timeout that was killing every attempt (§17). Re-run pending.** |
| Docking reproducibility | Root cause found (PDBFixer/OpenMM randomness) and **fixed**; prep + Vina both verified bit-identical — §15 |
| Labeling calibration | `provisional`, blocked on a run that scores — §3 |
| meeko receptor-prep blocker (3 of 9 cases) | Fixed via Open Babel fallback — §2 |
| Trial source (ClinicalTrials.gov / CTRI) | Not built — §5 |
| Mode 3 semantic retrieval | Never verified end to end — §10 |

## Where to start

| If you have | Work on |
|---|---|
| Half a day and want the biggest scientific win | **§3** — calibrate, now that §2 unblocks three more cases |
| A day and want the project defensible | **§15** — measure the docking noise floor at exhaustiveness 32, then re-run the gold standard |
| Backend/API interest | **§5** — build the ClinicalTrials.gov loader |
| ML interest | **§12 step 5** — MdrDB training run |

---

## 1. The headline: the pipeline does not currently predict binding change

### What the numbers say

From `validation/results.md`, regenerated from a clean full run:

```
Pass rate: 1/9 (11%) — threshold 70% — NOT MET
Hard-required positive controls: FAILED
Demo-ready: NO
scored_by_method: {'docking': 2, 'null': 7}
```

**Two of nine pairs scored**, both by docking, after the §2 receptor fix
unblocked them:

| Case | Known | Predicted | Delta | Result |
|---|---|---|---|---|
| ALK G1202R / crizotinib | resistance | `likely_reduced_binding` | **+1.763** | **PASS** |
| ALK I1171T / crizotinib | resistance | `uncertain` | +0.027 | uncertain |

That is the first non-zero pass rate this project has produced, and where the
delta cleared the noise floor the direction was right. The threshold was not
touched.

The other seven remain `signal_type="proximity_only"` — a measured distance in
Ångströms with no binding-change prediction attached.

So for most cases SecondLook Tier 2 still answers *"how far is this mutation
from the drug pocket?"* and not *"will this drug still work?"* Those are
different products and the difference must not be blurred in a presentation.

**Resolved (2026-08-23):** the "identical +1.763 for two unrelated cases" was
coincidence — ABL1's delta ranges 1.45–1.86 at exhaustiveness 8. Chasing it
uncovered a genuine reproducibility bug, which is now fixed and verified (§15).
The numbers in this section were nonetheless produced *before* that fix, so
**`validation/results.md` needs regenerating** before the 1/9 is quoted as the
current pass rate.

### Why — and this part is not a bug

Docking's published error is **1.0–2.19 kcal/mol**. The deltas we measured are
**0.02–0.56 kcal/mol**. The instrument is 20–100× less precise than the effect
being measured. That gap does not close with better code or more compute.

Worse, the measured distances show *why* most cases are unreachable:

| Case | Drug | Distance | Band |
|---|---|---|---|
| EGFR T790M | gefitinib | 3.5 Å | `in_contact` |
| EGFR T790M | osimertinib | 3.5 Å | `in_contact` |
| ABL1 T315I | imatinib | 3.5 Å | `in_contact` |
| ALK G1202R | crizotinib | 3.7 Å | `in_contact` |
| ALK I1171T | crizotinib | 6.9 Å | `pocket_adjacent` |
| EGFR C797S | osimertinib | 10.1 Å | `distant` |
| BRAF V600E | vemurafenib | 7.2 Å (chain B) / 12.2 Å (chain A) | `pocket_adjacent` |
| KIT D816V | imatinib | 16.6 Å | `distant` |
| KIT V560G | imatinib | 20.0 Å | `distant` |

Three of nine sit 10–20 Å from the drug. Those are **conformational and
allosteric** mechanisms — the mutation changes the protein's shape or activation
state, not its contact with the drug. A pocket-geometry method cannot see them.
Not with a fix, not ever.

**BRAF V600E is the instructive exception.** Structure 8C7X holds two copies of
the protein and they disagree: residue 600 is 12.2 Å from ligand TXV in chain A
and 7.2 Å in chain B — a 5 Å spread across chemically identical protomers,
straddling the `pocket_adjacent` / `distant` boundary. The pipeline originally
reported whichever chain the *docking* chain-selector happened to return. It now
measures every protomer, reports the closest approach, and says so in the
description when they disagree. The lesson generalises: a single distance from a
single chain of a multimeric structure is an arbitrary choice presented as a
fact.

### What to do about it

Not "try harder." **Narrow the claim to what can be proven, then prove it
completely.** Proximity is a *measurement*, not a prediction. But note what the
BRAF case above proves: reproducible is not the same as correct. Those nine
values reproduced exactly across two full runs and two pipeline versions, and
one of them was measuring the wrong protomer the entire time. Validating a
measurement means re-deriving it independently, not re-running it. See
`docs/proximity-validation-plan.md`.

Critical methodological rule:

- You may **not** move the 70% threshold on the binding-prediction claim. That
  is fitting the test to the data, and `docs/validation-plan.md` forbids it in
  writing.
- You **may** pre-commit a *new* criterion for the *new, narrower* claim —
  written down **before** running it. For example: *"separates contact-mediated
  from distant mechanisms at a 5 Å cutoff across all nine cases."* Current data
  suggests that passes, but write it down first or it doesn't count.

---

## 2. meeko mutant-receptor prep — FIXED via Open Babel fallback

**Resolved 2026-08-22.** Kept in full because the diagnosis took a while and the
trap at the end is worth not falling into.

### The symptom

```
meeko receptor prep failed: Updated 1 H positions but deleted 4
```

Blocked **ABL1 T315I (5HU9)**, **ALK G1202R (4Z55)** and **ALK I1171T (4Z55)** —
three of nine cases, each with a perfectly good ligand-bound structure. The
*wild-type* receptor of the identical structure prepared cleanly, which is what
marked it as a receptor-prep bug rather than a coordinate problem.

### The cause

`meeko/polymer.py::update_H_positions` deletes hydrogens needing repositioning,
re-adds them with `Chem.AddHs`, matches old to new by parent atom, and raises
when the counts differ. A THR->ILE swap changes hydrogen topology at that
residue — THR carries a hydroxyl hydrogen, ILE does not — so the match comes up
short.

Ruled out along the way: not an RDKit sanitisation failure (`SanitizeMol` passes
clean, zero over-valent atoms); not fixable by re-protonating; not fixable by
supplying heavy atoms only (that breaks the wild-type too); not a version issue
(0.7.1 is current).

### The fix

`OpenBabelReceptorPreparer` in `vina_dock.py`. Open Babel does not use residue
templates, so it is unaffected. Verified on the exact structure meeko refuses:

```
5HU9 chain A, T315I
  meeko  -> ReceptorPrepError: Updated 1 H positions but deleted 4
  obabel -> 211,010 byte PDBQT, exit 0
```

It is a **fallback, not the default**. meeko's AutoDock atom typing is the
reference implementation Vina was parameterised against; Open Babel's is a
reimplementation. Preferring it everywhere would trade known-good typing for
merely adequate typing on every structure in order to fix a minority.

Requires `brew install open-babel`. If it is missing, the error says so and
gives the command rather than surfacing a bare `FileNotFoundError`.

### The trap this opened, and how it is closed

Adding a fallback creates a subtler bug than the one it fixes: the wild-type
prepares under meeko, the mutant falls back to Open Babel, and **the delta
between them folds in a preparer difference**. Different atom typing, different
partial charges — and the contamination is the same order of magnitude as the
0.02-0.56 kcal/mol signal being measured. The number would look fine and mean
nothing.

`prepare_receptor_pair()` closes it: the mutant is prepared first (it is the one
that fails), and if either receptor needs the fallback **both** are re-prepared
with it. The pair always travels the same path. Eight tests cover this,
including that the wild-type is never prepared by the primary once the mutant
has refused it.

### What this does not fix

Three more cases can now be *scored*. Whether they produce a **meaningful**
delta is a separate question, and §1 still applies: the measured deltas sit far
below docking's own error bar. This removes a tooling blocker; it does not move
the method's ceiling.

## 3. Labeling calibration is provisional and cannot be calibrated yet

`src/secondlook/labeling.py`:

```python
LABELING_VERSION  = "0.1.0-provisional"
CALIBRATION_STATUS = "provisional"
```

`CALIBRATION_STATUS` rides on every label and every graph node, so nothing can
be mistaken for a validated threshold. Calibration means: get a run that
actually scores, sweep cutoffs with `validation/sweep_thresholds.py` (cheap — it
re-labels cached deltas without re-docking), pick values, flip to `"calibrated"`,
bump `LABELING_VERSION`, record which run produced them.

**Blocked on §2.** Nothing scores, so there is nothing to calibrate against.

---

## 4. One gold-standard case is mislabelled *for this method*

`EGFR T790M / gefitinib` is listed as `resistance` in `GOLD_STANDARD_CASES`.
MdrDB's measured binding change is **−1.21 kcal/mol** — gefitinib binds the
T790M mutant **more tightly**, not less.

Both statements are true. T790M is beyond dispute a clinical gefitinib-resistance
mutation, and its measured effect on gefitinib *binding affinity* is favourable.
Clinical resistance and ΔΔG are different quantities, and here they point
opposite ways.

**Consequence:** on this case a *perfectly accurate* ΔΔG predictor is scored as
wrong. That is a property of the validation set, not of the predictor.

**Do not fix this by flipping the label.** The mechanism usually cited is that
T790M restores ATP affinity, so resistance arises through competition rather
than loss of drug binding — but that explanation has **not been retrieved and
cited**, and under this project's own rule it must be before anyone repeats it.
The defensible action is to mark the case out of scope for a binding delta, with
the measured ΔΔG as the evidence, and say so openly.

Sign convention, established from MdrDB data rather than assumed: `KinaseMD` (a
curated resistance set) is **100% positive**, n=68; TKI 81%; AIMMS 83%. Three
individually known directions agree — ABL1 T315I with dasatinib (**+3.29**) and
nilotinib (**+3.00**), both clinically resistant, versus axitinib (**−1.26**),
clinically active against T315I. **Positive ΔΔG = weaker binding = resistance.**

---

## 5. No clinical-trial source — half the activation rule is unevaluated

`tier1-retrieval.md` §Activation defines a strong hit as *"≥1 CIViC evidence
item at level A or B for this exact variant, **OR** ≥1 recruiting trial matching
this exact biomarker."*

Tier 1 loads CIViC and PubMed. It loads **no** trial source — no
ClinicalTrials.gov, no CTRI. So the second half of that clause cannot be
evaluated.

Treating "no trial data" as "no trial" would silently convert an unimplemented
source into a negative clinical finding. `tier1_adapter._classify` therefore
appends to every strong/weak rationale:

> *Trial half of the strong-hit rule not evaluated (no trial source loaded).*

**Fix:** implement a ClinicalTrials.gov loader in Tier 1 (endpoints and rate
limits are in `docs/data-sources.md`), then drop the caveat string. CTRI access
mode still needs verifying and is tracked as Tier 1 work.

---

## 6. All nine gold-standard cases return `no_hit` from Tier 1 — this is correct

Tier 1's loaded scope is 17 rare cancers (`civic_scope.yaml`): NTRK fusions,
EWSR1::FLI1, PAX3::FOXO1 and similar. **None of EGFR / ABL1 / KIT / BRAF / ALK
are in it.**

Those nine cases were chosen because they have well-established ground truth for
validating the *structural method* — not because they are the target population.
So on the real graph every one returns `no_hit` and hands off to Tier 2.

**That is the two-tier architecture working exactly as designed**, and it is
worth saying plainly rather than letting someone discover it and read it as a
gap. Verified live alongside three in-scope variants (FGFR1 N546K, TP53 R273H,
KDR A1065T) which correctly return `weak_hit` with real CIViC evidence attached.

---

## 7. Repo layout — the two-repo split is gone

**Resolved by this repo.** Recorded because the failure mode is instructive and
the old repos still exist.

Athena was originally two repos — `Athena_tier_2` and `Athena_Tier_1` — sharing
one `secondlook` package as a PEP 420 implicit namespace. That arrangement had
three fragile rules (neither repo may ship `src/secondlook/__init__.py`; both
must install with `--config-settings editable_mode=compat`; they must be
siblings on disk), and a cross-repo git pin that had drifted to the wrong fork
at a stale commit.

**This repo is one package.** `secondlook` with `secondlook.tier1` as a normal
subpackage. One `pip install -e .`, one `pytest`, no namespace rules, no pin.
Every failure mode above is structurally impossible here.

### The old repos

| Repo | Status |
|---|---|
| `KrishG7/Athena_tier_2` | Superseded by this repo. History preserved. |
| `KrishG7/Athena_Tier_1` | Fork of Gagan's, with the §8a schema fix. Superseded. |
| `Gagansharma-code/Athena_Tier_1` | **Still has the §8a bug.** Gagan should take the fix. |

If you have an old clone, switch to this repo. If you're Gagan: the
`StructuralSignal` schema fix in §8a should land in your repo too, or anyone
cloning from you gets a schema that silently drops the only field Tier 2
currently outputs.

## 8. Two integration bugs found only by running against real data

Both passed every unit test with a mock. Recording them because the lesson
generalises: **the seam is where the bugs are.**

### (a) StructuralSignal schema drift — fixed

Tier 2 emitted 13 node properties; Tier 1's `STRUCTURAL_SIGNAL` NodeType
declared 7. Nothing crashed — `NodeType` is documentation, not an ORM — so a
query built from the declaration would simply never have read the missing six,
one of which is `binding_site_distance_angstrom`, *the entire output* of every
current gold-standard case.

Fixed in Tier 1's schema, plus a test there comparing the declaration against
what Tier 2 actually emits, so it cannot drift silently again.

### (b) The `hgvs_p` variant-join mismatch — fixed

CIViC writes `NP_005148.2:p.Thr315Ile`; Tier 2 computes `p.Thr315Ile`. The graph
sink's original `MERGE (v:Variant {hgvs_p: ...})` on the bare form matched **0 of
1** existing TP53 R273H nodes on the live graph. It would have minted a second,
orphan Variant for a variant the graph already had — unreachable from its Gene
and invisible to every Tier 1 traversal.

The sink now resolves find-then-create against
`(Gene)-[:HAS_VARIANT]->(Variant)` with an accession-anchored suffix match.
Verified: 9 signals → 9 full `Gene→Variant→Signal→Drug` paths, 0 orphans, EGFR
T790M correctly resolving to one Variant carrying two drug signals.

---

## 9. Environment issues that will bite you

### Docker wedges when the disk fills

During development the disk hit **263 MiB free of 228 GB**. `sentence-transformers`
(torch, ~2 GB) failed mid-install, and Docker Desktop's backend hung — `docker ps`
stopped responding entirely and would not restart from the CLI.

**Keep ≥3 GB free.** Reclaim with `docker builder prune` (safe) and by removing
stale images (`docker images` — check what you actually still need first).

### The test suite must skip, not error, when FalkorDB is down

`tests/test_tier1_integration.py`'s skip guard originally caught
`falkordb.exceptions.ConnectionError` — an attribute that **does not exist**;
`falkordb.exceptions` defines no error classes, the client raises redis-py's
types. And `redis.exceptions.ConnectionError` does **not** subclass the builtin
`ConnectionError`/`OSError`, so the fallback missed it too.

Result: with the database down the suite **errored** instead of skipping. Fixed
to catch `redis.exceptions.ConnectionError`/`TimeoutError`. If you add live
tests, get this right — a suite that errors when a service is down looks broken
in a demo.

### Live services do go down

An Ensembl VEP outage and a UniProt HTTP 500 both produced failures during
development that looked like regressions and were not. `http_retry.py` handles
429/5xx with backoff. Before diagnosing a failure as a code problem, check the
service:

```bash
curl -o /dev/null -w "%{http_code}\n" https://rest.uniprot.org/uniprotkb/P00533.fasta
curl -o /dev/null -w "%{http_code}\n" https://files.rcsb.org/download/2ITY.pdb
```

---

## 10. Smaller open items

| Item | Status |
|---|---|
| SIFTS UniProt↔PDB residue mapping | **Implemented** — `src/secondlook/sifts.py`, 14 tests. See below. |
| Demo cases locked and cached | Cache layer exists (`cache.py`); **no case is actually locked**. Needs a run that scores first. |
| Vina-vs-mCSM-lig coverage note | **Written** — `docs/method-coverage-note.md`. Measured split: 100% docking, 0% mCSM-lig. |
| Mode 3 semantic retrieval | Never verified end to end — 7 Tier 1 integration tests blocked on `sentence-transformers`, blocked on disk. |
| MdrDB full set | Only the 4,292-row CoreSet is in hand. The 100,537-row full set is behind a registration form needing an academic email. |

---

## 11. MdrDB licence — read before shipping anything

Retrieved from <https://quantum.tencent.com/mdrdb/download> (the repo has no
LICENSE file; terms live on the download page):

- **"Academic, or Educational Use Only."**
- **"a non-exclusive, non-transferable, limited license for academic research."**
- **"Not use the database or its contents for commercial purposes unless
  explicitly approved in writing by us"** — commercial licensing via
  `qlab@tencent.com`.
- Users must correctly cite the data source.

| Action | Allowed |
|---|---|
| Download and use as a university student | **Yes** |
| Train and evaluate a model for this project | **Yes** |
| Report metrics computed from it, citing MdrDB | **Yes** |
| Commit the TSV into this repo | **No** — the licence is *non-transferable* |
| Ship SecondLook commercially with MdrDB-derived data or weights | **Not without written approval** |

The paper itself is CC BY 4.0 and open access — DOI `10.1038/s42004-023-00920-7`,
*Communications Chemistry* 6:123, 2023. No university network needed.

**If SecondLook is pitched as a product, this needs a commercial licence.** SIH
submissions routinely gesture at deployment — decide early which track you are
in. One email settles it, and it is far cheaper to send now than after a model
is trained.

---

## 12. Suggested order of work

1. **Commit and push.** Everything else is worthless if the work is lost.
2. **Environment** — free disk, restart Docker, install `sentence-transformers`,
   verify Mode 3. Fix the cross-repo pin (§7).
3. **Pre-commit the proximity criterion**, then evaluate it. This is the claim
   you can actually make bulletproof.
4. **Timebox the meeko fix** (§2) — half a day via Open Babel, not three days.
   It is the only remaining work with real scientific payoff.
5. **MdrDB training run** to reproduce the published baseline (Scenario 3,
   ExtraTrees: RMSE 0.656, Pearson 0.607, R² 0.362). Gives you an accuracy
   figure you own, and it is the honest roadmap slide.

Note on that baseline before treating it as solved: in the paper's Scenario 1,
**all ten models have negative R²** — worse than predicting the mean. Scenario 3
is the only configuration that clears zero, and even there RMSE 0.656 kcal/mol
is the same order of magnitude as docking's error bar.

---

## 13. Proximity measurement — independently verified, one real error found

Criteria pre-registered in `docs/proximity-validation-plan.md` **before** any
evaluation, then evaluated with `validation/verify_proximity.py`, which
re-derives every distance from raw PDB coordinates and deliberately does not
import `proximity.py`.

| Criterion | Result |
|---|---|
| C1 geometric accuracy (≤ 0.05 Å) | **FAIL** — 8 of 9 |
| C2 residue identity | PASS — 9 of 9 |
| C3 ligand identity stated and present | PASS — 9 of 9 |
| C4 reproducibility | PASS |

### What failed

**BRAF V600E**: pipeline reported **12.2452 Å**, independent recomputation gave
**7.2077 Å**. That is 100x the tolerance, and it crosses a band boundary —
`distant` becomes `pocket_adjacent`.

Not an arithmetic error. 8C7X is a homodimer and the two protomers genuinely
differ: 12.2452 Å in chain A, 7.2077 Å in chain B. Both correct. The pipeline
measured in whichever chain the *docking* chain-selector returned — a rule
written for a task needing exactly one receptor — and reported it as the
structure's answer, with nothing recording that the other copy disagreed by 5 Å.

The other eight were exact to four decimal places.

### Fixed

`measure_proximity` now measures every chain holding both the residue and a
ligand, reports the closest approach, keeps all per-chain values in
`ProximitySignal.protomer_distances`, and states the disagreement in the
description when protomers land in different bands. Seven tests cover it.

### Re-verified — all four criteria now pass

Re-run on 2026-08-23 against PDB files re-downloaded from RCSB (the cache was
deleted first, so this is not a check against the same artifacts):

| Criterion | Result |
|---|---|
| C1 geometric accuracy (<= 0.05 Å) | **PASS** — 9/9 |
| C2 residue identity | PASS — 9/9 |
| C3 ligand identity stated and present | PASS — 9/9 |
| C4 cross-run reproducibility (<= 0.01 Å) | **PASS** — 9/9 |

BRAF V600E now reports 7.2077 Å, matching independent recomputation exactly. All
nine agree to four decimal places. `verify_proximity.py` exits 0.

**C4 had never actually run.** `C4_TOLERANCE_ANGSTROM` was defined but unused, so
the first evaluation's "C4 PASS" had no check behind it. The verifier now records
a baseline of independently recomputed distances in
`validation/proximity-reference.json` and compares against it on every later run,
and C4 counts toward the exit status. The 9/9 above is the first real C4 result.

This validates the proximity measurement only — not the docking deltas, which
have their own problem (§15).

### Why this matters beyond one number

Those nine distances reproduced **exactly** across two full pipeline runs at two
different pipeline versions. One of them was measuring the wrong protomer the
entire time.

Reproducible is not the same as correct. Any check built on re-running and
comparing would have confirmed the wrong number every time. That is the argument
for independent re-derivation, and it is the reason C1 was written the way it
was.

---

## 14. SIFTS residue mapping — implemented, and it confirms the assumption

The pipeline numbers residues by **UniProt** position (`T315I` = position 315 in
P00519) but indexes into **PDB** files, which number residues however the
depositor chose. The two agree often enough to be dangerous and not often enough
to assume.

When they disagree, nothing fails. The pipeline measures from a real residue at
a real position that is simply the wrong one, and reports a confident distance
for it. The existing residue-identity check catches most of that — but not a
THR-to-THR collision, which passes silently.

`src/secondlook/sifts.py` reads EBI PDBe's curated mapping so the assumption is
checked rather than trusted.

### Result on our structures

| Structure | UniProt | Position | Crystallised | Author numbering |
|---|---|---|---|---|
| 5HU9 | P00519 | 315 | yes | 315, offset +0 |
| 8C7X | P15056 | 600 | yes | 600, offset +0 |
| 8PQD | P10721 | 816 | yes | 816, offset +0 |
| 8S14 | P10721 | 560 | yes | 560, offset +0 |
| 8A27 | P00533 | 790, 797 | yes | not supplied by PDBe |
| 4Z55 | Q9UM73 | 1202, 1171 | yes | not supplied by PDBe |

**Every residue is genuinely present, and every entry that supplies numbering
confirms offset +0.** The assumption held — but it held by luck until now, and
this is what turns that into a checked fact.

### The one design decision worth knowing

PDBe does not always populate `author_residue_number`; EGFR 8A27 and ALK 4Z55
return null for both segment ends. The module keeps those segments rather than
discarding them, because *"present but I cannot tell you its number"* and
*"absent from this structure"* are different facts and only the second should
stop a structure being used. `is_crystallised()` answers the first,
`map_position()` the second, and `numbering_agrees()` returns `None` rather than
`False` when it cannot check.

Where only `author_end` is populated (BRAF 8C7X), the start is derived
backwards. The segment is contiguous, so that is exact, not a guess.

### Not yet wired into the pipeline

The module and its tests exist and are verified live. `structure.py` and
`vina_dock.py` still use the UniProt number directly. Wiring it in is the next
step — the natural place is a check at structure selection that refuses a
structure whose numbering SIFTS says disagrees, rather than silently measuring
the wrong residue.

---

## 15. Scored deltas were not reproducible at a fixed seed — cause found and fixed

**This supersedes the "identical delta" question, which turned out to be
coincidence. The reproducibility problem it uncovered was real, and has now been
traced to a specific cause and fixed.**

### Measured

ABL1 T315I / imatinib, 5HU9, exhaustiveness 8:

| Run | Wild-type | Mutant | Delta |
|---|---|---|---|
| seed 1 (first) | — | — | **+1.763** |
| seed 1 (repeat) | -11.440 | -9.987 | **+1.453** |
| seed 7 | -11.615 | -9.755 | **+1.860** |

Two observations, and the second is the serious one:

1. **The coincidence is resolved.** ABL1's delta is not pinned at +1.763; it
   ranges roughly 1.45-1.86. Matching ALK G1202R's +1.763 exactly once was
   chance. No bug.

2. **The same seed did not reproduce the same delta.** +1.763 and +1.453 came
   from identical settings — same structure, ligand, seed, exhaustiveness. The
   run-to-run spread was about **0.4 kcal/mol**, the same order as the effects
   being measured.

The component scores were healthy (-11.4 to -11.6 wild-type, -9.8 to -10.0
mutant), so these were real docked poses, not failed runs. `BindingScore`
carrying `.wildtype_score` / `.mutant_score` is what made this measurable.

### Diagnosis

Each preparation stage was run twice on identical input and the outputs compared
byte for byte. That isolated the variation immediately:

| Stage | Reproducible? |
|---|---|
| Grid box derivation | yes |
| Meeko ligand prep (SMILES → PDBQT) | yes |
| Receptor prep given identical PDB input (meeko / Open Babel) | yes |
| Vina itself (already `cpu=1`, fixed seed) | yes |
| **PDBFixer protonation and mutant side-chain placement** | **no** |

So Vina was never the problem, and neither was charge assignment. **PDBFixer was
returning different coordinates on every call for the same input** — and not only
for hydrogens: backbone heavy atoms moved by up to ~0.07 Å. Those coordinates
become the docking receptor, so the drift reached the score.

Three independent sources of randomness, all of which had to be pinned:

1. **`addMissingAtoms`** minimizes newly-placed heavy atoms with a Langevin
   integrator seeded from OS entropy unless a `seed` argument is passed.
2. **`addMissingHydrogens`** (OpenMM's `Modeller`) places every added hydrogen at
   a `random.random()` offset and then minimizes — it consumes the process-wide
   Python RNG. This was the dominant term.
3. **OpenMM's CPU platform is multi-threaded.** Forces are summed in whatever
   order threads finish, and floating-point addition is not associative, so each
   minimization takes a slightly different path.

Seeding alone does not work — with the threading and `addMissingAtoms` sources
still live, a seeded run still diverges. All three fixes are required together.

### Fix

In `src/secondlook/vina_dock.py`:

- `_pdbfixer_complete` now passes `seed=DEFAULT_SEED` to `addMissingAtoms`,
  seeds Python's RNG around `addMissingHydrogens` (saving and restoring the
  caller's RNG state, so structure prep does not silently reset anyone else's
  random stream), and runs on a single-threaded CPU platform.
- `minimize_mutated_sidechain` sets its integrator's random seed and uses the
  same single-threaded platform.
- New `deterministic_cpu_platform()` returns the pinned platform, or `None` when
  OpenMM is absent so callers degrade instead of failing.

### Verified

Every preparation stage, run three times through the patched code path on 5HU9
chain A — all byte-identical, where they differed on every run before:

```
stage            run1             run2             run3             verdict
protonate        4c85a6e69d7fa8ff 4c85a6e69d7fa8ff 4c85a6e69d7fa8ff IDENTICAL
mutant           f30597f60d5f09a9 f30597f60d5f09a9 f30597f60d5f09a9 IDENTICAL
rec_wt           ffb4d19045206e95 ffb4d19045206e95 ffb4d19045206e95 IDENTICAL
rec_mut          71f54a8fc9e0e0c8 71f54a8fc9e0e0c8 71f54a8fc9e0e0c8 IDENTICAL
ligand           d879e65ff5439f15 d879e65ff5439f15 d879e65ff5439f15 IDENTICAL
```

And Vina itself, three runs at **production exhaustiveness 32** with the
prepared files held constant:

```
run           wildtype    mutant     delta
seed 1 #1      -11.622    -9.697     1.925
seed 1 #2      -11.622    -9.697     1.925
seed 1 #3      -11.622    -9.697     1.925
```

Bit-identical. That closes the question definitively: **Vina was never the
source of the variance** — with `cpu=1` and a fixed seed it reproduces exactly,
even at exhaustiveness 32. The whole ~0.4 kcal/mol spread came from structure
preparation, and structure preparation is now deterministic.

Prep and search are both reproducible, so a re-run of the gold standard at a
fixed seed should now return the same deltas every time.

Regression tests are in `tests/test_vina_dock.py` — two unit tests (thread
pinning, RNG-state restoration) plus an `integration`-marked test asserting
byte-identical prep on a real structure.

### Still outstanding

Reproducibility is settled. What is not yet done:

- **Re-run the gold standard** with the fix in place. `validation/results.md`
  was produced before it, so its numbers came off the non-deterministic prep.
  This is the one remaining item.
### Measured noise floor (exhaustiveness 32, prep held constant)

> **Provenance note.** These eight runs reused a receptor pair prepared *before*
> the determinism fix landed. Prep was held constant across all eight, so they
> measure exactly what they claim — Vina's sensitivity to its seed — but the
> absolute delta is specific to that conformer. Re-timed later against the
> patched, canonical prep, the same case gives -11.622 / -9.915 = **+1.707**.
> The spread figures below are unaffected; only the absolute value moves.


ABL1 T315I / imatinib, 5HU9, eight runs — three at seed 1, then seeds 2–6:

```
same seed, 3 runs : 1.9250, 1.9250, 1.9250   -> DETERMINISTIC
all 8 runs        : mean 1.947  sd 0.048  min 1.887  max 2.058  range 0.171
```

Two separate numbers, and they must not be confused:

- **Run-to-run precision at a fixed seed: exactly zero.** Production pins the
  seed, so this is the figure that governs whether a reported delta repeats.
- **Seed-to-seed spread: sd 0.048, range 0.171 kcal/mol.** A measure of how
  converged Vina's search is at exhaustiveness 32. Much tighter than the ~0.4
  seen at exhaustiveness 8, which is what made the earlier probes look noisy.

**Precision is not accuracy.** These figures say a delta *repeats*, not that it
is *right*. Docking's published error against measured binding affinity is still
1.0–2.19 kcal/mol (§1), which is 20–100× the effects in the validation set. A
number can be reproducible to four decimal places and still be wrong by 2
kcal/mol. §1's argument is unaffected by anything in this section.

Until that is done, **`validation/results.md` should still be read as one
sample, not a result**, and the single PASS should not be quoted as a finding.
The numbers in it were produced before this fix.

---

## 16. Tier 1 variant-key routing downgrades real strong hits

**Found while building `ask.py`, then corrected against the live graph — the
first version of this section was based on a wrong assumption about how CIViC
stores variants. What follows is what the loaded graph actually contains,
queried on 2026-08-23.**

### The data has three conventions, not one

```
gene          v.name            v.hgvs_p
ABL1          (empty)           p.Thr315Ile
ALK           (empty)           p.Ile1171Thr
BRAF          (empty)           p.Val600Glu
ERCC2         K751Q             (empty)
FGFR1         N546K             NP_075598.2:p.Asn546Lys
ETV6::NTRK3   Fusion            (empty)
CDK4          Amplification     (empty)
```

`name` is populated for essentially every variant that carries evidence.
`hgvs_p` is populated for only a handful, and when it is, CIViC writes it
RefSeq-prefixed — so exact equality against a bare `p.Asn546Lys` fails too.

### Correction: the ABL1/ALK/BRAF/EGFR rows are not CIViC data

Those four have an empty `name` and a bare, unprefixed `hgvs_p`. They were
**created by Tier 2's own graph sink**, not loaded from CIViC. Confirmed:

```
MATCH (g:Gene {symbol:'ABL1'})-[:HAS_VARIANT]->(v:Variant) ...
p.Thr315Ile | (no name) | HAS_COMPUTATIONAL_SIGNAL -> StructuralSignal | 0 evidence items
```

Zero evidence items, and the only edge is to a `StructuralSignal`. So Tier 1
reporting `no_hit` for ABL1 T315I is **correct** — imatinib resistance in CML is
outside the loaded 17-disease scope. That is a coverage limit (§5), not a bug,
and the earlier claim that a key change would fix it was wrong.

### The actual bug, demonstrated live

`run_tier2` passes its `mutation` argument through as `hgvs_p`. For any variant
whose identity lives in `name` — every fusion, amplification, expression call,
and most missense — the exact-match query finds nothing and Tier 1 falls through
to Mode 2's relaxed gene-level query. It still returns items, so nothing looks
broken. But the verdict is downgraded:

```
tried mutation='Fusion'      -> weak_hit   (4 items)   "No exact-variant evidence"
tried variant_name='Fusion'  -> strong_hit (4 items)   "1 CIViC evidence item at level A"
```

Same variant, same four items, opposite decisions. And the decision is the one
that matters: `strong_hit` means *published Level A evidence answers this, do
not run Tier 2*; `weak_hit` means *run Tier 2*. So the bug both misreports
evidence strength and spends a full structural run on a question already
answered by a phase 1/2 trial.

### The workaround in place

`ask.py` routes on the shape of the string:

| Input looks like | Key used |
|---|---|
| contains `p.` (`p.Thr315Ile`, `NP_005148.2:p.Thr315Ile`) | `hgvs_p` |
| shorthand (`T315I`) | `variant_name`, then `p.Thr315Ile` as `hgvs_p` |
| anything else (`Fusion`, `Amplification`, `Loss`) | `variant_name` |

`Tier1RetrievalPolicy.decide` already accepts both keywords, so the CLI needed
no library change.

### The real fix

Move that routing into `run_tier2` before it calls `activation_policy.decide`,
so every caller benefits. `ask.py`'s `tier1_lookups()` is the reference
implementation and its tests (`tests/test_ask_cli.py`) are the reference cases.

Worth considering alongside it: `retrieve_exact` could match `hgvs_p` with
`ENDS WITH` as well as equality, the way `FalkorDBGraphSink` already does for
exactly this RefSeq-prefix reason. That would make `p.Asn546Lys` match
`NP_075598.2:p.Asn546Lys` instead of silently relaxing to Mode 2.

---

## 17. Docking failures are swallowed, so nobody can tell why 7 of 9 declined

**This is a diagnostic gap, and it may be hiding a plain configuration bug.
Read this before accepting §1's "docking cannot answer these cases" as settled.**

### The swallow

`binding._vina_or_unavailable` ends:

```python
    except vina_error:
        return _unavailable(BINDING_UNAVAILABLE_MESSAGE)
```

`vina_error` is `VinaError`, the base class for **every** docking failure:
`ReceptorPrepError`, `LigandPrepError`, `MutantPlacementError`,
`NoBindingSiteError`, `VinaTimeoutError`, `VinaRunError`. All six collapse into
one generic sentence:

> "No structure-based binding-affinity signal could be computed for this
> candidate (mCSM-lig and docking were both unavailable or inapplicable)."

The exception carries the actual cause; it is discarded rather than recorded.
`MutationOutsidePocketError` is handled properly right above it — with its own
message and a `reason_code="outside_pocket"` — which shows the intended pattern.
Every other failure mode just vanishes.

This is the same class of problem §13 and §14 were written to prevent: an output
that cannot distinguish between causes that mean very different things. A
receptor-prep failure is a bug to fix. A timeout is a setting to raise. "Docking
is inapplicable to this chemistry" is a scientific limit. They are
indistinguishable in the current output, and §1's headline rests on assuming the
third.

### Why this matters right now

Run through the full pipeline, **ABL1 T315I / imatinib returns
`proximity_only`** — "docking unavailable or inapplicable" — on structure
**5HU9**.

Docking that exact structure directly works fine, and scores it:

```
wildtype -11.622   mutant -9.697   delta +1.925    (exhaustiveness 32, seed 1)
```

So for at least one case, the "docking could not answer" verdict is false.
Docking answers it in about six minutes when called directly.

### The cause — measured, not inferred

`DEFAULT_TIMEOUT_SECONDS = 120.0`, applied **per docking run** (so a
wild-type/mutant pair costs two of them). Timing a single run directly:

```
pipeline default timeout: 120.0s
rec_wt   score  -11.622     377.3s     3.1x the timeout   WOULD TIME OUT
rec_mut  score   -9.915     371.2s     3.1x the timeout   WOULD TIME OUT
```

**Both runs ~375 seconds against a 120-second limit.** Every docking attempt on
this structure raised `VinaTimeoutError`, was swallowed by the clause above, and
came back as `proximity_only` with a message saying docking was "unavailable or
inapplicable".

That is the 7-of-9 decline rate explained as a configuration bug, not a limit of
the method. The score it produces when allowed to finish (-11.622) is the same
value to three decimals as the standalone runs in §15, so the calculation was
always working — it was just being killed.

The two cases that *did* score (both ALK/crizotinib on 4Z55) were the ones whose
search happened to finish inside 120 s — a smaller box or fewer rotatable bonds,
not a scientific difference. That also means the surviving deltas were the
*fastest* dockings rather than a representative sample, which is its own
selection effect on the reported pass rate.

### Fixed

1. **The swallow is gone.** Each `VinaError` subclass now carries its own message
   and `reason_code`, and the underlying exception detail is interpolated into
   the message rather than discarded:

   | Failure | `reason_code` |
   |---|---|
   | `VinaTimeoutError` | `docking_timeout` |
   | `ReceptorPrepError` | `receptor_prep_failed` |
   | `LigandPrepError` | `ligand_prep_failed` |
   | `MutantPlacementError` | `mutant_placement_failed` |
   | `NoBindingSiteError` | `no_binding_site` |
   | anything else deriving from `VinaError` | `docking_failed` |

   The timeout message says explicitly that the calculation was *cut off, not
   refused* — the distinction that was lost before.

2. **`DEFAULT_TIMEOUT_SECONDS` raised from 120 s to 900 s**, with the measurement
   and reasoning recorded at the constant. It is a runaway guard, not a budget.

   Two tests that asserted the old generic message were updated rather than
   deleted — they had encoded the swallow as intended behaviour. A new test
   asserts a timeout is reported as a limit and not as a refusal.

### Still to do

3. **Re-run the gold standard** and see what the pass rate actually is now that
   docking is allowed to finish. This is the open item, and until it is done the
   1-of-9 in §1 describes a run in which docking was being killed at 120 s.

### Related, unverified: is mCSM-lig's 0% coverage also a timeout?

`docs/method-coverage-note.md` records the measured split as **100% docking, 0%
mCSM-lig** and (correctly) warns against quoting mCSM-lig's published rho=0.67
near Vina results. But `mcsm_lig.DEFAULT_TIMEOUT_MS = 180_000` — three minutes
for a Playwright-driven web submission that queues server-side — was chosen the
same way the 120 s docking timeout was: not measured.

An mCSM-lig timeout is caught and falls through to Vina, which is correct
behaviour and means no result is lost. But it makes "mCSM-lig scored 0 of 9"
ambiguous in exactly the way §17 is about: it could mean the tool declined
these cases, or it could mean the scrape never finished. Those support very
different claims, and the coverage note currently assumes the first.

**Worth one measurement before that split is quoted anywhere:** time a single
mCSM-lig submission end to end and compare it with the 180 s limit. If it
exceeds, raise it and re-measure the split. Cheap to check, and it decides
whether a documented deliverable is reporting a property of the tool or a
property of our timeout.

Until 1–3 are done, **§1's "the instrument is 20–100x less precise than the
effect" argument stands on its own merits, but the 7-of-9 decline rate is not
evidence for it.** Those are separate claims and only the first is currently
supported.

---

## 18. The protomer fix reached the proximity signal but not the docking gate

**Found by reading the reason codes §17 added, on the gold-standard re-run.
Fixed. Affected cases need re-running.**

### The contradiction

§13 fixed `measure_proximity` to report the closest approach across protomers,
because 8C7X is a BRAF homodimer whose copies disagree: V600 sits **12.2 Å**
from the ligand in chain A and **7.2 Å** in chain B.

The docking path was never given the same fix. `select_docking_chain` picked a
chain by **ligand size**, which on a homodimer is effectively arbitrary — it
chose chain A. The 8 Å contact gate then measured 12.2 Å and refused:

```
reason_code = outside_pocket
reason      = This mutation lies 12.2 A from the drug-binding site, beyond the
              8 A contact range within which a docking score can detect a
              binding change.
```

While, in the same result object, the proximity signal read:

```
distance_angstrom = 7.207651281797702
description = "...This structure holds 2 copies of the protein and they do not
               agree (chain A 12.2 A, chain B 7.2 A); the closest approach is
               reported..."
```

**One result, two contradictory distances, and the wrong one drove the
decision.** 7.2 Å is inside the 8 Å gate, so BRAF V600E should have been docked.
It was refused for the entire gold-standard run.

### Fix

`select_docking_chain` now ranks candidate chains by **how close the mutated
residue is to that chain's ligand**, falling back to ligand size as a tie-break
so the original buffer-component reasoning still applies. That is also the right
question on its own terms: the protomer where the mutation can actually
influence binding is the one docking is being asked about.

A new `min_residue_ligand_distance()` gives the same closest-approach number to
the refusal message, so the distance quoted always matches the one the decision
used.

Verified against the real structure:

```
BRAF V600E / 8C7X
  chain selected for docking : B
  distance in that chain     : 7.2077 A
  contact gate               : 8.0 A
  -> WILL DOCK
  agrees with proximity 7.2077: True
```

Four tests cover it, including one asserting the invariant the bug violated:
the distance the gate uses must equal the distance the proximity signal reports.

### Consequence for the current re-run

The gold-standard run in flight imported the pre-fix code, so its **BRAF V600E
result is stale** — cached as `outside_pocket` from chain A. Delete that cache
entry and re-run the case.

Any homodimer whose protomers disagree across the 8 Å boundary is affected. The
other `outside_pocket` refusals in the same run (KIT D816V at 16.6 Å, KIT V560G
at 20.0 Å) are far enough out that no protomer disagreement could bring them
inside the gate, so those stand.

### The general lesson, again

This is the third time a fix has been applied to one code path and not to its
twin (§8's schema drift, §16's variant keys, now this). The proximity and
docking paths both answer "how far is this residue from the drug?" and had
drifted to different answers. The regression test now pins them together rather
than trusting that a future fix reaches both.

---

## 19. A None key crashed the gold-standard runner on EGFR T790M

**Found because the run kept dying at case 1 and I assumed memory pressure. It
was a plain bug, reproducible every time. Fixed.**

### The crash

```
TypeError: '<' not supported between instances of 'str' and 'NoneType'
  File ".../json/encoder.py", line 354, in _iterencode_dict
    items = sorted(dct.items())
```

`Tier2Output.scored_by_method` is annotated `dict[str, int]` but was built as:

```python
counts[item.method] = counts.get(item.method, 0) + 1
```

`Tier2ResultItem.method` is `ScoringMethod | None` — None on every result that
produced no delta. So the dict could hold a None key, and `cache.save_payload`
serialises with `sort_keys=True`, which calls `sorted()` over the keys.

### Why only EGFR T790M

`sorted()` succeeds when the keys are *all* strings or *all* None. It raises only
on a mix. A payload holds one entry per drug, so the mix requires a single case
with one scored and one unscored drug.

**EGFR T790M is the only gold-standard case with two drugs** — gefitinib
(docked, `method="docking"`) and osimertinib (covalent, refused,
`method=None`). Exactly the mixed-key combination, every run.

### What it cost

The runner is ordered with EGFR T790M first, so the crash killed the process
before any other case ran. Three separate "the run was killed" incidents were
this, not the memory pressure I attributed them to. The run only ever got past
case 1 while a stale `EGFR_T790M.json` existed and the case was skipped — which
is why clearing that cache entry, meant to force a clean re-run, made the
failure permanent instead.

Worth noting how it hid: the caching layer is designed so a failed re-run cannot
destroy a good cached result, and that worked. But it also meant the only case
that could trigger the bug was usually skipped.

### Fix

Unscored results now count under the string key `"null"`:

```python
key = item.method if item.method is not None else "null"
```

`"null"` rather than `"unscored"` because `json.dumps` already rendered the None
key that way — every payload written before this fix uses it, so readers do not
have to handle two spellings. The annotation `dict[str, int]` is now true.

A regression test asserts no key is ever None and that
`json.dumps(..., sort_keys=True)` succeeds on a mixed scored/unscored payload.

### Lesson

The type annotation said `dict[str, int]` and was wrong; nothing checked it. The
same class of mismatch as §16's `dict[str, int]`-shaped assumptions about the
graph. A type checker over `pipeline.py` would have caught this before it cost
three runs — worth considering as a CI step.
