import pytest

from secondlook.binding import BindingScore
from secondlook.vina_dock import (
    LigandPrepError,
    MutantPlacementError,
    NoBindingSiteError,
    ReceptorPrepError,
    VinaDockClient,
    VinaError,
    VinaRunError,
    VinaTimeoutError,
    grid_box_from_pdb,
    ligand_hetatm_coords,
    parse_mutation_shorthand,
    protein_atoms_only,
)

LIGAND_PDB = """\
ATOM      1  CA  ASP A  30      1.000   0.000   0.000  1.00 20.00           C
HETATM    2  C1  LIG A  99      0.000   0.000   0.000  1.00 20.00           C
HETATM    3  C2  LIG A  99     10.000   0.000   0.000  1.00 20.00           C
HETATM    4  O   HOH A 100     99.000  99.000  99.000  1.00 20.00           O
END
"""

APO_PDB = """\
ATOM      1  CA  ASP A  30      1.000   0.000   0.000  1.00 20.00           C
END
"""

PEPTIDE_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       2.009  -0.771  -1.214  1.00 20.00           C
ATOM      6  N   ASP A   2       3.458   1.654   0.000  1.00 20.00           N
ATOM      7  CA  ASP A   2       4.200   2.900   0.000  1.00 20.00           C
ATOM      8  C   ASP A   2       5.700   2.700   0.000  1.00 20.00           C
ATOM      9  O   ASP A   2       6.200   1.600   0.000  1.00 20.00           O
ATOM     10  CB  ASP A   2       3.700   3.700   1.200  1.00 20.00           C
ATOM     11  CG  ASP A   2       2.200   4.000   1.200  1.00 20.00           C
ATOM     12  OD1 ASP A   2       1.400   3.200   0.600  1.00 20.00           O
ATOM     13  OD2 ASP A   2       1.900   5.100   1.700  1.00 20.00           O
ATOM     14  N   ALA A   3       6.400   3.800   0.000  1.00 20.00           N
ATOM     15  CA  ALA A   3       7.850   3.800   0.000  1.00 20.00           C
ATOM     16  C   ALA A   3       8.400   5.200   0.000  1.00 20.00           C
ATOM     17  O   ALA A   3       7.650   6.200   0.000  1.00 20.00           O
ATOM     18  CB  ALA A   3       8.400   3.000  -1.200  1.00 20.00           C
HETATM   19  C1  LIG A  99      4.200   2.900   5.000  1.00 20.00           C
END
"""


class FakePrep:
    def __init__(self, error: Exception | None = None) -> None:
        self.error = error
        self.calls: list[str] = []

    def to_pdbqt(self, pdb_text: str) -> str:
        self.calls.append(pdb_text)
        if self.error:
            raise self.error
        return f"PDBQT:{hash(pdb_text)}"


class FakeMutant:
    def __init__(self, error: Exception | None = None, pdb: str = "MUTANT") -> None:
        self.error = error
        self.pdb = pdb
        self.calls: list[tuple[str, str, int]] = []

    def place_sidechain(self, pdb_text: str, mutation: str, position: int) -> str:
        self.calls.append((pdb_text, mutation, position))
        if self.error:
            raise self.error
        return self.pdb


class FakeLigand:
    def __init__(self, error: Exception | None = None) -> None:
        self.error = error
        self.calls: list[str] = []

    def to_pdbqt(self, smiles: str) -> str:
        self.calls.append(smiles)
        if self.error:
            raise self.error
        return "LIGAND_PDBQT"


class FakeDock:
    def __init__(self, scores: list[float] | None = None, error: Exception | None = None) -> None:
        self.scores = list(scores or [])
        self.error = error
        self.calls: list[dict] = []

    def dock(self, receptor_pdbqt: str, ligand_pdbqt: str, box, timeout_seconds: float, seed: int) -> float:
        self.calls.append(
            {
                "receptor": receptor_pdbqt,
                "ligand": ligand_pdbqt,
                "box": box,
                "timeout": timeout_seconds,
                "seed": seed,
            }
        )
        if self.error:
            raise self.error
        return self.scores.pop(0)


def test_parse_mutation_shorthand():
    assert parse_mutation_shorthand("D30N") == ("D", 30, "N")


def test_ligand_coords_ignore_water_and_use_hetatm():
    coords = ligand_hetatm_coords(LIGAND_PDB)
    assert coords == ((0.0, 0.0, 0.0), (10.0, 0.0, 0.0))


def test_grid_box_is_ligand_bounds_plus_five_angstrom_padding():
    box = grid_box_from_pdb(LIGAND_PDB, padding=5.0)
    assert box.center == pytest.approx((5.0, 0.0, 0.0))
    assert box.size == pytest.approx((20.0, 10.0, 10.0))


def test_grid_box_without_ligand_raises_no_binding_site():
    with pytest.raises(NoBindingSiteError):
        grid_box_from_pdb(APO_PDB)


def test_protein_atoms_only_strips_heteroatoms():
    protein = protein_atoms_only(LIGAND_PDB)
    assert "HETATM" not in protein
    assert "CA" in protein


def test_score_reports_mutant_minus_wildtype_delta():
    dock = FakeDock(scores=[-8.0, -6.5])
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=dock,
    )
    result = client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)
    assert result.status == "scored"
    assert result.method == "docking"
    assert result.delta_score == pytest.approx(1.5)
    assert len(dock.calls) == 2
    assert dock.calls[0]["receptor"] != dock.calls[1]["receptor"]
    assert dock.calls[0]["ligand"] == "LIGAND_PDBQT"


def test_apo_structure_raises_no_binding_site_not_a_fake_box():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(scores=[-1.0, -1.0]),
    )
    with pytest.raises(NoBindingSiteError):
        client.score(pdb_text=APO_PDB, smiles="CCO", mutation="D30N", position=30)


def test_receptor_prep_failure_raises_receptor_prep_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(error=ReceptorPrepError("bad protonation")),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(scores=[-1.0, -1.0]),
    )
    with pytest.raises(ReceptorPrepError):
        client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)


def test_mutant_placement_failure_raises_mutant_placement_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(error=MutantPlacementError("no rotamer")),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(scores=[-1.0, -1.0]),
    )
    with pytest.raises(MutantPlacementError):
        client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)


def test_unembeddable_smiles_raises_ligand_prep_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(error=LigandPrepError("embed failed")),
        dock_engine=FakeDock(scores=[-1.0, -1.0]),
    )
    with pytest.raises(LigandPrepError):
        client.score(pdb_text=LIGAND_PDB, smiles="not-a-smiles", mutation="D30N", position=30)


def test_vina_timeout_raises_timeout_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(error=VinaTimeoutError("dock hung")),
    )
    with pytest.raises(VinaTimeoutError):
        client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)


def test_vina_run_failure_raises_vina_run_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(error=VinaRunError("maps failed")),
    )
    with pytest.raises(VinaRunError):
        client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)


def test_specific_errors_are_vina_errors_for_orchestrator_catch():
    for exc in (
        ReceptorPrepError("x"),
        MutantPlacementError("x"),
        LigandPrepError("x"),
        NoBindingSiteError("x"),
        VinaTimeoutError("x"),
        VinaRunError("x"),
    ):
        assert isinstance(exc, VinaError)


def test_pdbfixer_places_mutant_sidechain_without_moving_backbone():
    from secondlook.vina_dock import PdbFixerMutantBuilder

    builder = PdbFixerMutantBuilder()
    mutant = builder.place_sidechain(PEPTIDE_PDB, "D2N", 2)
    ca_wt = [line[30:54] for line in PEPTIDE_PDB.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CA" and int(line[22:26]) == 2][0]
    ca_mut = [line[30:54] for line in mutant.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CA" and int(line[22:26]) == 2][0]
    assert ca_mut == ca_wt
    res2 = [line for line in mutant.splitlines() if line.startswith("ATOM") and int(line[22:26]) == 2]
    assert any("ASN" in line[17:20] for line in res2)
    assert any(line[12:16].strip() == "ND2" for line in res2)


def test_meeko_ligand_prep_rejects_unembeddable_smiles():
    from secondlook.vina_dock import MeekoLigandPreparer

    with pytest.raises(LigandPrepError):
        MeekoLigandPreparer().to_pdbqt("not-a-smiles")


def test_score_binding_maps_vina_errors_to_binding_unavailable_message():
    from secondlook.binding import BINDING_UNAVAILABLE_MESSAGE, score_binding
    from secondlook.candidates import DrugCandidate
    from secondlook.mutation_validation import MutationValidationResult
    from secondlook.structure import StructureResult

    validation = MutationValidationResult(
        status="valid",
        gene="TEST",
        hgvs_normalized="p.Asp30Asn",
        uniprot_accession="P00000",
        isoform_note="canonical",
        position=30,
        reference_residue_expected="D",
        reference_residue_actual="D",
        mutation_type="missense",
        error_message=None,
        wildtype_sequence="D" + "A" * 40,
        mutant_sequence="N" + "A" * 40,
    )
    structure = StructureResult(
        status="found",
        source="PDB",
        id="2Z4O",
        plddt_at_residue=90.0,
        plddt_global=90.0,
        reliability_flag="high",
        ligand_bound=True,
        annotated_position=30,
        pdb_text=LIGAND_PDB,
    )
    candidate = DrugCandidate(
        name="EXAMPLE",
        source="DGIdb",
        target_tier="exact_protein",
        approved=True,
        smiles="CCO",
        smiles_source="PubChem",
    )

    class RaisingVina:
        def score(self, **kwargs) -> BindingScore:
            raise ReceptorPrepError("bad atoms")

    result = score_binding(
        validation,
        structure,
        candidate,
        mcsm_client=type("M", (), {"submit": staticmethod(lambda **k: (_ for _ in ()).throw(Exception("unused")))})(),
        vina_client=RaisingVina(),
        het_resolver=type("H", (), {"resolve": staticmethod(lambda *a: None)})(),
        sleeper=lambda _: None,
    )
    assert result.status == "unavailable"
    # A receptor-prep failure must say so. It used to collapse into the generic
    # "mCSM-lig and docking were both unavailable or inapplicable", which reads
    # as a statement about the molecule rather than a structure-handling bug.
    assert result.reason_code == "receptor_prep_failed"
    assert "receptor structure could not be prepared" in result.error_message
    assert "bad atoms" in result.error_message, "the underlying cause must survive"
    assert result.error_message != BINDING_UNAVAILABLE_MESSAGE


def test_docking_timeout_is_reported_as_a_limit_not_a_refusal():
    """A timeout is a setting, not a property of the chemistry.

    This is the failure that made 7 of 9 gold-standard cases look like docking
    had refused them, when the calculation was simply cut off. See ISSUES.md #17.
    """
    from secondlook.binding import score_binding
    from secondlook.candidates import DrugCandidate
    from secondlook.mutation_validation import MutationValidationResult
    from secondlook.structure import StructureResult
    from secondlook.vina_dock import VinaTimeoutError

    validation = MutationValidationResult(
        status="valid", gene="TEST", hgvs_normalized="p.Asp30Asn",
        uniprot_accession="P00000", isoform_note="canonical", position=30,
        reference_residue_expected="D", reference_residue_actual="D",
        mutation_type="missense", error_message=None,
        wildtype_sequence="D" + "A" * 40, mutant_sequence="N" + "A" * 40,
    )
    structure = StructureResult(
        status="found", source="PDB", id="2Z4O", plddt_at_residue=90.0,
        plddt_global=90.0, reliability_flag="high", ligand_bound=True,
        annotated_position=30, pdb_text=LIGAND_PDB,
    )
    candidate = DrugCandidate(
        name="EXAMPLE", source="DGIdb", target_tier="exact_protein",
        approved=True, smiles="CCO", smiles_source="PubChem",
    )

    class TimingOutVina:
        timeout_seconds = 900.0

        def score(self, **kwargs) -> BindingScore:
            raise VinaTimeoutError("Vina docking exceeded 900.0s")

    result = score_binding(
        validation, structure, candidate,
        mcsm_client=type("M", (), {"submit": staticmethod(
            lambda **k: (_ for _ in ()).throw(Exception("unused")))})(),
        vina_client=TimingOutVina(),
        het_resolver=type("H", (), {"resolve": staticmethod(lambda *a: None)})(),
        sleeper=lambda _: None,
    )
    assert result.status == "unavailable"
    assert result.reason_code == "docking_timeout"
    assert "900s limit" in result.error_message
    assert "NOT a statement about the molecule" in result.error_message


@pytest.mark.integration
def test_9c5s_tp53_r175h_cannot_place_uniprot_residue_on_fragment():
    import httpx

    pdb_text = httpx.get("https://files.rcsb.org/download/9C5S.pdb", timeout=60.0).text
    client = VinaDockClient(timeout_seconds=30.0, seed=1, exhaustiveness=1)
    with pytest.raises((MutantPlacementError, NoBindingSiteError)):
        client.score(pdb_text=pdb_text, smiles="Fc1c[nH]c(=O)[nH]c1=O", mutation="R175H", position=175)


@pytest.mark.integration
def test_live_vina_delta_on_2z4o_with_small_ligand():
    import httpx

    # 9C5S (Step 3's TP53 hit) is a 14-residue fragment without UniProt 175 and
    # without a drug-like HET ligand, so it cannot host a defensible R175H dock.
    # 2Z4O is the verified ligand-bound complex (Prompt E example) and is the
    # smallest real structure where WT-vs-mutant Vina can actually run.
    #
    # Residue 130, NOT 30. 2Z4O is a two-chain structure: chain A is residues
    # 1-99, chain B is 101-199, and ligand 065 sits in chain B. This test
    # originally used D30N, which was fine until select_docking_chain() began
    # requiring the mutated residue and the ligand to share a chain — after
    # which D30N is correctly *refused* (see the companion test below), because
    # a grid box around a chain-B ligand says nothing about a chain-A mutation.
    # D130N is 2.76 A from ligand 065 in the same chain, so it is a real dock.
    pdb_text = httpx.get("https://files.rcsb.org/download/2Z4O.pdb", timeout=60.0).text
    client = VinaDockClient(timeout_seconds=300.0, seed=1, exhaustiveness=1)
    result = client.score(
        pdb_text=pdb_text, smiles="Fc1c[nH]c(=O)[nH]c1=O", mutation="D130N", position=130
    )
    assert result.status == "scored"
    assert result.method == "docking"
    assert result.delta_score is not None
    assert result.delta_score == result.delta_score  # not NaN

    # Deliberately no assertion on magnitude or sign. The measured value here is
    # ~0.05 kcal/mol, far below docking's own 1.0-2.19 error bar; asserting a
    # direction would be encoding noise as an expectation. What this test
    # proves is that the WT-vs-mutant path runs end to end on a real structure
    # — which matters, because the gold-standard run currently scores nothing,
    # and this distinguishes "the docking integration is broken" (it is not)
    # from "those particular structures fail meeko receptor prep" (they do).


@pytest.mark.integration
def test_live_vina_refuses_a_mutation_in_a_different_chain_from_the_ligand():
    """The chain rule, exercised against a real two-chain structure.

    2Z4O residue 30 is in chain A; ligand 065 is in chain B. Docking there
    would build a grid box around a pocket the mutation cannot reach and return
    a confident-looking delta for it. Refusing is the correct answer, and it
    must stay refused — this is the exact case the rule was written for.
    """
    import httpx

    pdb_text = httpx.get("https://files.rcsb.org/download/2Z4O.pdb", timeout=60.0).text
    client = VinaDockClient(timeout_seconds=180.0, seed=1, exhaustiveness=1)
    with pytest.raises(NoBindingSiteError) as exc:
        client.score(
            pdb_text=pdb_text, smiles="Fc1c[nH]c(=O)[nH]c1=O", mutation="D30N", position=30
        )
    # The message must name both sides, or a reader cannot tell why it refused.
    assert "residue in ['A']" in str(exc.value)
    assert "ligand in ['B']" in str(exc.value)


# --- Regressions from the first successful real docking run -------------------
# Four independent defects, all found by running BRAF V600E end to end. Each one
# alone made Vina scoring impossible on real structures, and none was reachable
# by the existing unit tests because they used small hand-written PDB fixtures.

from secondlook.vina_dock import (  # noqa: E402
    chains_with_ligand,
    chains_with_residue,
    extract_chain,
    is_primary_altloc,
    select_docking_chain,
    strip_terminal_oxt,
)

TWO_CHAIN_SPLIT_PDB = """\
ATOM      1  CA  VAL A 600       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CA  ALA B 600       5.000   0.000   0.000  1.00 20.00           C
HETATM    3  C1  LIG B 900       6.000   0.000   0.000  1.00 20.00           C
HETATM    4  C2  LIG B 900       7.000   0.000   0.000  1.00 20.00           C
END
"""

ALTLOC_PDB = """\
ATOM      1  CA AVAL A 600       0.000   0.000   0.000  0.60 20.00           C
ATOM      2  CA BVAL A 600       0.200   0.000   0.000  0.40 20.00           C
ATOM      3  CB  VAL A 600       1.500   0.000   0.000  1.00 20.00           C
HETATM    4  C1  LIG A 900       6.000   0.000   0.000  1.00 20.00           C
END
"""

OXT_PDB = """\
ATOM      1  CA  SER A 722       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  C   SER A 722       1.520   0.000   0.000  1.00 20.00           C
ATOM      3  OXT SER A 722       1.830   1.000   0.000  1.00 20.00           O
END
"""


def test_alternate_conformations_are_dropped():
    """Duplicate atoms make RDKit see an over-bonded carbon and break meeko prep."""
    assert is_primary_altloc("ATOM      1  CA AVAL A 600") is True
    assert is_primary_altloc("ATOM      2  CA BVAL A 600") is False
    out = protein_atoms_only(ALTLOC_PDB)
    ca_lines = [l for l in out.splitlines() if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    assert len(ca_lines) == 1


def test_altloc_column_is_blanked_not_just_filtered():
    """Downstream parsers must not see a stale altLoc code on a kept atom."""
    out = protein_atoms_only(ALTLOC_PDB)
    for line in out.splitlines():
        if line.startswith("ATOM"):
            assert line[16] == " "


def test_altloc_duplicates_do_not_skew_the_grid_box():
    coords = ligand_hetatm_coords(ALTLOC_PDB)
    assert len(coords) == 1


def test_terminal_oxt_is_stripped_before_receptor_prep():
    """PDBFixer's OXT lands ~1.83 A from CA, inside RDKit's proximity-bond cutoff."""
    out = strip_terminal_oxt(OXT_PDB)
    assert "OXT" not in out
    # Everything else survives.
    assert "CA  SER A 722" in out
    assert "C   SER A 722" in out


def test_oxt_stripping_leaves_non_terminal_atoms_alone():
    assert strip_terminal_oxt(PEPTIDE_PDB).count("ATOM") == PEPTIDE_PDB.count("ATOM")


# --- Chain selection ---------------------------------------------------------


def test_chains_with_residue_finds_every_copy():
    assert chains_with_residue(TWO_CHAIN_SPLIT_PDB, 600) == ["A", "B"]
    assert chains_with_residue(TWO_CHAIN_SPLIT_PDB, 999) == []


def test_chains_with_ligand_ignores_waters_and_ions():
    assert chains_with_ligand(LIGAND_PDB) == {"A": 2}  # HOH excluded
    assert chains_with_ligand(APO_PDB) == {}


def test_selects_the_chain_holding_both_residue_and_ligand():
    assert select_docking_chain(TWO_CHAIN_SPLIT_PDB, 600) == "B"


def test_refuses_when_residue_and_ligand_are_in_different_chains():
    """The 3OG7 case: BRAF 600 resolved only in chain B, ligand only in chain A."""
    split = """\
ATOM      1  CA  GLU B 600       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  C1  LIG A 900       6.000   0.000   0.000  1.00 20.00           C
END
"""
    with pytest.raises(NoBindingSiteError, match="both residue 600"):
        select_docking_chain(split, 600)


def test_refuses_when_no_chain_resolves_the_residue():
    with pytest.raises(NoBindingSiteError, match="resolves residue"):
        select_docking_chain(LIGAND_PDB, 600)


def test_extract_chain_keeps_only_that_chain():
    out = extract_chain(TWO_CHAIN_SPLIT_PDB, "B")
    assert "ALA B 600" in out
    assert "VAL A 600" not in out
    assert "LIG B 900" in out


def test_grid_box_uses_only_the_selected_chain_ligand():
    """Otherwise the box can sit on another chain's pocket entirely."""
    box = grid_box_from_pdb(extract_chain(TWO_CHAIN_SPLIT_PDB, "B"))
    assert box.center[0] == pytest.approx(6.5)


# --- Structures deposited as the mutant --------------------------------------
# Well-studied oncogenic mutations are frequently crystallized AS the mutant
# (3OG7 is BRAF with GLU already at 600). Assuming wild-type input would refuse
# exactly the best-characterised structures for the cases that matter most.

MUT_FORM_PDB = """\
ATOM      1  N   GLU A 600       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  GLU A 600       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   GLU A 600       2.009   1.420   0.000  1.00 20.00           C
HETATM    4  C1  LIG A 900       6.000   0.000   0.000  1.00 20.00           C
HETATM    5  C2  LIG A 900       7.000   0.000   0.000  1.00 20.00           C
END
"""

WT_FORM_PDB = MUT_FORM_PDB.replace("GLU A 600", "VAL A 600")


class _RecordingBuilder:
    def __init__(self):
        self.calls = []

    def place_sidechain(self, pdb_text, mutation, position):
        self.calls.append(mutation)
        return pdb_text


class _CountingDock:
    def __init__(self, scores):
        self.scores = list(scores)

    def dock(self, receptor, ligand, box, timeout, seed):
        return self.scores.pop(0)


def _client(builder, scores):
    return VinaDockClient(
        receptor_preparer=type("R", (), {"to_pdbqt": lambda self, t: "REC"})(),
        mutant_builder=builder,
        ligand_preparer=type("L", (), {"to_pdbqt": lambda self, s: "LIG"})(),
        dock_engine=_CountingDock(scores),
        protonator=None,
    )


def test_wildtype_structure_builds_the_mutant_forward():
    builder = _RecordingBuilder()
    result = _client(builder, [-8.0, -7.0]).score(
        pdb_text=WT_FORM_PDB, smiles="CCO", mutation="V600E", position=600
    )
    assert builder.calls == ["V600E"]
    assert result.delta_score == pytest.approx(1.0)  # mutant minus wild-type


def test_mutant_structure_builds_the_wildtype_in_reverse():
    """The structure is already V600E, so the wild-type is what must be built."""
    builder = _RecordingBuilder()
    result = _client(builder, [-8.0, -7.0]).score(
        pdb_text=MUT_FORM_PDB, smiles="CCO", mutation="V600E", position=600
    )
    assert builder.calls == ["E600V"]
    assert result.status == "scored"


def test_delta_stays_mutant_minus_wildtype_whichever_form_was_deposited():
    """Sign convention must not flip based on which form happened to be crystallized.

    Scores are keyed off the residue actually present in each prepared receptor
    (GLU = mutant, VAL = wild-type), so both paths describe the same physical
    pair and must produce the same delta.
    """

    class SubstitutingBuilder:
        """Rewrites residue 600 to the mutation's target, as PDBFixer would."""

        def place_sidechain(self, pdb_text, mutation, position):
            target = {"E": "GLU", "V": "VAL"}[mutation[-1]]
            out = []
            for line in pdb_text.splitlines():
                if line.startswith("ATOM") and line[22:26].strip() == str(position):
                    line = line[:17] + target + line[20:]
                out.append(line)
            return "\n".join(out) + "\n"

    class ResidueAwareDock:
        """GLU600 receptor docks at -7.0, VAL600 at -8.0."""

        def dock(self, receptor, ligand, box, timeout, seed):
            return -7.0 if "GLU A 600" in receptor else -8.0

    def client():
        return VinaDockClient(
            receptor_preparer=type("R", (), {"to_pdbqt": lambda self, t: t})(),
            mutant_builder=SubstitutingBuilder(),
            ligand_preparer=type("L", (), {"to_pdbqt": lambda self, s: "LIG"})(),
            dock_engine=ResidueAwareDock(),
            protonator=None,
        )

    from_wt = client().score(
        pdb_text=WT_FORM_PDB, smiles="CCO", mutation="V600E", position=600
    )
    from_mut = client().score(
        pdb_text=MUT_FORM_PDB, smiles="CCO", mutation="V600E", position=600
    )
    # mutant (-7.0) minus wild-type (-8.0) = +1.0, from either deposited form.
    assert from_wt.delta_score == pytest.approx(1.0)
    assert from_mut.delta_score == pytest.approx(1.0)


def test_residue_matching_neither_form_is_refused():
    with pytest.raises(MutantPlacementError, match="neither the"):
        _client(_RecordingBuilder(), [-8.0, -7.0]).score(
            pdb_text=MUT_FORM_PDB.replace("GLU A 600", "TRP A 600"),
            smiles="CCO",
            mutation="V600E",
            position=600,
        )


def test_selected_chain_is_reported_on_the_score():
    result = _client(_RecordingBuilder(), [-8.0, -7.0]).score(
        pdb_text=WT_FORM_PDB, smiles="CCO", mutation="V600E", position=600
    )
    assert result.chain == "A"


# --- Mutations too far from the pocket to score ------------------------------
# The gold-standard run showed docking deltas collapse to noise beyond direct
# contact range: EGFR T790M at 3.5 A gave the largest delta and the correct
# direction, while BRAF V600E (12.2 A) and KIT D816V (16.6 A) gave |delta| < 0.06
# regardless of their known clinical direction. Emitting a number there would be
# reporting noise as signal.

from secondlook.vina_dock import (  # noqa: E402
    MUTATION_CONTACT_MAX_ANGSTROM,
    MutationOutsidePocketError,
    residue_min_distance_to_ligand,
    resolve_overvalence,
)

FAR_MUTATION_PDB = """\
ATOM      1  N   VAL A 600      50.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  VAL A 600      51.000   0.000   0.000  1.00 20.00           C
HETATM    3  C1  LIG A 900       0.000   0.000   0.000  1.00 20.00           C
HETATM    4  C2  LIG A 900       1.000   0.000   0.000  1.00 20.00           C
END
"""


def test_distance_to_ligand_is_measured():
    assert residue_min_distance_to_ligand(FAR_MUTATION_PDB, 600) == pytest.approx(49.0)
    assert residue_min_distance_to_ligand(FAR_MUTATION_PDB, 999) is None
    assert residue_min_distance_to_ligand(APO_PDB, 30) is None


def test_mutation_beyond_contact_range_is_refused_before_docking():
    docked = []

    class CountingDock:
        def dock(self, *a, **k):
            docked.append(1)
            return -8.0

    client = VinaDockClient(
        receptor_preparer=type("R", (), {"to_pdbqt": lambda self, t: "REC"})(),
        mutant_builder=type("M", (), {"place_sidechain": lambda self, t, m, p: t})(),
        ligand_preparer=type("L", (), {"to_pdbqt": lambda self, s: "LIG"})(),
        dock_engine=CountingDock(),
        protonator=None,
    )
    with pytest.raises(MutationOutsidePocketError, match="beyond the"):
        client.score(pdb_text=FAR_MUTATION_PDB, smiles="CCO", mutation="V600E", position=600)
    assert docked == [], "must refuse before spending two docking runs"


def test_out_of_pocket_reports_its_own_message_not_a_failure_message():
    """An informative negative result must not be described as a scoring failure."""
    from secondlook.binding import BINDING_UNAVAILABLE_MESSAGE, MUTATION_OUTSIDE_POCKET_MESSAGE

    rendered = MUTATION_OUTSIDE_POCKET_MESSAGE.format(distance=16.6, threshold=8.0)
    assert "16.6" in rendered
    assert "allosteric" in rendered
    assert rendered != BINDING_UNAVAILABLE_MESSAGE
    assert "unavailable or inapplicable" not in rendered


def test_contact_threshold_admits_direct_contact_distances():
    """EGFR T790M sits at 3.5 A and must remain scoreable."""
    assert MUTATION_CONTACT_MAX_ANGSTROM > 3.5


def test_resolve_overvalence_is_a_noop_on_clean_structures():
    assert resolve_overvalence(PEPTIDE_PDB).count("ATOM") == PEPTIDE_PDB.count("ATOM")


# --- Receptor preparer fallback ----------------------------------------------
# meeko's residue templates reject PDBFixer-built mutant side chains
# ("Updated 1 H positions but deleted 4"), which blocked ABL1 T315I (5HU9) and
# both ALK cases (4Z55) — three of the nine gold-standard cases. Open Babel does
# not use residue templates and prepares them.

from secondlook.vina_dock import (  # noqa: E402
    FallbackReceptorPreparer,
    OpenBabelReceptorPreparer,
    prepare_receptor_pair,
)


class _Preparer:
    """Records what it was asked to prepare; optionally refuses some inputs."""

    def __init__(self, name, refuse_containing=None):
        self.name = name
        self.calls = []
        self._refuse = refuse_containing

    def to_pdbqt(self, pdb_text):
        self.calls.append(pdb_text)
        if self._refuse and self._refuse in pdb_text:
            raise ReceptorPrepError(f"{self.name} refuses this structure")
        return f"{self.name}:{pdb_text}"


def test_fallback_is_not_used_when_the_primary_succeeds():
    """meeko's atom typing is the reference Vina was parameterised against, so
    it must win wherever it works."""
    primary, fallback = _Preparer("meeko"), _Preparer("obabel")
    prep = FallbackReceptorPreparer(primary, fallback)
    assert prep.to_pdbqt("good") == "meeko:good"
    assert prep.used_fallback is False
    assert fallback.calls == []


def test_fallback_is_used_when_the_primary_refuses():
    primary, fallback = _Preparer("meeko", refuse_containing="MUT"), _Preparer("obabel")
    prep = FallbackReceptorPreparer(primary, fallback)
    assert prep.to_pdbqt("MUT") == "obabel:MUT"
    assert prep.used_fallback is True


def test_both_failing_reports_both_reasons():
    """Hiding the primary's message behind the fallback's would lose the actual
    diagnosis — the primary is the one whose failure mode is understood."""
    primary = _Preparer("meeko", refuse_containing="BAD")
    fallback = _Preparer("obabel", refuse_containing="BAD")
    with pytest.raises(ReceptorPrepError) as exc:
        FallbackReceptorPreparer(primary, fallback).to_pdbqt("BAD")
    assert "meeko refuses" in str(exc.value)
    assert "obabel refuses" in str(exc.value)


# --- The pair must travel the same path --------------------------------------


def test_a_mutant_needing_fallback_drags_the_wildtype_with_it():
    """This is the whole point. meeko and Open Babel assign different AutoDock
    types and different partial charges. A wild-type typed by one and a mutant
    typed by the other yields a delta contaminated by the preparer difference —
    and that contamination is the same order as the 0.02-0.56 kcal/mol signal
    being measured."""
    primary = _Preparer("meeko", refuse_containing="MUT")
    fallback = _Preparer("obabel")
    wt, mt = prepare_receptor_pair(FallbackReceptorPreparer(primary, fallback), "WT", "MUT")
    assert wt == "obabel:WT"
    assert mt == "obabel:MUT"
    assert primary.calls == ["MUT"], "wild-type must not be prepared by the primary"


def test_both_use_the_primary_when_neither_needs_the_fallback():
    primary, fallback = _Preparer("meeko"), _Preparer("obabel")
    wt, mt = prepare_receptor_pair(FallbackReceptorPreparer(primary, fallback), "WT", "MUT")
    assert (wt, mt) == ("meeko:WT", "meeko:MUT")
    assert fallback.calls == []


def test_the_mutant_is_attempted_first():
    """The mutant is the one that fails, so trying it first decides the path
    before any work is spent on a wild-type that may have to be redone."""
    primary, fallback = _Preparer("meeko"), _Preparer("obabel")
    prepare_receptor_pair(FallbackReceptorPreparer(primary, fallback), "WT", "MUT")
    assert primary.calls[0] == "MUT"


def test_a_wildtype_only_failure_also_forces_the_pair_onto_the_fallback():
    primary = _Preparer("meeko", refuse_containing="WT")
    fallback = _Preparer("obabel")
    wt, mt = prepare_receptor_pair(FallbackReceptorPreparer(primary, fallback), "WT", "MUT")
    assert wt == "obabel:WT" and mt == "obabel:MUT"


def test_an_explicitly_injected_preparer_is_used_as_given():
    """A caller who supplies a preparer means it; this is not the place to
    second-guess an explicit injection."""
    p = _Preparer("custom")
    wt, mt = prepare_receptor_pair(p, "WT", "MUT")
    assert (wt, mt) == ("custom:WT", "custom:MUT")


def test_open_babel_receptor_strips_heteroatoms():
    """The ligand is docked separately. Leaving the crystallised copy in the
    receptor docks into an already-occupied pocket — wrong answer, no error."""
    pdb = (
        "ATOM      1  N   THR A 315      10.000  10.000  10.000  1.00 20.00           N\n"
        "HETATM 5000  C1  66K A 601      12.000  12.000  12.000  1.00 20.00           C\n"
        "HETATM 5001  O   HOH A 700      14.000  14.000  14.000  1.00 20.00           O\n"
    )
    out = OpenBabelReceptorPreparer._protein_only(pdb, OpenBabelReceptorPreparer._DISCARD_HET)
    assert "THR" in out
    assert "66K" not in out
    assert "HOH" not in out


def test_missing_open_babel_names_the_install_command():
    """A fallback that fails with a bare FileNotFoundError is worse than no
    fallback — the reader cannot tell what to do about it."""
    import shutil

    real = shutil.which
    shutil.which = lambda name: None if name == "obabel" else real(name)
    try:
        with pytest.raises(ReceptorPrepError) as exc:
            OpenBabelReceptorPreparer().to_pdbqt("ATOM      1  N   THR A 315\n")
        assert "brew install open-babel" in str(exc.value)
    finally:
        shutil.which = real


# --- Reproducibility of structure preparation ---------------------------------
#
# Vina is invoked with cpu=1 and a fixed seed, so it is deterministic given
# identical input files. Docking deltas nonetheless failed to reproduce at a
# fixed seed (~0.4 kcal/mol spread, the same order as the effects being
# measured). The variation came from structure preparation, upstream of Vina:
# PDBFixer's addMissingAtoms minimization, OpenMM Modeller's random hydrogen
# placement, and multi-threaded force summation on the CPU platform.

from secondlook.vina_dock import deterministic_cpu_platform  # noqa: E402


def test_deterministic_cpu_platform_pins_single_thread():
    platform = deterministic_cpu_platform()
    if platform is None:
        pytest.skip("openmm is not installed")
    assert platform.getPropertyDefaultValue("Threads") == "1"


def test_pdbfixer_complete_restores_the_global_rng_state():
    """Seeding is required for reproducible hydrogen placement, but the process-wide
    RNG must be left as it was found -- callers do not expect a structure-prep call
    to reset their random stream."""
    import random

    pytest.importorskip("pdbfixer")
    from secondlook.vina_dock import _pdbfixer_complete

    random.seed(12345)
    before = random.random()
    random.seed(12345)
    try:
        _pdbfixer_complete(SINGLE_RESIDUE_PDB, None, None)
    except ReceptorPrepError:
        pytest.skip("pdbfixer could not complete the minimal test structure")
    assert random.random() == before


SINGLE_RESIDUE_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       1.988  -0.773  -1.199  1.00 20.00           C
END
"""


@pytest.mark.integration
def test_structure_preparation_is_byte_identical_across_runs():
    """The regression test for the fixed-seed reproducibility bug.

    Marked integration because it runs PDBFixer and OpenMM on a real structure,
    which takes minutes. Requires 5HU9 in the local PDB cache.
    """
    from pathlib import Path

    pytest.importorskip("pdbfixer")
    pdb_path = Path("/tmp/athena-pdb-cache/5HU9.pdb")
    if not pdb_path.exists():
        pytest.skip("5HU9 not cached; run validation/verify_proximity.py first")

    from secondlook.vina_dock import (
        PdbFixerMutantBuilder,
        PdbFixerProtonator,
        extract_chain,
        select_docking_chain,
    )

    text = pdb_path.read_text()
    chain = select_docking_chain(text, 315)
    protein = protein_atoms_only(extract_chain(text, chain))

    protonator = PdbFixerProtonator()
    assert protonator.protonate(protein) == protonator.protonate(protein)

    builder = PdbFixerMutantBuilder()
    first = builder.place_sidechain(protein, "T315I", 315)
    second = builder.place_sidechain(protein, "T315I", 315)
    assert first == second


# --- protomer selection for docking ------------------------------------------
#
# ISSUES.md #18. `measure_proximity` was fixed to report the closest approach
# across protomers; the docking path was not, and still picked a chain by ligand
# size. On 8C7X that chose BRAF V600's far copy (12.2 A) over its near one
# (7.2 A), so the case was refused as "outside the 8 A contact range" while the
# proximity signal in the SAME result object reported 7.2 A.

TWO_PROTOMER_PDB = """\
ATOM      1  CA  VAL A 600      20.000   0.000   0.000  1.00 20.00           C
HETATM    2  C1  LIG A 700       0.000   0.000   0.000  1.00 20.00           C
HETATM    3  C2  LIG A 700       1.000   0.000   0.000  1.00 20.00           C
HETATM    4  C3  LIG A 700       2.000   0.000   0.000  1.00 20.00           C
ATOM      5  CA  VAL B 600       3.000   0.000   0.000  1.00 20.00           C
HETATM    6  C1  LIG B 700       0.000   0.000   0.000  1.00 20.00           C
HETATM    7  C2  LIG B 700       1.000   0.000   0.000  1.00 20.00           C
END
"""


def test_docking_chain_is_the_protomer_where_the_mutation_is_closest():
    """Chain A holds the bigger ligand but the residue is 18 A away; chain B's
    residue is 1 A from its ligand. The near copy must win, or the contact gate
    refuses a case it should have scored."""
    from secondlook.vina_dock import select_docking_chain

    assert select_docking_chain(TWO_PROTOMER_PDB, 600) == "B"


def test_min_residue_ligand_distance_takes_the_closest_protomer():
    from secondlook.vina_dock import min_residue_ligand_distance

    # chain B: residue at x=3.0, nearest chain-B ligand atom at x=1.0
    assert min_residue_ligand_distance(TWO_PROTOMER_PDB, 600) == pytest.approx(2.0)


def test_ligand_size_still_breaks_ties():
    """The original rule — prefer the biggest HET group, so a buffer component
    does not win — must survive when the distances are equal."""
    from secondlook.vina_dock import select_docking_chain

    tied = """\
ATOM      1  CA  VAL A 600       3.000   0.000   0.000  1.00 20.00           C
HETATM    2  C1  LIG A 700       0.000   0.000   0.000  1.00 20.00           C
HETATM    3  C2  LIG A 700      50.000   0.000   0.000  1.00 20.00           C
HETATM    4  C3  LIG A 700      51.000   0.000   0.000  1.00 20.00           C
ATOM      5  CA  VAL B 600       3.000   0.000   0.000  1.00 20.00           C
HETATM    6  C1  LIG B 700       0.000   0.000   0.000  1.00 20.00           C
END
"""
    assert select_docking_chain(tied, 600) == "A"


def test_docking_gate_and_proximity_cannot_disagree():
    """The invariant the bug violated: the distance the contact gate uses must be
    the distance the proximity signal reports."""
    from secondlook.proximity import build_proximity_signal
    from secondlook.vina_dock import (
        extract_chain,
        min_residue_ligand_distance,
        residue_min_distance_to_ligand,
        select_docking_chain,
    )

    chain = select_docking_chain(TWO_PROTOMER_PDB, 600)
    gate_distance = residue_min_distance_to_ligand(extract_chain(TWO_PROTOMER_PDB, chain), 600)
    across_protomers = min_residue_ligand_distance(TWO_PROTOMER_PDB, 600)
    assert gate_distance == pytest.approx(across_protomers)

    signal = build_proximity_signal(
        across_protomers,
        structure_id="TEST",
        structure_source="PDB",
        protomer_distances={"A": 18.0, "B": 2.0},
    )
    assert signal.distance_angstrom == pytest.approx(gate_distance)
