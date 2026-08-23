"""Tests for the ask.py CLI's Tier 1 key routing.

The bug these exist to prevent (ISSUES.md §16): Tier 1 matches `hgvs_p` by exact
equality, and CIViC stores it RefSeq-prefixed and three-letter
(`NP_005148.2:p.Thr315Ile`). Tier 2's mutation format is shorthand (`T315I`).
Sending shorthand as `hgvs_p` can never match, and the failure is silent -- Tier 1
returns `no_hit`, which reads as "nothing documented" rather than as an error.
"""

import pytest

from ask import HGVS_P, SHORTHAND, tier1_lookups, to_three_letter, variant_keys
from secondlook.tier1_adapter import Tier1RetrievalPolicy


class FakeResult:
    def __init__(self, rows):
        self.result_set = rows


class FakeGraph:
    def __init__(self, responses):
        self._responses = list(responses)
        self.calls: list[tuple[str, dict]] = []

    def query(self, cypher, params=None):
        self.calls.append((cypher, params or {}))
        return FakeResult(self._responses.pop(0) if self._responses else [])


EVIDENCE_ROW = (
    "Imatinib",
    "resistant",
    "B",
    "https://civicdb.org/evidence/12",
    "T315I confers resistance to imatinib.",
    12,
)


@pytest.mark.parametrize("mutation", ["T315I", "V600E", "G1202R", "R273*"])
def test_shorthand_is_routed_to_variant_name(mutation):
    assert SHORTHAND.match(mutation)
    assert variant_keys(mutation) == {"variant_name": mutation}


@pytest.mark.parametrize(
    "mutation",
    ["NP_005148.2:p.Thr315Ile", "p.Arg273His", "NP_000537.3:p.Arg273His"],
)
def test_hgvs_is_routed_to_the_mutation_keyword(mutation):
    assert not SHORTHAND.match(mutation)
    assert variant_keys(mutation) == {"mutation": mutation}


def test_shorthand_reaches_tier1_as_a_variant_name_query():
    """The regression assertion: the Cypher that runs must bind variant_name."""
    graph = FakeGraph(responses=[[EVIDENCE_ROW]])
    decision = Tier1RetrievalPolicy(graph=graph).decide(
        gene="ABL1", **variant_keys("T315I")
    )
    cypher, params = graph.calls[0]
    assert "$variant_name" in cypher
    assert params["variant_name"] == "T315I"
    assert decision.tier1_results, "evidence should have been found"


def test_shorthand_sent_as_hgvs_p_finds_nothing():
    """Documents the broken path, so the fix in run_tier2 has a failing case to
    turn green. This is what run_tier2 does today."""
    graph = FakeGraph(responses=[[], []])
    decision = Tier1RetrievalPolicy(graph=graph).decide(gene="ABL1", mutation="T315I")
    cypher, params = graph.calls[0]
    assert "$hgvs_p" in cypher
    assert params["hgvs_p"] == "T315I"
    assert decision.state == "no_hit"
    assert not decision.tier1_results


# --- routing for the variant types that are NOT point substitutions ----------
#
# CIViC keeps "Fusion", "Amplification", "Loss", "Overexpression" and friends in
# `name`, and leaves `hgvs_p` empty for all of them. Sending one as `hgvs_p`
# matches nothing and silently relaxes to Mode 2's gene-level query, which
# downgrades a `strong_hit` to a `weak_hit` -- verified live against the loaded
# graph: ETV6::NTRK3 "Fusion" returns strong_hit(4) by name and weak_hit(4) by
# hgvs_p, from the same four evidence items.


@pytest.mark.parametrize(
    "descriptor",
    ["Fusion", "Amplification", "Loss", "Overexpression", "EXPRESSION", "UNDEREXPRESSION"],
)
def test_non_missense_descriptors_route_to_variant_name(descriptor):
    assert not SHORTHAND.match(descriptor)
    assert not HGVS_P.search(descriptor)
    assert tier1_lookups(descriptor) == [{"variant_name": descriptor}]


@pytest.mark.parametrize(
    "mutation",
    ["p.Thr315Ile", "NP_005148.2:p.Thr315Ile", "NP_000537.3:p.Arg273His"],
)
def test_hgvs_is_recognised_and_used_verbatim(mutation):
    assert HGVS_P.search(mutation)
    assert tier1_lookups(mutation) == [{"mutation": mutation}]


def test_shorthand_tries_variant_name_then_three_letter_hgvs():
    assert tier1_lookups("T315I") == [
        {"variant_name": "T315I"},
        {"mutation": "p.Thr315Ile"},
    ]


@pytest.mark.parametrize(
    ("shorthand", "expected"),
    [("T315I", "p.Thr315Ile"), ("V600E", "p.Val600Glu"), ("G1202R", "p.Gly1202Arg")],
)
def test_three_letter_conversion(shorthand, expected):
    assert to_three_letter(shorthand) == expected


def test_nonsense_shorthand_has_no_three_letter_form():
    """`R273*` is a stop gain -- `*` is not a standard residue, so there is no
    three-letter HGVS to fall back to. It must still route by name."""
    assert to_three_letter("R273*") is None
    assert tier1_lookups("R273*") == [{"variant_name": "R273*"}]


def test_fusion_reaches_tier1_as_an_exact_name_match():
    """The regression case for the strong_hit downgrade."""
    graph = FakeGraph(responses=[[EVIDENCE_ROW]])
    decision = Tier1RetrievalPolicy(graph=graph).decide(
        gene="ETV6::NTRK3", **tier1_lookups("Fusion")[0]
    )
    cypher, params = graph.calls[0]
    assert "$variant_name" in cypher
    assert params["variant_name"] == "Fusion"
    assert decision.tier1_results
