"""Own tests for the variant triage service (SECB3203 Assignment 1)."""
import copy
import json
from pathlib import Path

import pytest
from src.a1_triage import (
    classify_variant,
    triage,
    decision_counts,
    impact_counts,
    gene_summary,
    rank_variants,
    results_to_csv,
    build_report,
)
from copy import deepcopy

#Parameters
PARAMS_PATH = Path(__file__).resolve().parents[1] / "data" / "params.json"
TEST_PARAMS = json.loads(PARAMS_PATH.read_text())
MIN_DEPTH = int(TEST_PARAMS["min_depth"])
MIN_QUALITY = int(TEST_PARAMS["min_quality"])
PRIORITY_GENE = TEST_PARAMS["priority_genes"][0]
OTHER_GENE = "NOTAGENE"            
GOOD_DEPTH = MIN_DEPTH + 50
GOOD_QUALITY = MIN_QUALITY + 50


def make_row(variant="v1", gene=PRIORITY_GENE, impact="high",
             depth=GOOD_DEPTH, quality=GOOD_QUALITY, pos="100"):
    """Build one variant row of text values, like the CSV gives."""
    return {"variant": variant, "gene": gene, "chrom": "1", "pos": str(pos),
            "ref": "A", "alt": "T", "depth": str(depth),
            "quality": str(quality), "impact": impact}


# Test 1 and 2: classify_variant() and triage()   (Fatini)

@pytest.mark.parametrize("gene, impact, expected", [
    (PRIORITY_GENE, "high", "prioritize"),    # rule 3
    (PRIORITY_GENE, "low", "review"),         # rule 4
    (OTHER_GENE, "moderate", "review"),       # rule 5
    (OTHER_GENE, "modifier", "exclude"),      # rule 6
])
def test_classify_variant_rule_matrix(gene, impact, expected):
    """Normal: rules 3 to 6."""
    decision, reason = classify_variant(make_row(gene=gene, impact=impact), TEST_PARAMS)
    assert decision == expected
    assert reason and "," not in reason


@pytest.mark.parametrize("depth, expected", [
    (MIN_DEPTH - 1, "exclude"), (MIN_DEPTH, "prioritize"), (MIN_DEPTH + 1, "prioritize"),
])
def test_classify_variant_depth_boundary(depth, expected):
    """Boundary: depth one below, at, one above."""
    decision, reason = classify_variant(make_row(depth=depth), TEST_PARAMS)
    assert decision == expected
    if expected == "exclude":
        assert "depth" in reason.lower()
        assert str(depth) in reason and str(MIN_DEPTH) in reason


@pytest.mark.parametrize("quality, expected", [
    (MIN_QUALITY - 1, "exclude"), (MIN_QUALITY, "prioritize"), (MIN_QUALITY + 1, "prioritize"),
])
def test_classify_variant_quality_boundary(quality, expected):
    """Boundary: quality one below, at, one above."""
    decision, reason = classify_variant(make_row(quality=quality), TEST_PARAMS)
    assert decision == expected
    if expected == "exclude":
        assert "quality" in reason.lower()
        assert str(quality) in reason and str(MIN_QUALITY) in reason


def test_classify_variant_rule_order_depth_first():
    """Awkward: depth and quality both fail, depth reports."""
    decision, reason = classify_variant(make_row(depth=1, quality=1), TEST_PARAMS)
    assert decision == "exclude"
    assert "depth" in reason.lower() and "quality" not in reason.lower()


def test_classify_variant_messy_text():
    """Awkward: spaces and case in gene and impact."""
    row = make_row(gene=f"  {PRIORITY_GENE} ",impact=" HIGH ")
    assert classify_variant(row, TEST_PARAMS)[0] == "prioritize"


def test_classify_variant_respects_params():
    """Other thresholds and genes give other answers."""
    stricter = dict(TEST_PARAMS, min_depth=GOOD_DEPTH + 1)
    assert classify_variant(make_row(), stricter)[0] == "exclude"
    other_genes = dict(TEST_PARAMS, priority_genes=[OTHER_GENE])
    assert classify_variant(make_row(), other_genes)[0] == "review"


def test_triage_empty_list():
    assert triage([], TEST_PARAMS) == []


def test_triage_order_fields_and_no_side_effects():
    rows = [make_row(variant="a"), make_row(variant="b", depth=1)]
    rows_before, params_before = copy.deepcopy(rows), copy.deepcopy(TEST_PARAMS)
    out = triage(rows, TEST_PARAMS)
    assert [r["variant"] for r in out] == ["a", "b"]
    assert out[0]["decision"] == "prioritize" and out[1]["decision"] == "exclude"
    assert all("reason" in r for r in out)
    assert out[0] is not rows[0]
    assert rows == rows_before and TEST_PARAMS == params_before


def test_classify_variant_invalid_depth_raises():
    """Invalid input."""
    with pytest.raises(ValueError):
        classify_variant(make_row(depth="abc"), TEST_PARAMS)



# Test 3 and 4: decision_counts(), impact_counts(), gene_summary()

DEPTH_LOW = 10
DEPTH_MID = 20
DEPTH_HIGH = 30
EXPECTED_DECISION_KEYS = ["prioritize", "review", "exclude"]
EXPECTED_IMPACT_KEYS = ["high", "moderate", "low", "modifier"]
EXPECTED_MEAN_DEPTH_BRCA2 = 25.0
EXPECTED_MEAN_DEPTH_TP53 = 10.0


def make_results():
    """A small hand-made results list (rows plus decision and reason)."""
    return [
        {
            **make_row(variant="r1", gene="BRCA2", impact="high", depth=DEPTH_HIGH),
            "decision": "prioritize",
            "reason": "priority gene brca2 with high impact",
        },
        {
            **make_row(variant="r2", gene="BRCA2", impact="low", depth=DEPTH_MID),
            "decision": "review",
            "reason": "priority gene brca2 with low impact",
        },
        {
            **make_row(variant="r3", gene="TP53", impact="modifier", depth=DEPTH_LOW),
            "decision": "exclude",
            "reason": "depth 10 below 12",
        },
        {
            **make_row(variant="r4", gene="EGFR", impact="moderate", depth=DEPTH_MID),
            "decision": "review",
            "reason": "non-priority gene egfr with moderate impact",
        },
    ]


def test_decision_counts_normal():
    result = decision_counts(make_results())
    assert list(result.keys()) == EXPECTED_DECISION_KEYS
    assert result == {"prioritize": 1, "review": 2, "exclude": 1}


def test_decision_counts_empty_has_zeros():
    result = decision_counts([])
    assert list(result.keys()) == EXPECTED_DECISION_KEYS
    assert result == {"prioritize": 0, "review": 0, "exclude": 0}


def test_impact_counts_ignores_excluded():
    rows = [
        {
            **make_row(variant="r1", gene="BRCA2", impact="high"),
            "decision": "prioritize",
            "reason": "priority gene brca2 with high impact",
        },
        {
            **make_row(variant="r2", gene="TP53", impact="modifier", depth=DEPTH_LOW),
            "decision": "exclude",
            "reason": "depth 10 below 12",
        },
        {
            **make_row(variant="r3", gene="EGFR", impact=" MODERATE "),
            "decision": "review",
            "reason": "non-priority gene egfr with moderate impact",
        },
    ]
    result = impact_counts(rows)
    assert list(result.keys()) == EXPECTED_IMPACT_KEYS
    assert result == {"high": 1, "moderate": 1, "low": 0, "modifier": 0}


def test_gene_summary_alphabetical_and_mean_depth():
    rows = [
        {
            **make_row(variant="r1", gene="TP53", impact="modifier", depth=DEPTH_LOW),
            "decision": "exclude",
            "reason": "depth 10 below 12",
        },
        {
            **make_row(variant="r2", gene="BRCA2", impact="high", depth=DEPTH_HIGH),
            "decision": "prioritize",
            "reason": "priority gene brca2 with high impact",
        },
        {
            **make_row(variant="r3", gene="BRCA2", impact="low", depth=DEPTH_MID),
            "decision": "review",
            "reason": "priority gene brca2 with low impact",
        },
    ]
    summary = gene_summary(rows)
    assert list(summary.keys()) == ["BRCA2", "TP53"]
    assert summary["BRCA2"] == {
        "total": 2,
        "prioritize": 1,
        "review": 1,
        "exclude": 0,
        "mean_depth": EXPECTED_MEAN_DEPTH_BRCA2,
    }
    assert summary["TP53"] == {
        "total": 1,
        "prioritize": 0,
        "review": 0,
        "exclude": 1,
        "mean_depth": EXPECTED_MEAN_DEPTH_TP53,
    }


def test_gene_summary_empty():
    assert gene_summary([]) == {}


# Test 5 and 6: rank_variants(), results_to_csv() - DASNEEM

def test_rank_variants_order_and_excludes():
    # TODO: prioritize before review, then impact, then higher quality, higher depth, smaller pos
    results = [
        {"variant": "s1", "gene": "BRCA2", "chrom": "13", "pos": "100", "depth": "30", "quality": "90", "impact": "high", "decision": "prioritize", "reason": "r1"},
        {"variant": "s2", "gene": "MYC", "chrom": "8", "pos": "200", "depth": "20", "quality": "95", "impact": "high", "decision": "prioritize", "reason": "r2"},
        {"variant": "s3", "gene": "CFTR", "chrom": "7", "pos": "150", "depth": "25", "quality": "90", "impact": "moderate", "decision": "prioritize", "reason": "r3"},
        {"variant": "s4", "gene": "TP53", "chrom": "17", "pos": "50", "depth": "40", "quality": "99", "impact": "high", "decision": "review", "reason": "r4"},
        {"variant": "s5", "gene": "RB1", "chrom": "13", "pos": "300", "depth": "18", "quality": "80", "impact": "low", "decision": "review", "reason": "r5"},
        {"variant": "s6", "gene": "BRCA2", "chrom": "13", "pos": "80", "depth": "10", "quality": "95", "impact": "high", "decision": "exclude", "reason": "r6"},
    ]
    ranked = rank_variants(results,5)

    names = []
    for row in ranked:
        names.append(row["variant"])

    assert names == ["s2", "s1", "s3", "s4", "s5"]
    assert "s6" not in names

def test_rank_variants_ties():
    # TODO: two variants equal on decision/impact/quality/depth: smaller pos comes first
    results = [
        {"variant": "A", 
                "gene": "BRCA2", 
                "chrom": "13", 
                "pos": "200", 
                "depth": "30", 
                "quality": "90", 
                "impact":"high", 
                "decision": "prioritize", 
                "reason":"rA"},

         {"variant": "B", 
                  "gene": "BRCA2", 
                  "chrom": "13", 
                  "pos": "100", 
                  "depth": "30", 
                  "quality": "90", 
                  "impact":"high", 
                  "decision": "prioritize", 
                  "reason":"rB"},
    ]

    ranked = rank_variants(results,2)

    assert ranked[0]["variant"] == "B"
    assert ranked[1]["variant"] == "A"

def test_rank_variants_k_zero_and_k_too_large():
    # TODO: k=0 gives []; k bigger than candidates gives all candidates
    results = [
        {"variant": "s1", "gene": "BRCA2", "chrom": "13", "pos": "100", "depth": "30", "quality": "90", "impact": "high", "decision": "prioritize", "reason": "r1"},
        {"variant": "s2", "gene": "MYC", "chrom": "8", "pos": "200", "depth": "20", "quality": "95", "impact": "high", "decision": "prioritize", "reason": "r2"},
    ]

    #k=0 will return the empty list
    assert rank_variants(results, 0) == []

    #k=5 when only 2 qualify returns both in ranked order
    ranked_all = rank_variants(results,5)
    assert len(ranked_all) == 2
    assert ranked_all[0]["variant"] == "s2"
    assert ranked_all[1]["variant"] == "s1"

def test_rank_variants_negative_k_raises():
    with pytest.raises(ValueError):
        rank_variants([], -1)


def test_rank_variants_does_not_modify_input():
    # TODO: deepcopy before, compare after
    results = [
            {"variant": "s1", "gene": "BRCA2", "chrom": "13", "pos": "100", "depth": "30", "quality": "90", "impact": "high", "decision": "prioritize", "reason": "r1"},
            {"variant": "s2", "gene": "MYC", "chrom": "8", "pos": "200", "depth": "20", "quality": "95", "impact": "high", "decision": "prioritize", "reason": "r2"},
        ]
    original = deepcopy(results)
    rank_variants(results,5)
    assert results == original

def test_results_to_csv_format():
    # TODO: first line == "variant,gene,chrom,pos,depth,quality,impact,decision,reason"
    #       every line ends with "\n"; one line per result in the same order
    results = [
                {"variant": "s1", "gene": "BRCA2", "chrom": "13", "pos": "100", "depth": "30", "quality": "90", "impact": "high", "decision": "prioritize", "reason": "r1"},
                {"variant": "s2", "gene": "MYC", "chrom": "8", "pos": "200", "depth": "20", "quality": "95", "impact": "high", "decision": "prioritize", "reason": "r2"},
        ]

    csv_text = results_to_csv(results)
    lines = csv_text.splitlines()

    assert csv_text.endswith("\n")
    assert len(lines) == 3
    assert lines[0] == "variant,gene,chrom,pos,depth,quality,impact,decision,reason"
    assert lines[1] == "s1,BRCA2,13,100,30,90,high,prioritize,r1"
    assert lines[2] == "s2,MYC,8,200,20,95,high,prioritize,r2"

    #empty list returns only the header line and the trailing newline
    assert results_to_csv([]) == "variant,gene,chrom,pos,depth,quality,impact,decision,reason\n"

# Test 7 and 7: build_report()

def test_build_report_lines():
    # TODO: check "Thresholds: min_depth=... min_quality=..." line, "Priority genes: A, B, ...",
    #       "Variants: N", "prioritize: n", "review: n", "exclude: n"
    pass


def test_build_report_top_k_listing():
    # TODO: lines start with "1. v12 | GENE1 | high | prioritize" (check the separator in the brief)
    pass


def test_build_report_fewer_than_k_and_empty():
    # TODO: fewer than K qualifying variants lists only those; empty results do not crash
    pass


# ==================================================================
# Integration (Member 1): triage output feeds the other functions
# ==================================================================
def test_integration_triage_into_summaries():
    rows = [make_row(variant="a"), make_row(variant="b", depth=1),
            make_row(variant="c", gene=OTHER_GENE, impact="low")]
    results = triage(rows, TEST_PARAMS)
    counts = decision_counts(results)
    assert sum(counts.values()) == len(rows)

    impacts = impact_counts(results)
    assert sum(impacts.values()) == counts["prioritize"] + counts["review"]
    assert impacts["high"] == 1
    assert impacts["moderate"] == 0
    assert impacts["low"] == 0
    assert impacts["modifier"] == 0

    summary = gene_summary(results)
    assert sum(gene_info["total"] for gene_info in summary.values()) == len(rows)
    assert list(summary.keys()) == sorted(summary.keys())
    # TODO: also call gene_summary, rank_variants and results_to_csv on `results`