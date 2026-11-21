"""Assignment 1 - clinical variant triage service (synthetic teaching data only).

Reads variant calls, applies ordered triage rules with named thresholds, records a
reason for every decision and produces summaries, a ranked top-k list, a CSV text
and a short text report. Nothing here is clinical advice.
"""
import csv
import json
import sys
from pathlib import Path

DECISIONS = ("prioritize", "review", "exclude")
IMPACT_ORDER = ("high", "moderate", "low", "modifier")
PRIORITY_IMPACTS = frozenset({"high", "moderate"})
CSV_COLUMNS = [
    "variant", "gene", "chrom", "pos", "depth", "quality", "impact", "decision", "reason",
]
CLI_ARGUMENTS = 3
MEAN_DEPTH_DECIMALS = 1


def load_variants(path) -> list[dict]:
    """Read variants.csv into a list of dicts. Every value stays text (Week 6 explains how)."""
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_params(path) -> dict:
    """Read the group parameters (data/params.json) into a dict."""
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def classify_variant(row: dict, params: dict) -> tuple[str, str]:
    """Return (decision, reason) for one variant row.

    Rules run in this order and the first rule that fires decides:
      1. depth < params["min_depth"]                      -> exclude (low depth)
      2. quality < params["min_quality"]                  -> exclude (low quality)
      3. priority gene and impact high or moderate        -> prioritize
      4. priority gene and impact low or modifier         -> review
      5. other gene and impact high or moderate           -> review
      6. other gene and impact low or modifier            -> exclude
    Gene symbols and impact labels are compared after strip() and case folding.
    depth, quality and pos may arrive as text.
    """
    #Convert text to numbers and clean text
    depth = int(row["depth"])
    quality = int(row["quality"])
    gene = row["gene"].strip().casefold()
    impact = row["impact"].strip().casefold()

    #Clean priority list from params
    priority_genes = {
        gene_name.strip().casefold()
        for gene_name in params["priority_genes"]
    }

    #Rule 1: Low depth
    if depth < params["min_depth"]:
        return "exclude", f"depth {depth} below {params['min_depth']}"

    #Rule 2: Low quality
    if quality < params["min_quality"]:
        return "exclude", f"quality {quality} below {params['min_quality']}"

    #Check condition firsy for rule 3-6
    is_priority = gene in priority_genes
    is_high_mod = impact in ("high", "moderate")

    #Rule 3: Priority gene + high/moderate impact
    if is_priority and is_high_mod:
        return "prioritize", f"priority gene {gene} with {impact} impact"

    #Rule 4: Priority gene + low/modifier impact
    if is_priority and not is_high_mod:
        return "review", f"priority gene {gene} with {impact} impact"

    #Rule 5: Non-priority gene + high/moderate impact
    elif not is_priority and is_high_mod:
        return "review", f"non-priority gene {gene} with {impact} impact"

    #Rule 6: Non-priority gene + low/modifier impact
    else:
        return "exclude", f"non-priority gene {gene} with {impact} impact"


def triage(rows: list[dict], params: dict) -> list[dict]:
    """Classify every row. Returns NEW dicts (row fields plus decision and reason), same order.

    The input rows and params are never modified.
    """
    result = []

    for row in rows:
        #get decision and reason
        decision, reason = classify_variant(row, params)

        #copy of the row
        new_row = row.copy()

        new_row["decision"] = decision
        new_row["reason"] = reason

        result.append(new_row)

    return result


def decision_counts(results: list[dict]) -> dict[str, int]:
    """Number of results per decision; keys are always prioritize, review, exclude (in order)."""
    counts = {decision: 0 for decision in DECISIONS}
    for row in results:
        decision = row["decision"].strip().casefold()
        if decision in counts:
            counts[decision] += 1
    return counts

def impact_counts(results: list[dict]) -> dict[str, int]:
    """Impact counts over variants that were NOT excluded.

    Keys are always high, moderate, low, modifier (in that order), zeros included.
    """
    counts = {impact: 0 for impact in IMPACT_ORDER}
    for row in results:
        if row["decision"].strip().casefold() == "exclude":
            continue
        impact = row["impact"].strip().casefold()
        if impact in counts:
            counts[impact] += 1
    return counts


def gene_summary(results: list[dict]) -> dict[str, dict]:
    """Per-gene summary, genes in alphabetical order.

    Each value is {"total": n, "prioritize": a, "review": b, "exclude": c, "mean_depth": d}
    where mean_depth is the mean depth of ALL the gene's variants, rounded to 1 decimal.
    """
    if not results:
        return {}

    grouped: dict[str, list[dict]] = {}
    for row in results:
        gene = row["gene"].strip()
        grouped.setdefault(gene, []).append(row)

    summary: dict[str, dict] = {}
    for gene in sorted(grouped):
        gene_rows = grouped[gene]
        per_decision = decision_counts(gene_rows)
        total_depth = sum(int(row["depth"]) for row in gene_rows)
        mean_depth = round(total_depth / len(gene_rows), MEAN_DEPTH_DECIMALS)
        summary[gene] = {
            "total": len(gene_rows),
            "prioritize": per_decision["prioritize"],
            "review": per_decision["review"],
            "exclude": per_decision["exclude"],
            "mean_depth": mean_depth,
        }

    return summary

def rank_variants(results: list[dict], k: int) -> list[dict]:
    """Top k variants that were not excluded, best first.

    Sort by: prioritize before review, then impact (high first), then higher quality,
    then higher depth, then smaller pos. Fewer than k are returned if fewer qualify.

    Raises:
        ValueError: if k is negative.
    """
    if k < 0 :
        raise ValueError("k must be a non negative number.")

    #this is to keep only the variant that is not excluded
    eligible = []
    for row in results:
                 decision = row["decision"].strip().casefold()
                 if decision != "exclude":
                     eligible.append(row)

    #to map the text label to the rank no. - small no = high priority
    decision_rank = {"prioritize":0, "review":1}
    impact_rank = {"high":0, "moderate":1, "low":2, "modifier":3}

    #the helper func will tell python how to compare 2 rows
    def sort_key(row:dict) -> tuple[int, int, int, int, int]:
        decision = row["decision"].strip().casefold()
        impact = row["impact"].strip().casefold()
        return(
            decision_rank[decision],
            impact_rank[impact],
            -int(row["quality"]), # -ve so big one come first
            -int(row["depth"]),    # -ve so big one come first
            int(row["pos"]),      # +ve so smaller pos comes first 
        )

    #sort & return the top k
    ranked = sorted(eligible, key=sort_key)
    return ranked[:k]        


def results_to_csv(results: list[dict]) -> str:
    """CSV text: a header line (CSV_COLUMNS) and one line per result in the given order.

    Every line ends with a newline. Reason text contains no commas, so no quoting is needed.
    """
    raise NotImplementedError("TODO: implement results_to_csv")


def build_report(results: list[dict], params: dict) -> str:
    """Short text report; the line formats are fixed so the report can be checked.

    Lines: title, 'Thresholds: min_depth=D min_quality=Q', 'Priority genes: A, B, ...',
    'Variants: N', one 'decision: count' line per decision, a blank line,
    'Top K variants' and then 'rank. variant | gene | impact | decision' lines.
    """
    raise NotImplementedError("TODO: implement build_report")


def main(argv: list[str] | None = None) -> int:
    """Command line: python src/a1_triage.py VARIANTS_CSV PARAMS_JSON OUT_DIR."""
    args = sys.argv[1:] if argv is None else argv
    if len(args) != CLI_ARGUMENTS:
        print("usage: python src/a1_triage.py data/variants.csv data/params.json results")
        return 2
    variants_path, params_path, out_dir = args
    rows = load_variants(variants_path)
    params = load_params(params_path)
    results = triage(rows, params)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "triage.csv").write_text(results_to_csv(results), encoding="utf-8")
    report = build_report(results, params)
    (out / "report.txt").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
