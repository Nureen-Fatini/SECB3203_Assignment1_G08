# Design - Variant Triage Service

Write this BEFORE you write the code, then keep it up to date. It is marked by the lecturer.
Use your own words, your own thresholds and your own data. Delete the guidance in italics.

Group ID: <G08>      Thresholds from data/params.json: min_depth = <12>, min_quality = <75>, top_k = <5>

## 1. Decomposition

| Question | Our answer |
|---|---|
| Input | A `row` dictionary (variant, gene, depth, quality, impact as text) and a `params` dictionary (min_depth, min_quality, priority_genes, top_k). |
| Output | A tuple `(decision, reason)`. Decision is 'exclude', 'review' or 'prioritize'. Reason is a string. `triage` returns a new list of new dictionaries (row fields plus decision and reason) in the same order. |
| Constraints (named thresholds, what happens at exactly the threshold) | Named thresholds come from `params`: `min_depth` (12), `min_quality` (75) and `priority_genes`.None are typed into the rules. The rules are checked in order 1 to 6 and the first one that fires decides. Boundaries are inclusive: `depth == min_depth` and `quality == min_quality` pass and are NOT excluded. depth and quality are converted with `int()` because CSV values arrive as text. Gene and impact text is stripped and lowercased before comparing. Reasons contain no commas. The input rows and params are never modified. |
| Rule order and why the order matters | 1. depth < min_depth: exclude <br> 2. quality < min_quality: exclude <br> 3. priority gene + high/moderate: prioritize <br> 4. priority gene + low/modifier: review <br> 5. other gene + high/moderate: review <br> 6. other gene + low/modifier: exclude. <br> The two quality-control rules run first, so a variant that fails QC can never be prioritized or reviewed by its gene and impact. Depth is checked before quality, so a variant with depth 1 and quality 1 reports low depth. |
| Out of scope | Real clinical diagnosis or guideline-based variant classification, reading the CSV and JSON files (done by `main()`), and the summaries, ranking, CSV export and text report (done by other functions). |

## 2. Pseudocode or flowchart of classify_variant

```text
FUNCTION classify_variant(row, params):
    depth = int(row["depth"])
    quality = int(row["quality"])
    gene = casefold(strip(row["gene"]))
    impact = casefold(strip(row["impact"]))
    priority_genes = set(casefold(strip(g)) for g in params["priority_genes"])
    
    IF depth < params["min_depth"] THEN RETURN ("exclude", "depth {depth} below min_depth {params['min_depth']}")
    IF quality < params["min_quality"] THEN RETURN ("exclude", "quality {quality} below min_quality {params['min_quality']}")
    
    is_priority = gene IN priority_genes
    is_high_mod = impact IN ["high", "moderate"]
    
    IF is_priority AND is_high_mod THEN RETURN ("prioritize", "priority gene {gene} with {impact} impact")
    ELSE IF is_priority AND NOT is_high_mod THEN RETURN ("review", "priority gene {gene} with {impact} impact")
    ELSE IF NOT is_priority AND is_high_mod THEN RETURN ("review", "non-priority gene {gene} with {impact} impact")
    ELSE RETURN ("exclude", "non-priority gene {gene} with {impact} impact")
END FUNCTION
```
## 3. Hand trace

| Variant | Gene | Depth | Quality | Impact | Rule that fires | Decision | Reason |
|---|---|---|---|---|---|---|---|
| v1 | RB1 | 9 | 90 | high | 1 (Depth) | exclude | depth 9 below 12 |
| v2 | TP53 | 10 | 90 | modifier | 1 (Depth) | exclude | depth 10 below 12 |
| v3 | MLH1 | 11 | 90 | high | 1 (Depth) | exclude | depth 11 below 12 |
| v4 | CFTR | 30 | 79 | high | 3 (Priority + High/Moderate) | prioritize | priority gene cftr with high impact |
| v5 | MYC | 30 | 80 | moderate | 3 (Priority + High/Moderate) | prioritize | priority gene myc with moderate impact |
| v6 | MLH1 | 30 | 81 | modifier | 6 (Non-priority + Low/Modifier) | exclude | non-priority gene mlh1 with modifier impact |
| v9 | EGFR | 33 | 68 | low | 2 (Quality) | exclude | quality 68 below 75 |
| v16 | BRCA2 | 32 | 88 | low | 4 (Priority + Low/Modifier) | review | priority gene brca2 with low impact |


## 4. Test plan

*   **Normal:** A standard variant with a priority gene and high impact returns `prioritize`. The other gene and impact combinations return the decision of their rule: priority gene + low impact returns `review` (Rule 4), non-priority gene + high impact returns `review` (Rule 5), and non-priority gene + modifier impact returns `exclude` (Rule 6).
*   **Boundary (Exactly at threshold):** Depth = 12 and Quality = 75. Must NOT be excluded (boundaries are inclusive).
*   **Boundary (One below threshold):** Depth = 11 returns `exclude` (Rule 1); Quality = 74 returns `exclude` (Rule 2). The reason contains the word, the observed value and the threshold (for example `quality 74 below 75`).
*   **Boundary (One above threshold):** Depth = 13 and Quality = 76. Must pass the QC rules and evaluate the impact/gene rules.
*   **Rule Order:** A variant with Depth=5 and Quality=5 must fire Rule 1 (Depth) and NOT Rule 2 (Quality) in the reason string.
*   **Awkward input (Messy text):** Gene `" brca2 "` and impact `" HIGH "` must be stripped and lowercased to match correctly.
*   **Awkward input (Empty list):** Calling `triage([])` must return an empty list `[]` without crashing.
*   **Invalid input:** Calling `classify_variant()` with `depth = "abc"` must raise `ValueError`. Calling `rank_variants()` with a negative `k` must raise `ValueError`.

## 5. Work plan

**Order of execution:**
1. Fatini completes the core engine (`classify_variant`, `triage`) first, as the other functions depend on its output.
2. Once pushed,other member work in parallel on their specific functions.
3. All members write their own tests.

**Roles:**
*   **Fatini:** `classify_variant`, `triage`, and tests.
*   **Yaash:** `decision_counts`, `impact_counts`, `gene_summary`, and tests.
*   **Dasneem:** `rank_variants`, `results_to_csv`, and tests.
*   **Sarves:** `build_report`, runs `main()` integration, writes `README.md`, and tests.

**Review process:**
All code is pushed to GitHub. No one merges to the `main` branch without at least one other member reviewing the code via a Pull Request (or by running the code locally to confirm it passes `pytest` and `ruff`).
