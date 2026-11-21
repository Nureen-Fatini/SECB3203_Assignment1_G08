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

### 3. Hand trace

Using G08 parameters (`min_depth = 12`, `min_quality = 75`, `priority_genes = ["BRCA2", "RB1", "CFTR", "MYC", "KRAS"]`) on `data/variants.csv` (`v1`–`v6`, plus `v7` and `v13` to reach Rules 5 and 4)[cite: 1]:

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

**Group G08 Parameters (`data/params.json`):** `min_depth = 12`, `min_quality = 75`, `priority_genes = ["BRCA2", "RB1", "CFTR", "MYC", "KRAS"]`, `top_k = 5`.

| ID | Category | Function(s) | Input | Expected Output |
| :--- | :--- | :--- | :--- | :--- |
| **TC-01** | Normal | `classify_variant` | `BRCA2`, depth `30`, qual `90`, `high` | `("prioritize", "priority gene brca2 with high impact")` (Rule 3) |
| **TC-02** | Normal | `classify_variant` | `BRCA2`, depth `30`, qual `90`, `low` | `("review", "priority gene brca2 with low impact")` (Rule 4) |
| **TC-03** | Normal | `classify_variant` | `TP53`, depth `30`, qual `90`, `moderate` | `("review", "non-priority gene tp53 with moderate impact")` (Rule 5) |
| **TC-04** | Normal | `classify_variant` | `TP53`, depth `30`, qual `90`, `modifier` | `("exclude", "non-priority gene tp53 with modifier impact")` (Rule 6) |
| **TC-05** | Normal | `triage`, summaries, `rank_variants`, `results_to_csv` | **Shared Dataset (`s1`–`s6`)** below (`k=5`) | Input unchanged; `decision_counts`: `prioritize=3`, `review=2`, `exclude=1`; `impact_counts`: `high=3`, `moderate=1`, `low=1`, `modifier=0`; alphabetical `gene_summary`; top 5 `["s2", "s1", "s3", "s4", "s5"]`; CSV has the correct header, column order, input row order, and trailing newline |
| **TC-06** | Boundary | `classify_variant` | Tested separately (`BRCA2`, `high`): `depth="11"` vs `"12"` (qual `90`); `quality="74"` vs `"75"` (depth `30`) | Depth `11` -> `"exclude"`, `12` -> `"prioritize"`; Quality `74` -> `"exclude"`, `75` -> `"prioritize"` |
| **TC-07** | Awkward | `classify_variant` | `gene="  bRcA2 "`, `impact=" HIGH "`, depth `30`, qual `90` | `("prioritize", "priority gene brca2 with high impact")` |
| **TC-08** | Edge (Empty/Small) | `triage`, `decision_counts`, `impact_counts`, `gene_summary`, `rank_variants`, `results_to_csv` | Empty list `[]`; or 2 triaged results (`[s1, s4]`) with `k=5` | `triage([])` / `rank_variants([], 5)` -> `[]`; `decision_counts([])` / `impact_counts([])` -> all zeros; `gene_summary([])` -> `{}`; `results_to_csv([])` -> header only, ending with `\n`; `k=5` with `[s1, s4]` returns all 2 rows safely |
| **TC-09** | Invalid | `classify_variant`, `rank_variants` | Tested separately: `depth="abc"`, `quality="abc"`, `k=-1` | Raises `ValueError` for each separate invalid input |

### Shared Test Dataset (`s1`–`s6`) for TC-05 & TC-08
| Variant | Gene | Chrom | Pos | Depth | Quality | Impact | Expected Decision |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `s1` | BRCA2 | 13 | 100 | 30 | 90 | high | `prioritize` |
| `s2` | MYC | 8 | 200 | 20 | 95 | high | `prioritize` |
| `s3` | CFTR | 7 | 150 | 25 | 90 | moderate | `prioritize` |
| `s4` | TP53 | 17 | 50 | 40 | 99 | high | `review` |
| `s5` | RB1 | 13 | 300 | 18 | 80 | low | `review` |
| `s6` | BRCA2 | 13 | 80 | 10 | 95 | high | `exclude` (depth 10 < 12) |

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
