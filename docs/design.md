# Design - Variant Triage Service

Write this BEFORE you write the code, then keep it up to date. It is marked by the lecturer.
Use your own words, your own thresholds and your own data. Delete the guidance in italics.

Group ID: <G08>      Thresholds from data/params.json: min_depth = <12>, min_quality = <75>, top_k = <5>

## 1. Decomposition

| Question | Our answer |
|---|---|
| Input | A `row` dictionary (containing variant, gene, depth, quality, impact as text) and a `params` dictionary (containing min_depth, min_quality, priority_genes, top_k). |
| Output | A tuple `(decision, reason)`. Decision is 'exclude', 'review', or 'prioritize'. Reason is a string. |
| Constraints (named thresholds, what happens at exactly the threshold) | The rules are checked in exact order (1 to 6). Boundaries are inclusive (e.g., `depth == min_depth` passes and is NOT excluded). Gene and impact text must be stripped of spaces and lowercased. Reasons cannot contain commas. |
| Rule order and why the order matters | 1. Depth <br> 2. Quality <br> 3. Priority + high/mod <br> 4. Priority + low/mod <br> 5. Non-priority + high/mod <br> 6. Non-priority + low/mod. <br>  |
| Out of scope |  Real clinical diagnosis, reading the CSV file (this is handled by `main()`), and generating the final report. |

## 2. Pseudocode or flowchart of classify_variant

```text
FUNCTION classify_variant(row, params):
    depth = int(row["depth"])
    quality = int(row["quality"])
    gene = lower(strip(row["gene"]))
    impact = lower(strip(row["impact"]))
    priority_genes = lower(strip(params["priority_genes"]))
    
    IF depth < params["min_depth"] THEN RETURN ("exclude", "depth {depth} below {min_depth}")
    IF quality < params["min_quality"] THEN RETURN ("exclude", "quality {quality} below {min_quality}")
    
    is_priority = gene IN priority_genes
    is_high_mod = impact IN ["high", "moderate"]
    
    IF is_priority AND is_high_mod THEN RETURN ("prioritize", "priority gene {gene} with {impact} impact")
    ELSE IF is_priority AND NOT is_high_mod THEN RETURN ("review", "priority gene {gene} with {impact} impact")
    ELSE IF NOT is_priority AND is_high_mod THEN RETURN ("review", "non-priority gene {gene} with {impact} impact")
    ELSE RETURN ("exclude", "non-priority gene {gene} with {impact} impact")
END FUNCTION

## 3. Hand trace

### 3. Hand trace

Using G08 parameters (`min_depth = 12`, `min_quality = 75`, `priority_genes = ["BRCA2", "RB1", "CFTR", "MYC", "KRAS"]`) on `data/variants.csv` (`v1`–`v6`, plus `v7` and `v13` to reach Rules 5 and 4)[cite: 1]:

| Variant | Gene | Depth | Quality | Impact | Rule that fires | Decision | Reason |
| --- | --- | --- | --- | --- | --- | --- | --- |
| v1 | RB1 | 9 | 90 | high | Rule 1 (`depth < min_depth`: 9 < 12)[cite: 1] | exclude[cite: 1] | Rule 1 depth 9 below min_depth 12[cite: 1] |
| v2 | TP53 | 10 | 90 | modifier | Rule 1 (`depth < min_depth`: 10 < 12)[cite: 1] | exclude[cite: 1] | Rule 1 depth 10 below min_depth 12[cite: 1] |
| v3 | MLH1 | 11 | 90 | high | Rule 1 (`depth < min_depth`: 11 < 12)[cite: 1] | exclude[cite: 1] | Rule 1 depth 11 below min_depth 12[cite: 1] |
| v4 | CFTR | 30 | 79 | high | Rule 3 (priority gene + high impact)[cite: 1] | prioritize[cite: 1] | Rule 3 priority gene cftr with high impact[cite: 1] |
| v5 | MYC | 30 | 80 | moderate | Rule 3 (priority gene + moderate impact)[cite: 1] | prioritize[cite: 1] | Rule 3 priority gene myc with moderate impact[cite: 1] |
| v6 | MLH1 | 30 | 81 | modifier | Rule 6 (non-priority gene + modifier impact)[cite: 1] | exclude[cite: 1] | Rule 6 non-priority gene mlh1 with modifier impact[cite: 1] |
| v7 | PTEN | 36 | 98 | high | Rule 5 (non-priority gene + high impact)[cite: 1] | review[cite: 1] | Rule 5 non-priority gene pten with high impact[cite: 1] |
| v13 | CFTR | 28 | 80 | modifier | Rule 4 (priority gene + modifier impact)[cite: 1] | review[cite: 1] | Rule 4 priority gene cftr with modifier impact[cite: 1] |


## 4. Test plan

*   **Normal:** A standard variant with a priority gene and high impact returns `prioritize`.
*   **Boundary (Exactly at threshold):** Depth = 12 and Quality = 75. Must NOT be excluded (boundaries are inclusive).
*   **Boundary (One below threshold):** Depth = 11 returns `exclude` (Rule 1); Quality = 74 returns `exclude` (Rule 2).
*   **Boundary (One above threshold):** Depth = 13 and Quality = 76. Must pass the QC rules and evaluate the impact/gene rules.
*   **Rule Order:** A variant with Depth=5 and Quality=5 must fire Rule 1 (Depth) and NOT Rule 2 (Quality) in the reason string.
*   **Awkward input (Messy text):** Gene `" brca2 "` and impact `" HIGH "` must be stripped and lowercased to match correctly.
*   **Awkward input (Empty list):** Calling `triage([])` must return an empty list `[]` without crashing.
*   **Invalid input:** Calling `rank_variants()` with a negative `k` must raise `ValueError`.

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
