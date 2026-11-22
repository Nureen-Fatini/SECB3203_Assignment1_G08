# Variant Triage Service (SECB3203 Assignment 1)

> Replace the text in this template. The README must let a stranger install, run and test the
> project with the exact commands below. All data are synthetic teaching data, not clinical data.

## Group

- Student name: Nureen Fatini Binti Zulkefli, Yaaswiny Abbiramavally Thinaharan, Dasneem Banu Binti Haja, Sarveswari A/P Vasuthevan
- Matric no: A24CS0169, A24CS0213, A24CS0066, A24CS0185
- Group ID: G08
- Members: Nureen Fatini Binti Zulkefli (A24CS0169), Yaaswiny Abbiramavally Thinaharan (A24CS0213), Dasneem Banu Binti Haja (A24CS0066), Sarveswari A/P Vasuthevan (A24CS0185)

## What this project does

<A variant triage service applies ordered, threshold-based rules to filter and prioritize genetic variant calls consistently instead of relying on ad-hoc spreadsheets. It takes a CSV file of variant records (such as gene, chromosome, position, depth, quality, and impact) and a JSON parameter configuration as inputs. From those inputs, it produces triage decisions (`prioritize`, `review`, `exclude`) with recorded reasons, summary counts, ranked top-k variants, a CSV export, and a formatted report.>

## Set up

```
python -m venv .venv
<activate command for your system>
pip install -r requirements.txt
```

## Run

<Run the pipeline using `python src/a1_triage.py data/variants.csv data/params.json results`.>
<Running the command generates `results/triage.csv` and `results/report.txt` inside the `results/` folder.>

## Test and lint

```
python -m pytest
ruff check src
```

## Our parameters

Our thresholds, priority genes and top_k are in data/params.json (`min_depth=12`, `min_quality=75`, `priority_genes=["BRCA2", "RB1", "CFTR", "MYC", "KRAS"]`, `top_k=5`).
They are read at run time and passed to the functions; they are not typed into the rules.

## Design

[docs/design.md](docs/design.md) contains the problem decomposition table, six-rule variant classification pseudocode/flowchart, and hand-trace tables for sample rows. It also outlines the group test plan, work plan, and design review notes.

## Limitations

This triage logic relies only on read depth, base/variant quality, priority gene membership, and impact labels under simplified teaching rules rather than clinical variant-interpretation guidelines. In addition, `data/variants.csv` uses synthetic teaching data and is not intended for clinical use or real patient records.

## Who did what

| Member | Main contributions | Commits |
|---|---|---|
| Nureen Fatini Binti Zulkefli | T1 Part 1; T2 – Implemented `classify_variant()` and `triage()`; T6 – Tests for classification and triage | 5 |
| Yaaswiny Abbiramavally Thinaharan | T1 Part 3; T3 – Implemented the summary functions; T6 – Tests for summary functions | 5 |
| Dasneem Banu Binti Haja | T1 Part 4 – Prepared the test plan; T4 – Implemented `rank_variants()`; T5 – Implemented `results_to_csv()`; T6 – Tests for variant ranking and CSV export | 8 |
| Sarveswari A/P Vasuthevan | T1 Part 2; T5 – Implemented `build_report()`; T6 – Tests for report generation | [4] |

## AI use

The group used AI tools to support code development guidance, test planning hints and improvements to English wording. Details of AI usage are documented in `AI_USE_LOG.md`.
