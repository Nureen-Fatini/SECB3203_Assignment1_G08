# Variant Triage Service (SECB3203 Assignment 1)

> Replace the text in this template. The README must let a stranger install, run and test the
> project with the exact commands below. All data are synthetic teaching data, not clinical data.

## Group

- Group ID: <G08>
- Members: <full name (matric)>, <full name (matric)>, <full name (matric)>

## What this project does

<Three or four sentences: what a variant triage service is, what goes in, what comes out.>

## Set up

```
python -m venv .venv
<activate command for your system>
pip install -r requirements.txt
```

## Run

<Write the exact command. The course command is: python src/a1_triage.py data/variants.csv data/params.json results>
<Say which files appear in results/.>

## Test and lint

```
python -m pytest
ruff check src
```

## Our parameters

Our thresholds, priority genes and top_k are in data/params.json (<write the values here>).
They are read at run time and passed to the functions; they are not typed into the rules.

## Design

<Link to docs/design.md and say in two sentences what it contains.>

## Limitations

<At least two honest limitations, for example: rules use depth, quality and impact only; the data are synthetic.>

## Who did what

| Member | Main contributions | Commits |
|---|---|---|
| Nureen Fatini Binti Zulkefli | T1 Part 1; T2 – Implemented `classify_variant()` and `triage()`; T6 – Tests for classification and triage | [Count] |
| Yaaswiny Abbiramavally Thinaharan | T1 Part 3; T3 – Implemented the summary functions; T6 – Tests for summary functions | [Count] |
| Dasneem Banu Binti Haja | T1 Part 4 – Prepared the test plan; T4 – Implemented `rank_variants()`; T5 – Implemented `results_to_csv()`; T6 – Tests for variant ranking and CSV export | 8 |
| Sarveswari A/P Vasuthevan | T1 Part 2; T5 – Implemented `build_report()`; T6 – Tests for report generation | [Count] |

## AI use

The group used AI tools to support code development guidance, test planning hints and improvements to English wording. Details of AI usage are documented in `AI_USE_LOG.md`.
