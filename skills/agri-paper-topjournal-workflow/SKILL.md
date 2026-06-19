---
name: agri-paper-topjournal-workflow
description: >-
  End-to-end workflow for agricultural experimental papers: audit raw data and
  manuscripts, design defensible Python statistics, create publication-grade
  figures, prepare DOCX manuscripts and submission packages, and verify every
  result. Use whenever a user mentions crop/agronomy experiments, Excel data,
  significance letters, mixed models, manuscript figures, journal polishing,
  or submission readiness, even when they do not explicitly ask for this skill.
---

# Agricultural Paper Top-Journal Workflow

Use this skill for agricultural experimental papers from raw files through submission. The default template is `C:\Users\唐梓然\Desktop\lwxg`, but all scientific facts must be rediscovered from the active project.

## Coordinate the available skills

- Use `spreadsheets` for workbook inspection and structured spreadsheet edits.
- Use `nature-figure` for every scientific figure creation, revision or QA task. The default backend here is Python.
- Use `documents` for DOCX editing and render/inspect QA.
- Use `nature-polishing` for publication English.
- Use `humanizer-zh-academic` for Chinese academic de-template revision.
- Use `nature-citation` for citation addition and verification.
- Use `nature-data` for data-availability and repository planning.
- Use `systematic-debugging` when an audit, model, plot or document step fails.
- Use `verification-before-completion` before claiming readiness.

Read only the reference module needed for the current stage:

- Statistical design and model rules: [statistics.md](references/statistics.md)
- Figure rules: [figures.md](references/figures.md)
- Writing and DOCX rules: [writing-docx.md](references/writing-docx.md)
- Submission QA: [submission-qa.md](references/submission-qa.md)
- `lwxg` template map: [lwxg-template.md](references/lwxg-template.md)

## Non-negotiable defaults

1. Use Python unless the user explicitly selects another backend.
2. Inspect source files before proposing or making changes.
3. Preserve sources. Write to a clean versioned output directory unless replacement is explicitly requested.
4. Infer statistics from replicate-level experimental units, never from means alone.
5. Do not invent error bars, significance, methods, references, weather, sequences or data.
6. Do not claim causality, optimization or mechanism beyond the design and validation.
7. Keep confirmatory inference separate from exploratory analyses.
8. Verify outputs after every substantial stage.

## Stage 1: Inventory the project

If no path is supplied, use `C:\Users\唐梓然\Desktop\lwxg` as the example and ask before changing it.

Run the project audit before reanalysis:

```powershell
python scripts/audit_project.py --project "C:\path\to\project" --output "C:\path\to\audit"
```

Identify the authoritative manuscript, raw dataset, current analysis code, figure set, palette reference and target-journal files. Report duplicate versions, stale journal names, backups and missing sources. Do not delete them during an audit.

## Stage 2: Audit data and design

Extract factors, levels, seasons/years, cultivars/genotypes, blocks, experimental units, sampling hierarchy and timing from both manuscript and data.

```powershell
python scripts/audit_data.py --input "C:\path\to\data.xlsx" --output "C:\path\to\audit"
```

Confirm every sheet, variable, row, replicate and treatment combination was read. In a design with three field replicates, ensure each inferential combination has the intended three plot-level units and that subsamples are not treated as independent plots.

Stop and report a blocker when the experimental unit cannot be determined, raw replicates are absent, or the manuscript and data disagree materially.

## Stage 3: Design the statistical analysis

Read [statistics.md](references/statistics.md). Before fitting models, write a short analysis specification:

- primary responses and hypotheses;
- fixed and random effects;
- experimental unit and exact `n`;
- interaction hierarchy;
- multiplicity method;
- planned estimated means, contrasts, effect sizes and uncertainty;
- exploratory analyses and their limitations.

Use a design-correct mixed model for blocked, split-plot, multi-season or repeated-measure experiments. Use a compact-letter display only from the corresponding fitted model and comparison family.

Save analysis-ready data and every result table as CSV. Save a formatted XLSX only when the manuscript needs a styled table.

## Stage 4: Produce figures

Read [figures.md](references/figures.md) and invoke `nature-figure`.

Before plotting, state the figure's conclusion, evidence logic, source table, error-bar source, significance method, export dimensions and reviewer risk.

Historical defaults to preserve unless the user overrides them:

- Times New Roman;
- user-provided palette image, with documented same-family extensions only when necessary;
- black-edged points, bars, boxes and filled marks;
- real uncertainty bars;
- model-supported `a/b/c` labels in clear whitespace;
- separate panels plus combined figures;
- vector PDF and 600 dpi PNG, with SVG where useful.

Audit figures after export:

```powershell
python scripts/audit_figures.py --directory "C:\path\to\figures" --output "C:\path\to\audit"
```

Visually inspect every figure. Automated checks cannot prove that labels, legends or network leaders are correctly placed.

## Stage 5: Compare old and new analyses

Create a comparison table with original method, new method, assumptions, changed values/significance, conclusion impact, advantage, limitation and manuscript action. When a previous plot used incomplete data, identify exactly which variables/rows were omitted and whether the corrected result changes the scientific conclusion.

Do not preserve an old result merely to match the paper. Preserve the original method only when the user explicitly requires it and it remains scientifically defensible.

## Stage 6: Write and assemble the manuscript

Read [writing-docx.md](references/writing-docx.md). Keep the main story efficient: yield-aroma trade-off, canopy/light-use evidence, aroma-related metabolism, multivariate synthesis and bounded model interpretation.

Remove internal revision language, filenames, worksheet names and journal-quartile commentary. Use conservative mechanism language. Keep all values synchronized across abstract, Results, Discussion, Conclusions, tables and captions.

For DOCX work, edit a copy, maintain one authoritative submission version and run the document render gate when available.

Audit the manuscript:

```powershell
python scripts/audit_manuscript.py --input "C:\path\to\manuscript.docx" --figures "C:\path\to\figures" --output "C:\path\to\audit"
```

AI declarations are author-controlled. Never insert or remove one without explicit instruction.

## Stage 7: Prepare the submission package

Read [submission-qa.md](references/submission-qa.md). Browse the target journal's current official Guide for Authors before finalization.

Build a clean folder containing only intended upload files. Verify manuscript, highlights, graphical abstract, separate figures, supplementary material, declaration, cover letter, data/code and README as applicable. Generate a manifest, QA report and upload checklist. Validate ZIP entries after compression.

Keep journal-scope risk separate from file completeness. Never promise acceptance or a CAS quartile outcome.

## Output contract

Prefer this structure:

```text
<project>_topjournal_output/
├── analysis_ready_data.csv
├── results/
├── figures/panels/
├── figures/combined/
├── tables/
├── manuscript/
├── supplementary/
├── submission_package/
└── reports/
```

Reports should include project inventory, data/design audit, statistical audit, figure QA, manuscript QA and submission QA.

## Completion gate

Before calling the work complete, verify:

- source files are unchanged unless replacement was requested;
- replicate structure and model terms are documented;
- all requested variables were analyzed or explicitly excluded with reasons;
- error bars and significance labels are traceable to saved result tables;
- figure files open and match manuscript references;
- manuscript contains no internal residue or corrupted symbols;
- submission package contains one authoritative manuscript and no backup files;
- unresolved scientific, scope or visual-render risks are stated plainly.
