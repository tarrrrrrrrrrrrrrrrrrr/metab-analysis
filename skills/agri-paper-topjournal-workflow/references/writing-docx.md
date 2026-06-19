# Writing and DOCX Rules

## Scientific writing

- Use `nature-polishing` for English academic prose and `humanizer-zh-academic` for Chinese de-template revision.
- Use `nature-citation` when adding or verifying references. Verify title, journal, year, volume, pages/article number, DOI and claim fit.
- Results should answer what changed, by how much, uncertainty, significance and relevant interaction. Discussion explains mechanisms, limitations and literature context.
- Prefer “associated with”, “consistent with”, “may reflect” and “suggests possible involvement” unless manipulation and validation establish causality.
- Do not fabricate plot size, weather, dates, instruments, standards, internal standards, qPCR primers, gene accessions or formulas.

## Remove internal residue

Delete or flag manuscript-facing text such as “original manuscript”, “revised analysis”, “easy to reject”, “high-impact journal”, “CAS Q2”, `count.xlsx`, `Sheet1`, local paths, code filenames and workflow notes.

## Consistency

Abstract, Results, Discussion, Conclusions, captions, tables and supplementary files must use identical treatment names, stages, units, sample sizes and key values. Compare numeric claims programmatically where possible.

## DOCX editing

- Use the `documents` skill. Preserve the original and edit a copy unless the user explicitly requests replacement.
- Keep one authoritative manuscript in the clean submission folder.
- Insert figures in narrative order and keep captions with their figures.
- Use Word-safe scientific symbols: true superscript runs or verified Unicode glyphs for exponents, `μ`, `°C` and `×`.
- After edits, check for question-mark corruption and verify DOCX ZIP integrity.
- Render with `render_docx.py`, inspect page PNGs and iterate. If LibreOffice/soffice is unavailable, report that page-level QA was not completed.

## Supplement placement

Keep supplementary figures in the main text only when they materially advance the story or the user explicitly requests it. Avoid duplicate main/supplementary figures that make the same claim.

