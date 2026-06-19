from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

try:
    from scripts.common import ensure_output_dir, relative_or_absolute, resolve_project, write_csv, write_json
except ModuleNotFoundError:
    from common import ensure_output_dir, relative_or_absolute, resolve_project, write_csv, write_json


DATA_EXT = {".xlsx", ".xls", ".csv", ".tsv"}
MANUSCRIPT_EXT = {".docx", ".doc", ".pdf", ".md"}
FIGURE_EXT = {".pdf", ".svg", ".tif", ".tiff", ".png", ".jpg", ".jpeg"}
CODE_EXT = {".py", ".ipynb", ".r", ".rmd"}
VERSION_TOKENS = ("backup", "before", "old", "copy", "副本", "tmp", "restored", "updated_v")
JOURNAL_TOKENS = ("eja", "european_journal_of_agronomy", "journal_of_integrative_agriculture")


def classify(path: Path) -> str:
    lower = str(path).lower()
    suffix = path.suffix.lower()
    if suffix in DATA_EXT:
        return "data"
    if suffix in CODE_EXT:
        return "code"
    if suffix in FIGURE_EXT and any(x in lower for x in ("figure", "fig", "picture", "graphical")):
        return "figure"
    if suffix in MANUSCRIPT_EXT and any(x in lower for x in ("manuscript", "paper", "稿", "中文")):
        return "manuscript"
    if "cover" in lower and suffix in MANUSCRIPT_EXT:
        return "cover_letter"
    if "declaration" in lower or "competing" in lower:
        return "declaration"
    if "supp" in lower or "附表" in lower or "补充" in lower:
        return "supplementary"
    return "other"


def audit_project(project: Path) -> dict:
    project = project.resolve()
    records = []
    issues = []
    for path in sorted(p for p in project.rglob("*") if p.is_file()):
        category = classify(path)
        rel = relative_or_absolute(path, project)
        lower = path.name.lower()
        records.append({"path": rel, "category": category, "extension": path.suffix.lower(), "bytes": path.stat().st_size})
        if any(token in lower for token in VERSION_TOKENS):
            issues.append({"code": "version_or_backup_file", "path": rel, "message": "Version or backup marker in filename."})
        if any(token in lower for token in JOURNAL_TOKENS):
            issues.append({"code": "previous_journal_name", "path": rel, "message": "Possible previous-journal residue in filename."})
    counts_raw = Counter(r["category"] for r in records)
    manuscript_count = counts_raw["manuscript"]
    counts = {
        "all_files": len(records),
        "manuscript_files": manuscript_count,
        "data_files": counts_raw["data"],
        "code_files": counts_raw["code"],
        "figure_files": counts_raw["figure"],
        "cover_letters": counts_raw["cover_letter"],
        "declarations": counts_raw["declaration"],
        "supplementary_files": counts_raw["supplementary"],
    }
    if counts["data_files"] == 0:
        issues.append({"code": "missing_data", "path": "", "message": "No data file discovered."})
    if counts["code_files"] == 0:
        issues.append({"code": "missing_analysis_code", "path": "", "message": "No analysis code discovered."})
    return {
        "project": str(project),
        "counts": counts,
        "warnings": {"multiple_manuscripts": manuscript_count > 1},
        "issues": issues,
        "files": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Read-only inventory audit for an agricultural paper project.")
    parser.add_argument("--project")
    parser.add_argument("--output")
    args = parser.parse_args()
    default = Path.home() / "Desktop" / "lwxg"
    project = resolve_project(args.project, default)
    output = ensure_output_dir(project, args.output)
    report = audit_project(project)
    write_json(output / "project_audit.json", report)
    write_csv(output / "project_manifest.csv", report["files"], ["path", "category", "extension", "bytes"])
    print(output)


if __name__ == "__main__":
    main()
