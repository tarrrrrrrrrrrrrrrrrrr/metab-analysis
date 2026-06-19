from __future__ import annotations

import argparse
import html
import re
from pathlib import Path
from zipfile import ZipFile

try:
    from scripts.common import write_json
except ModuleNotFoundError:
    from common import write_json


INTERNAL_PATTERNS = ("revised analysis", "revised manuscript", "original manuscript", "easy to reject", "修改稿", "原稿", "二区投稿")
SOURCE_PATTERNS = ("count.xlsx", "sheet1", "sheet 1", ".csv file", "source code file")
OLD_JOURNALS = ("european journal of agronomy", "journal of integrative agriculture", "cover_letter_eja")


def extract_docx_text(path: Path) -> tuple[str, int]:
    with ZipFile(path) as zf:
        bad = zf.testzip()
        if bad:
            raise ValueError(f"Corrupt DOCX member: {bad}")
        chunks, media = [], 0
        for name in zf.namelist():
            if name.startswith("word/media/"):
                media += 1
            if name.startswith("word/") and name.endswith(".xml"):
                xml = zf.read(name).decode("utf-8", errors="ignore")
                xml = re.sub(r"</w:p>", "\n", xml)
                chunks.append(html.unescape(re.sub(r"<[^>]+>", " ", xml)))
    text = "\n".join(chunks)
    text = re.sub(r"[ \t\r\f\v]+", " ", text)
    text = re.sub(r"\n+", "\n", text)
    return text.strip(), media


def audit_manuscript(path: Path, figure_dir: Path | None = None) -> dict:
    path = path.resolve()
    text, media_count = extract_docx_text(path)
    lower = text.lower()
    issues = []
    if any(token in lower for token in INTERNAL_PATTERNS):
        issues.append({"code": "internal_revision_language", "message": "Internal revision language detected."})
    if any(token in lower for token in SOURCE_PATTERNS):
        issues.append({"code": "source_filename_residue", "message": "Source filename or worksheet residue detected."})
    if any(token in lower for token in OLD_JOURNALS):
        issues.append({"code": "previous_journal_name", "message": "Possible previous-journal name detected."})
    if re.search(r"(?:\d|m|g|min|mol|l)\?{1,3}(?:\b|c|m|l)", lower) or "?c" in lower:
        issues.append({"code": "corrupted_scientific_symbol", "message": "Question-mark corruption near a scientific unit."})
    refs = sorted(set(map(int, re.findall(r"(?:Fig\.|Figure)\s*(\d+)", text, flags=re.I))))
    captions = sorted(set(map(int, re.findall(r"(?mi)^\s*(?:Fig\.|Figure)\s*(\d+)\s*\.", text))))
    external_numbers = []
    if figure_dir and figure_dir.exists():
        for file in figure_dir.rglob("*"):
            match = re.search(r"(?:Figure|Fig)[_\- ]?(\d+)", file.stem, re.I)
            if match:
                external_numbers.append(int(match.group(1)))
    declarations = {
        "credit": "credit authorship" in lower or "author contributions" in lower,
        "funding": "funding" in lower,
        "competing_interest": "competing interest" in lower or "conflict of interest" in lower,
        "data_availability": "data availability" in lower,
        "references": "references" in lower,
        "ai_declaration": "declaration of generative ai" in lower,
    }
    replicate_mentions = re.findall(r"(?:three|3)\s+(?:field\s+)?replicates?", lower)
    return {
        "input": str(path),
        "characters": len(text),
        "question_marks": text.count("?"),
        "embedded_media": media_count,
        "issues": issues,
        "figures": {"references": refs, "captions": captions, "references_without_captions": sorted(set(refs) - set(captions)), "external_files": sorted(set(external_numbers)), "missing_external_files": sorted(set(refs) - set(external_numbers)) if figure_dir else []},
        "declarations": declarations,
        "replicate_statement_count": len(replicate_mentions),
    }


def to_markdown(report: dict) -> str:
    lines = ["# Manuscript QA", "", f"- Input: `{report['input']}`", f"- Question marks: {report['question_marks']}", f"- Embedded media: {report['embedded_media']}", f"- Replicate statements: {report['replicate_statement_count']}", "", "## Issues"]
    lines.extend(f"- **{item['code']}**: {item['message']}" for item in report["issues"])
    if not report["issues"]:
        lines.append("- No rule-based issues detected.")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit DOCX manuscript text, figures and submission sections.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--figures")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=True)
    report = audit_manuscript(Path(args.input), Path(args.figures) if args.figures else None)
    write_json(output / "manuscript_audit.json", report)
    (output / "manuscript_qa.md").write_text(to_markdown(report), encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
