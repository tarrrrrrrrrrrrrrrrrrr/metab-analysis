from __future__ import annotations

import csv
import json
from pathlib import Path


def resolve_project(path: str | None, default: Path) -> Path:
    project = Path(path).expanduser() if path else default
    project = project.resolve()
    if not project.is_dir():
        raise FileNotFoundError(f"Project directory not found: {project}")
    return project


def ensure_output_dir(project: Path, requested: str | None) -> Path:
    output = Path(requested).expanduser() if requested else project / "audit_reports"
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    return output


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())
