from __future__ import annotations

import argparse
import csv
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

from PIL import Image, ImageStat

try:
    from scripts.common import write_csv, write_json
except ModuleNotFoundError:
    from common import write_csv, write_json


RASTER = {".png", ".jpg", ".jpeg", ".tif", ".tiff"}
VECTOR = {".pdf", ".svg"}


def audit_figures(directory: Path) -> dict:
    directory = directory.resolve()
    records, issues, numbers = [], [], []
    for path in sorted(p for p in directory.rglob("*") if p.is_file() and p.suffix.lower() in RASTER | VECTOR):
        rec = {"path": str(path.relative_to(directory)), "extension": path.suffix.lower(), "bytes": path.stat().st_size, "width": None, "height": None, "dpi_x": None, "dpi_y": None, "mode": None}
        match = re.search(r"(?:figure|fig)[_\- ]?(\d+)", path.stem, re.I)
        if match:
            numbers.append(int(match.group(1)))
        try:
            if path.suffix.lower() in RASTER:
                with Image.open(path) as im:
                    dpi = im.info.get("dpi", (None, None))
                    rec.update({"width": im.width, "height": im.height, "dpi_x": round(float(dpi[0]), 1) if dpi and dpi[0] else None, "dpi_y": round(float(dpi[1]), 1) if dpi and len(dpi) > 1 and dpi[1] else None, "mode": im.mode})
                    if not dpi or not dpi[0] or float(dpi[0]) < 300:
                        issues.append({"code": "low_dpi", "path": rec["path"], "message": "Raster DPI is missing or below 300."})
                    if im.width < 600 or im.height < 400:
                        issues.append({"code": "small_pixel_dimensions", "path": rec["path"], "message": "Raster dimensions may be too small for publication."})
                    extrema = ImageStat.Stat(im.convert("L")).extrema[0]
                    if extrema[1] - extrema[0] < 2:
                        issues.append({"code": "nearly_blank_image", "path": rec["path"], "message": "Image appears nearly blank."})
            elif path.suffix.lower() == ".svg":
                root = ET.parse(path).getroot()
                rec["width"], rec["height"] = root.get("width"), root.get("height")
        except Exception as exc:
            issues.append({"code": "unreadable_figure", "path": rec["path"], "message": str(exc)})
        if any(token in path.name.lower() for token in ("preview", "backup", "before", "old", "副本")):
            issues.append({"code": "submission_backup_or_preview", "path": rec["path"], "message": "Backup or preview naming in figure directory."})
        records.append(rec)
    unique = sorted(set(numbers))
    missing = [n for n in range(min(unique), max(unique) + 1) if n not in unique] if unique else []
    bases = Counter(Path(r["path"]).stem for r in records)
    duplicates = sorted(k for k, v in bases.items() if v > 1)
    metadata = directory / "figure_metadata.csv"
    metadata_status = "not_provided"
    if metadata.exists():
        with metadata.open(encoding="utf-8-sig") as handle:
            list(csv.DictReader(handle))
        metadata_status = "provided"
    return {"directory": str(directory), "files": records, "issues": issues, "sequence": {"numbers": unique, "missing_numbers": missing}, "duplicate_base_names": duplicates, "metadata_status": metadata_status}


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit journal figure dimensions, DPI and numbering.")
    parser.add_argument("--directory", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=True)
    report = audit_figures(Path(args.directory))
    write_json(output / "figure_audit.json", report)
    write_csv(output / "figure_qa.csv", report["files"], ["path", "extension", "bytes", "width", "height", "dpi_x", "dpi_y", "mode"])
    print(output)


if __name__ == "__main__":
    main()
