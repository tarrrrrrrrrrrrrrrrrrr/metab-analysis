from __future__ import annotations

import argparse
import math
from itertools import product
from pathlib import Path

import pandas as pd

try:
    from scripts.common import write_csv, write_json
except ModuleNotFoundError:
    from common import write_csv, write_json


FACTOR_ALIASES = {
    "season": ("season", "季节", "year", "年份"),
    "cultivar": ("cultivar", "variety", "genotype", "material", "品种", "品系", "材料", "耕作模式"),
    "treatment": ("treatment", "shade", "shading", "nitrogen", "light", "处理", "遮光", "氮素", "光质"),
    "replicate": ("replicate", "rep", "block", "重复", "区组"),
}


def normalize(value: object) -> str:
    return str(value).strip().lower().replace(" ", "_")


def infer_column(columns: list[str], aliases: tuple[str, ...]) -> str | None:
    normalized = {normalize(c): c for c in columns}
    for alias in aliases:
        target = normalize(alias)
        if target in normalized:
            return normalized[target]
    return None


def audit_workbook(path: Path, factor_columns: list[str] | None = None, replicate_column: str | None = None) -> dict:
    path = path.resolve()
    if path.suffix.lower() in {".csv", ".tsv"}:
        frames = {path.stem: pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",")}
    else:
        with pd.ExcelFile(path) as book:
            frames = {name: pd.read_excel(book, sheet_name=name) for name in book.sheet_names}
    inventory = []
    for name, frame in frames.items():
        inventory.append({"sheet": name, "rows": int(frame.shape[0]), "columns": int(frame.shape[1]), "summary_only_candidate": bool(frame.shape[0] <= 10)})
    data_name, data = max(frames.items(), key=lambda item: item[1].shape[0] * max(item[1].shape[1], 1))
    columns = [str(c) for c in data.columns]
    if factor_columns is None:
        factor_columns = [c for key in ("season", "cultivar", "treatment") if (c := infer_column(columns, FACTOR_ALIASES[key]))]
    if replicate_column is None:
        replicate_column = infer_column(columns, FACTOR_ALIASES["replicate"])
    missingness = {str(c): int(data[c].isna().sum()) for c in data.columns}
    duplicate_count = int(data.duplicated().sum())
    factor_levels = {c: sorted(data[c].dropna().astype(str).unique().tolist()) for c in factor_columns if c in data.columns}
    replicate_levels = []
    if replicate_column and replicate_column in data.columns:
        replicate_levels = sorted(data[replicate_column].dropna().unique().tolist())
        replicate_levels = [int(x) if isinstance(x, float) and x.is_integer() else x for x in replicate_levels]
    expected_combinations = math.prod(len(v) for v in factor_levels.values()) if factor_levels else 0
    expected_units = expected_combinations * len(replicate_levels) if replicate_levels else expected_combinations
    observed_combinations = int(data[factor_columns].drop_duplicates().shape[0]) if factor_columns else 0
    numeric = data.select_dtypes(include="number")
    numeric_summary = []
    for col in numeric.columns:
        series = numeric[col].dropna()
        numeric_summary.append({"variable": str(col), "n": int(series.size), "min": float(series.min()) if len(series) else None, "max": float(series.max()) if len(series) else None, "mean": float(series.mean()) if len(series) else None})
    variables = [{"original_name": str(c), "normalized_name": normalize(c), "dtype": str(data[c].dtype), "missing": missingness[str(c)]} for c in data.columns]
    return {
        "input": str(path),
        "primary_sheet": data_name,
        "sheets": inventory,
        "design": {"factor_columns": factor_columns, "factor_levels": factor_levels, "replicate_column": replicate_column, "replicate_levels": replicate_levels},
        "duplicates": {"exact_rows": duplicate_count},
        "missingness": missingness,
        "coverage": {"expected_combinations": expected_combinations, "observed_combinations": observed_combinations, "expected_experimental_units": expected_units, "observed_rows": int(data.shape[0])},
        "variables": variables,
        "numeric_summary": numeric_summary,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit agricultural spreadsheet structure and replicate balance.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--factors", nargs="*")
    parser.add_argument("--replicate")
    args = parser.parse_args()
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=True)
    report = audit_workbook(Path(args.input), args.factors, args.replicate)
    write_json(output / "data_audit.json", report)
    write_csv(output / "sheet_inventory.csv", report["sheets"], ["sheet", "rows", "columns", "summary_only_candidate"])
    write_csv(output / "variable_dictionary.csv", report["variables"], ["original_name", "normalized_name", "dtype", "missing"])
    write_csv(output / "missingness.csv", [{"variable": k, "missing": v} for k, v in report["missingness"].items()], ["variable", "missing"])
    write_csv(output / "design_coverage.csv", [report["coverage"]], list(report["coverage"].keys()))
    print(output)


if __name__ == "__main__":
    main()
