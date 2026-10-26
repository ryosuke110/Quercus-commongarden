#!/usr/bin/env python3
# Retrieve daily climate data and prepare Quercus occurrence-climate input
# Author: Ryosuke Ito
# AMGSDS data and AMD_Tools4.py must be obtained separately from the provider.

from __future__ import annotations

import argparse
import calendar
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


### Input and climate variables ###
# File paths, period and completeness thresholds are supplied through parse_args().
# Output name: (AMGSDS element, annual aggregation).
# All seven candidate variables are retained in both the intermediate files and
# the final input so that alternative variable selection and collinearity
# screening can be reproduced downstream.
ELEMENTS = {
    "AMT": ("TMP_mea", "mean"),
    "Tmax": ("TMP_max", "max"),
    "Tmin": ("TMP_min", "min"),
    "APCP": ("APCP", "sum"),
    "SD": ("SD", "max"),
    "SSD": ("SSD", "sum"),
    "RH": ("RH", "mean"),
}
FINAL_CLIMATE = ["AMT", "Tmax", "Tmin", "APCP", "SD", "SSD", "RH"]
FINAL_COLUMNS = [
    "plot_id", "Qmon", "Qser", "longitude", "latitude", *FINAL_CLIMATE
]


### Load data access tools ###
def load_amd_tools():
    try:
        import AMD_Tools4 as amd  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "AMD_Tools4.py was not found. Obtain the current IDCS-compatible "
            "version from the AMGSDS user documentation and place it beside "
            "this script."
        ) from exc
    return amd


### Read input arguments ###
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Qmon_Qser_occurrence.csv")
    parser.add_argument("--outdir", default="data_work")
    parser.add_argument("--start-year", type=int, default=2016)
    parser.add_argument("--end-year", type=int, default=2025)
    parser.add_argument(
        "--min-coverage", type=float, default=0.95,
        help="Minimum fraction of daily values required for an annual value",
    )
    parser.add_argument(
        "--min-valid-years", type=int, default=10,
        help="Minimum annual values required for the multi-year mean",
    )
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--delay", type=float, default=0.2)
    parser.add_argument("--limit", type=int, default=None, help="Test first N plots")
    parser.add_argument(
        "--no-resume", action="store_true",
        help="Delete the existing annual file and start from the first plot",
    )
    return parser.parse_args()


### Validate occurrence data ###
def validate_occurrence(plots: pd.DataFrame) -> pd.DataFrame:
    required = {"plot_id", "Qmon", "Qser", "longitude", "latitude"}
    missing = required.difference(plots.columns)
    if missing:
        raise SystemExit(f"Occurrence input is missing columns: {sorted(missing)}")
    plots = plots[["plot_id", "Qmon", "Qser", "longitude", "latitude"]].copy()
    if plots["plot_id"].isna().any() or plots["plot_id"].duplicated().any():
        raise SystemExit("plot_id must be non-missing and unique")
    for col in ["Qmon", "Qser"]:
        plots[col] = pd.to_numeric(plots[col], errors="raise").astype("int64")
        if not plots[col].isin([0, 1]).all():
            raise SystemExit(f"{col} must contain only 0 and 1")
    for col in ["longitude", "latitude"]:
        plots[col] = pd.to_numeric(plots[col], errors="coerce")
    return plots


### Retrieve daily climate data ###
def to_dates(times: object) -> pd.DatetimeIndex:
    raw = np.asarray(times).ravel()
    return pd.DatetimeIndex(
        pd.to_datetime([str(value)[:10] for value in raw], errors="raise")
    )


def fetch_one(amd, element: str, start: str, end: str, lat: float, lon: float,
              retries: int) -> tuple[np.ndarray, pd.DatetimeIndex, float, float]:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            data, times, grid_lat, grid_lon = amd.GetMetData(
                element, [start, end], [lat, lat, lon, lon]
            )
            values = np.atleast_1d(
                np.asarray(np.ma.filled(data, np.nan), dtype=float).squeeze()
            )
            dates = to_dates(times)
            if len(values) != len(dates):
                raise ValueError(
                    f"AMGSDS returned {len(values)} values but {len(dates)} dates"
                )
            return (
                values,
                dates,
                float(np.asarray(grid_lat).ravel()[0]),
                float(np.asarray(grid_lon).ravel()[0]),
            )
        except Exception as exc:  # transient network/session failures
            last_error = exc
            if attempt < retries:
                time.sleep(5 * attempt)
    assert last_error is not None
    raise last_error


### Aggregate daily data by year ###
def aggregate_years(values: np.ndarray, dates: pd.DatetimeIndex, method: str,
                    years: list[int], min_coverage: float) -> dict[int, tuple[float, float]]:
    series = pd.Series(values, index=dates)
    result: dict[int, tuple[float, float]] = {}
    for year in years:
        annual = series[series.index.year == year]
        expected = 366 if calendar.isleap(year) else 365
        coverage = float(np.isfinite(annual.to_numpy()).sum() / expected)
        if coverage < min_coverage:
            value = np.nan
        elif method == "mean":
            value = float(annual.mean(skipna=True))
        elif method == "max":
            value = float(annual.max(skipna=True))
        elif method == "min":
            value = float(annual.min(skipna=True))
        elif method == "sum":
            value = float(annual.sum(skipna=True))
        else:
            raise ValueError(f"Unknown aggregation method: {method}")
        result[year] = (value, coverage)
    return result


### Manage annual data and resume retrieval ###
def read_annual(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"plot_id": "string"})


def completed_plot_ids(path: Path, years: list[int]) -> set[str]:
    if not path.exists():
        return set()
    annual = read_annual(path)
    required_columns = {"plot_id", "year", *ELEMENTS.keys()}
    missing = required_columns.difference(annual.columns)
    if missing:
        raise SystemExit(
            f"Existing annual file has an incompatible schema: missing {sorted(missing)}. "
            "Use --no-resume to rebuild it."
        )
    expected = set(years)
    complete: set[str] = set()
    for plot_id, group in annual.groupby("plot_id", dropna=False):
        got = set(pd.to_numeric(group["year"], errors="coerce").dropna().astype(int))
        if got == expected and len(group) == len(years):
            complete.add(str(plot_id))
    return complete


def append_site(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(
        path, mode="a", header=not path.exists(), index=False, encoding="utf-8-sig"
    )


### Calculate multi-year means and save output ###
def make_outputs(plots: pd.DataFrame, annual_path: Path, outdir: Path,
                 min_valid_years: int, start_year: int, end_year: int) -> None:
    annual = read_annual(annual_path)
    if annual.duplicated(["plot_id", "year"]).any():
        raise SystemExit(
            "Duplicated plot_id-year rows were found in the annual file. "
            "Inspect it or rerun with --no-resume."
        )

    rows: list[dict] = []
    for plot_id, group in annual.groupby("plot_id", sort=False):
        row: dict[str, object] = {
            "plot_id": plot_id,
            "grid_latitude": (
                group["grid_latitude"].dropna().iloc[0]
                if group["grid_latitude"].notna().any() else np.nan
            ),
            "grid_longitude": (
                group["grid_longitude"].dropna().iloc[0]
                if group["grid_longitude"].notna().any() else np.nan
            ),
        }
        for col in ELEMENTS:
            values = pd.to_numeric(group[col], errors="coerce")
            n_valid = int(values.notna().sum())
            row[f"{col}_n_years"] = n_valid
            row[col] = float(values.mean()) if n_valid >= min_valid_years else np.nan
        rows.append(row)

    climate = pd.DataFrame(rows)
    climate_path = outdir / f"climate_{start_year}_{end_year}.csv"
    climate.to_csv(climate_path, index=False, encoding="utf-8-sig")

    merged = plots.merge(climate, on="plot_id", how="left", validate="one_to_one")
    invalid_coordinates = (
        merged["longitude"].isna() | merged["latitude"].isna()
        | ~merged["longitude"].between(122, 154)
        | ~merged["latitude"].between(20, 46.5)
    )
    missing_climate = merged[FINAL_CLIMATE].isna().any(axis=1)
    exclusion_reason = pd.Series(np.select(
        [invalid_coordinates & missing_climate, invalid_coordinates, missing_climate],
        ["invalid_coordinate_and_missing_climate", "invalid_coordinate", "missing_climate"],
        default="",
    ), index=merged.index)
    excluded = merged.loc[invalid_coordinates | missing_climate, FINAL_COLUMNS].copy()
    excluded["exclusion_reason"] = exclusion_reason.loc[excluded.index]
    final = merged.loc[~invalid_coordinates & ~missing_climate, FINAL_COLUMNS].copy()

    if final["plot_id"].duplicated().any():
        raise SystemExit("Internal validation failed: duplicated plot_id in final output")
    if not np.isfinite(final[FINAL_CLIMATE].to_numpy(dtype=float)).all():
        raise SystemExit("Internal validation failed: non-finite climate value")
    if list(final.columns) != FINAL_COLUMNS:
        raise SystemExit("Internal validation failed: final column order")

    final.to_csv(outdir / "Quercus_climate.csv", index=False, encoding="utf-8-sig")
    excluded.to_csv(outdir / "excluded_plots.csv", index=False, encoding="utf-8-sig")

    audit = pd.DataFrame(
        {
            "metric": [
                "input_plots", "output_plots", "excluded_plots",
                "Qmon_presences", "Qser_presences", "cooccurrence_plots",
                "start_year", "end_year", "min_valid_years",
            ],
            "value": [
                len(plots), len(final), len(excluded),
                int(final["Qmon"].sum()), int(final["Qser"].sum()),
                int(((final["Qmon"] == 1) & (final["Qser"] == 1)).sum()),
                start_year, end_year, min_valid_years,
            ],
        }
    )
    audit.to_csv(outdir / "Quercus_climate_audit.csv", index=False)


### Run climate retrieval and input preparation ###
def main() -> int:
    args = parse_args()
    if args.start_year > args.end_year:
        raise SystemExit("--start-year must not exceed --end-year")
    if not 0 < args.min_coverage <= 1:
        raise SystemExit("--min-coverage must be in (0, 1]")
    years = list(range(args.start_year, args.end_year + 1))
    if not 1 <= args.min_valid_years <= len(years):
        raise SystemExit(f"--min-valid-years must be between 1 and {len(years)}")

    # Prepare input and output paths
    input_path = Path(args.input).expanduser()
    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)
    annual_path = outdir / f"annual_climate_{args.start_year}_{args.end_year}.csv"
    settings_path = outdir / "aggregation_settings.json"

    plots = validate_occurrence(
        pd.read_csv(input_path, dtype={"plot_id": "string"})
    )
    if args.limit is not None:
        plots = plots.head(args.limit).copy()
    if args.no_resume and annual_path.exists():
        annual_path.unlink()

    # Record aggregation settings and check compatibility with saved data
    metadata = {
        "start_year": args.start_year,
        "end_year": args.end_year,
        "min_coverage": args.min_coverage,
        "min_valid_years": args.min_valid_years,
        "elements": {name: list(spec) for name, spec in ELEMENTS.items()},
        "final_columns": FINAL_COLUMNS,
    }
    if annual_path.exists() and settings_path.exists() and not args.no_resume:
        previous = json.loads(settings_path.read_text(encoding="utf-8"))
        retrieval_keys = ["start_year", "end_year", "min_coverage", "elements"]
        changed = [key for key in retrieval_keys if previous.get(key) != metadata.get(key)]
        if changed:
            raise SystemExit(
                "Retrieval settings differ from the existing annual file for: "
                f"{changed}. Restore the previous settings or use --no-resume."
            )
    settings_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # Identify plots not yet processed
    done = completed_plot_ids(annual_path, years)
    todo = plots.loc[~plots["plot_id"].isin(done)].copy()
    print(f"Plots={len(plots)}; completed={len(done)}; remaining={len(todo)}")

    # Retrieve daily data and append annual summaries for each plot
    amd = load_amd_tools()
    start, end = f"{args.start_year}-01-01", f"{args.end_year}-12-31"
    for index, row in enumerate(todo.itertuples(index=False), start=1):
        plot_id = str(row.plot_id)
        lat, lon = float(row.latitude), float(row.longitude)
        print(f"[{index}/{len(todo)}] plot_id={plot_id}", flush=True)
        if not np.isfinite(lat) or not np.isfinite(lon):
            rows = [
                {
                    "plot_id": plot_id, "year": year,
                    "source_latitude": lat, "source_longitude": lon,
                    "grid_latitude": np.nan, "grid_longitude": np.nan,
                    **{name: np.nan for name in ELEMENTS},
                    **{f"{name}_coverage": 0.0 for name in ELEMENTS},
                }
                for year in years
            ]
            append_site(annual_path, rows)
            continue

        rows = [
            {
                "plot_id": plot_id, "year": year,
                "source_latitude": lat, "source_longitude": lon,
                "grid_latitude": np.nan, "grid_longitude": np.nan,
            }
            for year in years
        ]
        try:
            for output_name, (element, method) in ELEMENTS.items():
                values, dates, grid_lat, grid_lon = fetch_one(
                    amd, element, start, end, lat, lon, args.retries
                )
                annual_values = aggregate_years(
                    values, dates, method, years, args.min_coverage
                )
                for item in rows:
                    value, coverage = annual_values[int(item["year"])]
                    item[output_name] = value
                    item[f"{output_name}_coverage"] = coverage
                    item["grid_latitude"] = grid_lat
                    item["grid_longitude"] = grid_lon
                time.sleep(args.delay)
        except Exception as exc:
            print(
                f"ERROR plot_id={plot_id}: {type(exc).__name__}: {exc}",
                file=sys.stderr,
            )
            print(
                "No partial rows were saved for this plot. Fix the cause and rerun "
                "the same command to resume.",
                file=sys.stderr,
            )
            return 1
        append_site(annual_path, rows)

    # Create the final occurrence-climate CSV and audit files
    if not annual_path.exists():
        raise SystemExit("No annual climate file was created")
    make_outputs(
        plots, annual_path, outdir, args.min_valid_years,
        args.start_year, args.end_year,
    )
    print(f"Created: {outdir / 'Quercus_climate.csv'}")
    print(f"Audit:   {outdir / 'Quercus_climate_audit.csv'}")
    return 0


### Execute pipeline ###
if __name__ == "__main__":
    raise SystemExit(main())

