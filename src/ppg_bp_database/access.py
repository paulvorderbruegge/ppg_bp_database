"""Load the converted datasets, e.g. in notebooks.

    from ppg_bp_database import access
    access.available_datasets()                    # ['ppgbp', 'pttppg', ...]
    records = access.load_records(["ppgbp", "pttppg"])   # one or several datasets, concatenated
    bp = access.load_bp("ppgbp")
    waves = access.load_wave(records.iloc[0])      # {signal name: pd.Series indexed by time in s}

Every table has the columns described in ppg_bp_database.utils (SUBJECT_COLUMNS, RECORD_COLUMNS, BP_COLUMNS).
"""
from pathlib import Path

import pandas as pd

from ppg_bp_database.utils import BP_COLUMNS, STRUCTURED_DIR, read_wave

Datasets = str | list[str] | None  # one dataset, several, or None = all available


def available_datasets(structured_dir: Path = STRUCTURED_DIR) -> list[str]:
    """Datasets that have been converted (i.e. have a records.parquet)."""
    return sorted(p.parent.name for p in structured_dir.glob("*/records.parquet"))


def _names(datasets: Datasets, structured_dir: Path) -> list[str]:
    if datasets is None:
        return available_datasets(structured_dir)
    names = [datasets] if isinstance(datasets, str) else list(datasets)
    missing = [n for n in names if not (structured_dir / n / "records.parquet").is_file()]
    if missing:
        raise FileNotFoundError(f"not converted: {missing}. Available: {available_datasets(structured_dir)}. "
                                f"Run: uv run python scripts/convert.py {' '.join(missing)}")
    return names


def _load(table: str, datasets: Datasets, structured_dir: Path) -> pd.DataFrame:
    parts = []
    for name in _names(datasets, structured_dir):
        path = structured_dir / name / f"{table}.parquet"
        if path.is_file():
            parts.append(pd.read_parquet(path))
        elif table == "bp":  # dataset without BP values (yet)
            parts.append(pd.DataFrame({c: pd.Series(dtype=t) for c, t in BP_COLUMNS.items()}))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def load_subjects(datasets: Datasets = None, structured_dir: Path = STRUCTURED_DIR) -> pd.DataFrame:
    """Subjects table. Dataset-specific extra columns are NaN for the other datasets."""
    return _load("subjects", datasets, structured_dir)


def load_records(datasets: Datasets = None, structured_dir: Path = STRUCTURED_DIR) -> pd.DataFrame:
    """Records table. Filter on the `synthetic` column to separate measured from simulated data."""
    return _load("records", datasets, structured_dir)


def load_bp(datasets: Datasets = None, structured_dir: Path = STRUCTURED_DIR) -> pd.DataFrame:
    """BP values; empty for datasets that have none yet."""
    return _load("bp", datasets, structured_dir)


def wave_file(record: pd.Series, structured_dir: Path = STRUCTURED_DIR) -> Path:
    """Absolute path of a record's HDF5 file (record = one row of the records table)."""
    return structured_dir / record["dataset"] / record["wave_path"]


def load_wave(record: pd.Series, channels: list[str] | None = None, start_s: float | None = None,
              end_s: float | None = None, structured_dir: Path = STRUCTURED_DIR) -> dict[str, pd.Series]:
    """Signals of one record (one row of the records table), optionally only some channels or a time slice."""
    return read_wave(wave_file(record, structured_dir), channels, start_s, end_s)


def summary(structured_dir: Path = STRUCTURED_DIR) -> pd.DataFrame:
    """One row per converted dataset: synthetic or measured, counts, conditions and signals."""
    rows = []
    for name in available_datasets(structured_dir):
        records, bp = load_records(name, structured_dir), load_bp(name, structured_dir)
        rows.append({
            "dataset": name,
            "synthetic": bool(records["synthetic"].all()),
            "subjects": len(load_subjects(name, structured_dir)),
            "records": len(records),
            "bp_values": len(bp),
            "hours": round(records["duration_s"].sum() / 3600, 2),
            "conditions": ", ".join(sorted(records["condition"].dropna().unique())),
            "signals": ", ".join(sorted({s for sig in records["signals"] for s in sig})),
        })
    return pd.DataFrame(rows)
