"""Load the converted datasets, e.g. in notebooks.

    from ppg_bp_database import access
    access.available_datasets()                    # ['ppgbp', 'pttppg', ...]
    records = access.load_records(["ppgbp", "pttppg"])   # one or several datasets, concatenated
    bp = access.load_bp("ppgbp")
    waves = access.load_wave(records.iloc[0])      # {signal name: pd.Series indexed by time in s}

    # harmonized data (scripts/harmonize.py)
    windows = access.load_windows()                # one row per window
    bp = access.load_window_bp()                   # BP rows per window (time_s relative to the window start)
    arrays = access.load_window_waves(windows.iloc[:32])   # {"ppg": (32, n), "abp": (32, n)}

Every table has the columns described in ppg_bp_database.utils (SUBJECT_COLUMNS, RECORD_COLUMNS, BP_COLUMNS).
"""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from ppg_bp_database.utils import BP_COLUMNS, HARMONIZED_DIR, STRUCTURED_DIR, read_wave

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
                                f"Run: uv run python scripts/structure.py {' '.join(missing)}")
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


# ---------------------------------------------------------------------------------------------------------------------
# Harmonized data
# ---------------------------------------------------------------------------------------------------------------------
def _harmonized(config: str, harmonized_dir: Path) -> Path:
    path = harmonized_dir / config
    if not (path / "windows.parquet").is_file():
        raise FileNotFoundError(f"no harmonized data for config {config!r} in {harmonized_dir}. "
                                f"Run: uv run python scripts/harmonize.py --config {config}")
    return path


def _dataset_filter(datasets: Datasets) -> list | None:
    if datasets is None:
        return None
    return [("dataset", "in", [datasets] if isinstance(datasets, str) else list(datasets))]


def load_windows(config: str = "default", datasets: Datasets = None,
                 harmonized_dir: Path = HARMONIZED_DIR) -> pd.DataFrame:
    """Harmonized windows (one row per window, no BP values)."""
    return pd.read_parquet(_harmonized(config, harmonized_dir) / "windows.parquet", filters=_dataset_filter(datasets))


def load_window_bp(config: str = "default", datasets: Datasets = None,
                   harmonized_dir: Path = HARMONIZED_DIR) -> pd.DataFrame:
    """BP rows of the harmonized windows: per beat where possible, else cuff / whole-record value.
    time_s is relative to the window start (NaN = describes the whole window). Join on window_id."""
    return pd.read_parquet(_harmonized(config, harmonized_dir) / "bp.parquet", filters=_dataset_filter(datasets))


def load_harmonized_subjects(config: str = "default", harmonized_dir: Path = HARMONIZED_DIR) -> pd.DataFrame:
    """Subjects that have at least one harmonized window."""
    return pd.read_parquet(_harmonized(config, harmonized_dir) / "subjects.parquet")


def load_window_waves(windows: pd.DataFrame, config: str = "default",
                      harmonized_dir: Path = HARMONIZED_DIR) -> dict[str, np.ndarray]:
    """{"ppg": (n, samples), "abp": (n, samples)} for the given rows of the windows table, in the same order.
    PPG is z-scored, ABP in mmHg; NaN after n_samples (short records) and for the whole abp row if has_abp is False."""
    path = _harmonized(config, harmonized_dir)
    out: dict[str, np.ndarray] = {}
    position = np.arange(len(windows))
    for dataset, rows in windows.reset_index(drop=True).groupby("dataset"):
        order = np.argsort(rows["row"].to_numpy())
        h5_rows = rows["row"].to_numpy()[order]
        target = position[rows.index.to_numpy()][order]
        with h5py.File(path / "waves" / f"{dataset}.h5", "r") as f:
            for name in ("ppg", "abp"):
                data = f[name][h5_rows]  # type: ignore  (h5py needs increasing indices)
                if name not in out:
                    out[name] = np.full((len(windows), data.shape[1]), np.nan, dtype=np.float32)
                out[name][target] = data
    return out
