"""Shared helpers and the common output format used by all dataset modules.

Every dataset is written to ``data/structured/<dataset>/``:

    subjects.parquet   one row per subject      (SUBJECT_COLUMNS first, dataset-specific extras after)
    records.parquet    one row per recording    (RECORD_COLUMNS)
    bp.parquet         one row per BP value     (BP_COLUMNS); missing if the dataset has no BP values yet
    waves/<record_id>.h5   one HDF5 file per record, one dataset per signal (attrs: fs, unit, t0)
"""
import os
import shutil
from collections.abc import Iterable
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------------------------------------------------
# Paths. The data folder can be moved by setting the environment variable PPG_BP_DATA_DIR.
# ---------------------------------------------------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("PPG_BP_DATA_DIR", PROJECT_DIR / "data"))
RAW_DIR = DATA_DIR / "raw"
STRUCTURED_DIR = DATA_DIR / "structured"

# ---------------------------------------------------------------------------------------------------------------------
# Common schema
# ---------------------------------------------------------------------------------------------------------------------
SUBJECT_COLUMNS = {           # column -> dtype
    "subject_id": "string",
    "dataset": "string",
    "synthetic": "bool",      # True = simulated / virtual subject (e.g. pwdb), False = real person
    "age": "float64",         # years
    "sex": "string",          # "male" / "female" / <NA>
    "height_cm": "float64",
    "weight_kg": "float64",
}
RECORD_COLUMNS = {
    "record_id": "string",
    "subject_id": "string",
    "dataset": "string",
    "synthetic": "bool",      # True = simulated waves, False = measured; same value as the record's subject
    "session": "string",      # e.g. measurement day; <NA> if the dataset has only one session per subject
    "condition": "string",    # "rest", "sit", "walk", "icu", "surgery", ...
    "duration_s": "float64",
    "signals": "object",      # list of signal names in the HDF5 file
    "sampling_rates_hz": "object",  # list, same order as `signals`
    "wave_path": "string",    # relative to the dataset folder
}
BP_COLUMNS = {
    "record_id": "string",
    "time_s": "float64",      # time on the record's time axis; NaN = value describes the whole record
    "sbp": "float64",         # mmHg
    "dbp": "float64",
    "map": "float64",
    "method": "string",       # what measured the pressure: "cuff", "abp" (arterial line wave), "finapres"
    "source": "string",       # "dataset" = value given by the dataset, "calculated" = calculated by us from a wave
}
BP_SOURCES = {"dataset", "calculated"}


def normalize_sex(values: pd.Series) -> pd.Series:
    """Map different spellings (M/F, Male/Female, ...) to 'male' / 'female'. Unknown values become <NA>."""
    mapping = {"m": "male", "male": "male", "f": "female", "female": "female"}
    return values.astype(str).str.strip().str.lower().map(mapping).astype("string")


def mean_arterial_pressure(sbp: pd.Series, dbp: pd.Series) -> pd.Series:
    """Estimate mean arterial pressure from systolic and diastolic BP in mmHg: MAP = DBP + (SBP - DBP) / 3."""
    return dbp + (sbp - dbp) / 3


def _conform(df: pd.DataFrame, columns: dict[str, str], name: str, extra_allowed: bool) -> pd.DataFrame:
    """Add missing core columns (as NA), cast them to the shared dtypes and put them first."""
    df = df.copy()
    for col, dtype in columns.items():
        if col not in df:
            df[col] = pd.Series(None, index=df.index, dtype=dtype)
        elif dtype != "object":
            df[col] = df[col].astype(dtype)
    extra = [c for c in df.columns if c not in columns]
    if extra and not extra_allowed:
        raise ValueError(f"{name}: unexpected columns {extra}")
    return df[list(columns) + extra].reset_index(drop=True)


def prepare_out_dir(out_dir: Path) -> None:
    """Create the output folder and delete waves from an earlier run, so no stale files are left behind."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if (out_dir / "waves").is_dir():
        shutil.rmtree(out_dir / "waves")


def save_tables(out_dir: Path, subjects: pd.DataFrame, records: pd.DataFrame,
                bp: pd.DataFrame | Iterable[pd.DataFrame] | None, *, synthetic: bool) -> None:
    """Check the three tables against the shared schema and write them as parquet files.
    bp: one DataFrame, or an iterable of DataFrame chunks for tables too large for memory (written chunk by chunk,
    not sorted). Pass bp=None if the dataset has no BP values (yet); an old bp.parquet is then removed.
    synthetic: whether the dataset is simulated (written to the `synthetic` column of subjects and records)."""
    subjects = subjects.assign(synthetic=synthetic)
    records = records.assign(synthetic=synthetic)
    subjects = _conform(subjects, SUBJECT_COLUMNS, "subjects", extra_allowed=True)
    records = _conform(records, RECORD_COLUMNS, "records", extra_allowed=False)

    assert subjects["subject_id"].is_unique, "duplicate subject_id"
    assert records["record_id"].is_unique, "duplicate record_id"
    unknown = set(records["subject_id"]) - set(subjects["subject_id"])
    assert not unknown, f"records refer to unknown subjects: {sorted(unknown)[:5]}"
    missing = [p for p in records["wave_path"] if not (out_dir / p).is_file()]
    assert not missing, f"{len(missing)} wave files missing, e.g. {missing[0]}"

    out_dir.mkdir(parents=True, exist_ok=True)
    subjects.to_parquet(out_dir / "subjects.parquet", index=False)
    records.to_parquet(out_dir / "records.parquet", index=False)

    bp_path = out_dir / "bp.parquet"
    bp_path.unlink(missing_ok=True)
    n_bp = None
    if isinstance(bp, pd.DataFrame):
        bp = _check_bp(bp, set(records["record_id"]))
        bp.sort_values(["record_id", "source", "time_s"]).to_parquet(bp_path, index=False)
        n_bp = len(bp)
    elif bp is not None:
        n_bp = _write_bp_chunks(bp_path, bp, set(records["record_id"]))

    print(f"{out_dir.name}{' (synthetic)' if synthetic else ''}: {len(subjects)} subjects, {len(records)} records, "
          f"{'no' if n_bp is None else n_bp} BP values -> {out_dir}")


def _check_bp(bp: pd.DataFrame, record_ids: set[str]) -> pd.DataFrame:
    bp = _conform(bp, BP_COLUMNS, "bp", extra_allowed=False)
    unknown = set(bp["record_id"]) - record_ids
    assert not unknown, f"bp refers to unknown records: {sorted(unknown)[:5]}"
    bad = set(bp["source"].fillna("<NA>")) - BP_SOURCES
    assert not bad, f"bp: invalid source values {bad}, allowed: {BP_SOURCES}"
    return bp


def _write_bp_chunks(path: Path, chunks: Iterable[pd.DataFrame], record_ids: set[str]) -> int:
    """Write bp chunk by chunk, so the whole table never has to be in memory. Returns the number of rows."""
    schema = pa.Schema.from_pandas(_check_bp(pd.DataFrame(columns=list(BP_COLUMNS)), record_ids), preserve_index=False)
    n = 0
    with pq.ParquetWriter(path, schema) as writer:
        for chunk in chunks:
            chunk = _check_bp(chunk, record_ids)
            writer.write_table(pa.Table.from_pandas(chunk, schema=schema, preserve_index=False))
            n += len(chunk)
    return n


# ---------------------------------------------------------------------------------------------------------------------
# Waves
# ---------------------------------------------------------------------------------------------------------------------
def write_wave(path: Path, signals: dict[str, tuple[np.ndarray, float, str, float]]) -> None:
    """Write one HDF5 file. signals: name -> (data, fs_hz, unit, t0_s).
    Sample i of a signal is at time t0 + i / fs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        for name, (data, fs, unit, t0) in signals.items():
            ds = f.create_dataset(name, data=np.asarray(data, dtype="float32"), compression="gzip", chunks=True)
            ds.attrs["fs"] = float(fs)
            ds.attrs["unit"] = unit
            ds.attrs["t0"] = float(t0)


def read_wave(path: str | Path, channels: list[str] | None = None,
              start_s: float | None = None, end_s: float | None = None) -> dict[str, pd.Series]:
    """Read channels (default: all), optionally only the time slice [start_s, end_s).
    Each Series is indexed by time in s."""
    out = {}
    with h5py.File(path, "r") as f:
        for name in channels or list(f.keys()):
            ds = f[name]
            fs, t0 = float(ds.attrs["fs"]), float(ds.attrs["t0"])  # type: ignore
            n = len(ds)  # type: ignore
            i0 = min(n, max(0, int(np.ceil((start_s - t0) * fs)))) if start_s is not None else 0
            i1 = min(n, max(i0, int(np.ceil((end_s - t0) * fs)))) if end_s is not None else n
            values = ds[i0:i1]  # type: ignore
            index = pd.Index(t0 + np.arange(i0, i1) / fs, name="time_s")
            out[name] = pd.Series(values, index=index, name=name)  # type: ignore
    return out
