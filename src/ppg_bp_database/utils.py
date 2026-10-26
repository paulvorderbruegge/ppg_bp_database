"""Shared helpers used by several dataset modules."""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def normalize_sex(values: pd.Series) -> pd.Series:
    """Map different spellings (M/F, Male/Female, ...) to 'male' / 'female'. Unknown values become NaN."""
    mapping = {"m": "male", "male": "male", "f": "female", "female": "female"}
    return values.astype(str).str.strip().str.lower().map(mapping)


def mean_arterial_pressure(sbp: pd.Series, dbp: pd.Series) -> pd.Series:
    """Estimate mean arterial pressure from systolic and diastolic BP in mmHg."""
    return dbp + (sbp - dbp) / 3


def write_wave(path: Path, signals: dict[str, tuple[np.ndarray, float, str, float]]) -> None:
    """Write one HDF5 file. signals: name -> (data, fs_hz, unit, t0_s).
    Sample i of a signal is at time t0 + i / fs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        for name, (data, fs, unit, t0) in signals.items():
            ds = f.create_dataset(name, data=data.astype("float32"), compression="gzip", chunks=True)
            ds.attrs["fs"] = fs
            ds.attrs["unit"] = unit
            ds.attrs["t0"] = t0


def read_wave(path: Path, channels: list[str] | None = None,
              start_s: float | None = None, end_s: float | None = None) -> dict[str, pd.Series]:
    """Read channels (default: all), optionally only a time slice. Each Series is indexed by time in s."""
    out = {}
    with h5py.File(path, "r") as f:
        for name in channels or list(f.keys()):
            ds = f[name]
            fs, t0 = float(ds.attrs["fs"]), float(ds.attrs["t0"])  # type: ignore
            i0 = max(0, int((start_s - t0) * fs)) if start_s is not None else 0
            i1 = int((end_s - t0) * fs) if end_s is not None else None
            values = ds[i0:i1]  # type: ignore
            out[name] = pd.Series(values, index=t0 + (i0 + np.arange(len(values))) / fs, name=name) # type: ignore
    return out
