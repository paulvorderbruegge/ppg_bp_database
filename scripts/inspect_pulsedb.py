"""Run once on one downloaded PulseDB file to check the field names assumed in src/ppg_bp_database/datasets/pulsedb.py.
Usage: uv run python scripts/inspect_pulsedb.py data/raw/pulsedb/Segment_Files/PulseDB_MIMIC/<some file>.mat"""
import sys

import numpy as np
from mat73 import loadmat


def show(name, x, depth=0):
    pad = "  " * depth
    if isinstance(x, dict):
        print(f"{pad}{name}: dict with {len(x)} keys")
        for k, v in x.items():
            show(k, v, depth + 1)
    elif isinstance(x, list):
        first = np.asarray(x[0]) if len(x) else None
        print(f"{pad}{name}: list of {len(x)}, first: {None if first is None else (first.dtype, first.shape)}")
    else:
        a = np.asarray(x)
        print(f"{pad}{name}: array {a.dtype} {a.shape}  e.g. {a.reshape(-1)[:3]}")


data = loadmat(sys.argv[1])
show("file", data)

# T diagnostics: unit, sampling rate, clip starts
inner = next(iter(data.values())) if len(data) == 1 and isinstance(next(iter(data.values())), dict) else data
if "T" in inner:
    T = np.asarray(inner["T"], dtype=float).reshape(-1, 1250)
    print("\nT: segments:", len(T), "| dt of first segment:", np.diff(T[0])[:3], "| first starts:", T[:5, 0])
    print("number of times start time does not increase (= new clips):", int((np.diff(T[:, 0]) <= 0).sum()))
