"""Sanity check of the raw bp-graphene CSVs: columns, start time, measured vs. documented sampling rate."""
import sys
from pathlib import Path

import pandas as pd

from ppg_bp_database.datasets.bp_graphene_utils import SIGNALS

# optional argument: folder of a setup, e.g. data/raw/bpgraphene/subject1_day1/setup01_baseline
folder = Path(sys.argv[1] if len(sys.argv) > 1 else "data/raw/bpgraphene/subject1_day1/setup01_baseline")

for suffix, (name, fs_documented, unit) in SIGNALS.items():
    df = pd.read_csv(folder / f"data_trial01_{suffix}.csv")
    time = df.iloc[:, 0]
    fs_measured = 1 / time.diff().median()

    print(f"{name}")
    print(f"  columns:  {list(df.columns)}")
    print(f"  rows:     {len(df)}  ({len(df) / fs_documented:.1f} s at the documented rate)")
    print(f"  time:     {time.iloc[0]:.4f} s to {time.iloc[-1]:.4f} s")
    print(f"  rate:     documented {fs_documented} Hz, measured {fs_measured:.2f} Hz")