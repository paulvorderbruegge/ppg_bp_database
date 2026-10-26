"""Graphene bioimpedance BP dataset (PhysioNet, v1.0.0): Conversion to the shared format."""
from pathlib import Path

import numpy as np
import pandas as pd

from ppg_bp_database.utils import write_wave

DATASET = "bpgraphene"

# file suffix -> (signal name, sampling rate in Hz, unit); bioz is deliberately skipped
SIGNALS = {
    "ppg": ("ppg", 1250.0, "a.u."),
    "finapresPPG": ("finapres_ppg", 75.0, "a.u."),
    "finapresBP": ("finapres_bp", 200.0, "mmHg"),
}


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    """One row per subject (7). The dataset gives no demographics."""
    numbers = sorted({int(p.name.removeprefix("subject").split("_")[0])
                      for p in raw_dir.glob("subject*_day*")})
    return pd.DataFrame({
        "subject_id": [f"{DATASET}_s{n}" for n in numbers],
        "dataset": DATASET,
        "age": np.nan,
        "sex": None,
        "height_cm": np.nan,
        "weight_kg": np.nan,
    })


def write_records(raw_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Convert every trial to one HDF5 file and return the records table."""
    rows = []
    for bp_file in sorted(raw_dir.glob("subject*_day*/setup*/data_trial*_finapresBP.csv")):
        experiment, setup = bp_file.parent.parent.name, bp_file.parent.name  # "subject1_day1", "setup02_hgcp"
        subject, day = experiment.removeprefix("subject").split("_day")
        trial = bp_file.stem.split("_")[1].removeprefix("trial")             # "01"

        record_id = f"{DATASET}_s{subject}_day{day}_trial{trial}"
        wave_path = f"waves/{record_id}.h5"

        signals = {}
        for suffix, (name, fs, unit) in SIGNALS.items():
            df = pd.read_csv(bp_file.with_name(f"data_trial{trial}_{suffix}.csv"))
            # first column = time in s, second column = the signal
            signals[name] = (df.iloc[:, 1].to_numpy(), fs, unit, float(df.iloc[0, 0]))
        write_wave(out_dir / wave_path, signals)

        rows.append({
            "record_id": record_id,
            "subject_id": f"{DATASET}_s{subject}",
            "dataset": DATASET,
            "session": f"day{day}",
            "condition": setup.split("_", 1)[1],                             # "hgcp", "valsalva", ...
            "duration_s": max(len(d) / fs for d, fs, _, _ in signals.values()),
            "signals": list(signals),
            "sampling_rates_hz": [fs for _, fs, _, _ in signals.values()],
            "wave_path": wave_path,
        })
    return pd.DataFrame(rows)


def convert(raw_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    build_subjects(raw_dir).to_parquet(out_dir / "subjects.parquet", index=False)
    write_records(raw_dir, out_dir).to_parquet(out_dir / "records.parquet", index=False)

    # TODO: derive BP values (beat-wise SBP/DBP/MAP) from the finapres_bp wave and write bp.parquet