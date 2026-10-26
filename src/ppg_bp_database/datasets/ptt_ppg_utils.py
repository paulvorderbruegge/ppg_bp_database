"""PTT-PPG dataset (PhysioNet, v1.1.0): Conversion to the shared format."""
from pathlib import Path

import pandas as pd
import wfdb

from ppg_bp_database.utils import normalize_sex, mean_arterial_pressure, write_wave

DATASET = "pttppg"
SUBJECTS_INFO = "csv/subjects_info.csv"


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    info = pd.read_csv(raw_dir / SUBJECTS_INFO)

    # One row per subject: use the "walk" record (demographics are the same in all three)
    info = info[info["record"].str.endswith("_walk")].copy()

    subjects = pd.DataFrame({
        "subject_id": DATASET + "_" + info["record"].str.replace("_walk", "", regex=False),
        "dataset": DATASET,
        "age": info["age"],
        "sex": normalize_sex(info["gender"]),
        "height_cm": info["height"],
        "weight_kg": info["weight"],
    })
    return subjects.reset_index(drop=True)


def _filter_ppg_channels(h) -> list[str]:
    return [s for s in h.sig_name if s.startswith("pleth")]


def build_records(raw_dir: Path) -> pd.DataFrame:
    info = pd.read_csv(raw_dir / SUBJECTS_INFO)

    # Read the headers to get length, sampling rate and channel names
    headers = [wfdb.rdheader(str(raw_dir / name)) for name in info["record"]]
    signals = [_filter_ppg_channels(h) for h in headers]                      # CHANGED

    records = pd.DataFrame({
        "record_id": DATASET + "_" + info["record"],
        "subject_id": DATASET + "_" + info["record"].str.split("_").str[0],
        "dataset": DATASET,
        "condition": info["activity"],
        "duration_s": [h.sig_len / h.fs for h in headers],  # type: ignore
        "signals": signals,
        "sampling_rates_hz": [[h.fs] * len(s) for h, s in zip(headers, signals)],
        "wave_path": "waves/" + DATASET + "_" + info["record"] + ".h5",
    })
    return records.reset_index(drop=True)


def build_bp(raw_dir: Path, records: pd.DataFrame) -> pd.DataFrame:
    info = pd.read_csv(raw_dir / SUBJECTS_INFO)
    record_id = DATASET + "_" + info["record"]

    # Duration of each record, needed as the time of the "end" measurement
    durations = records.set_index("record_id")["duration_s"]

    start = pd.DataFrame({
        "record_id": record_id,
        "time_s": 0.0,
        "sbp": info["bp_sys_start"],
        "dbp": info["bp_dia_start"],
    })
    end = pd.DataFrame({
        "record_id": record_id,
        "time_s": record_id.map(durations),
        "sbp": info["bp_sys_end"],
        "dbp": info["bp_dia_end"],
    })

    bp = pd.concat([start, end], ignore_index=True)
    bp["map"] = mean_arterial_pressure(bp["sbp"], bp["dbp"])
    bp["method"] = "cuff"

    return bp.sort_values(["record_id", "time_s"]).reset_index(drop=True)


def write_waves(raw_dir: Path, out_dir: Path, records: pd.DataFrame) -> None:
    for row in records.itertuples():
        name = row.record_id.removeprefix(DATASET + "_") # type: ignore
        rec = wfdb.rdrecord(str(raw_dir / name), channel_names=list(row.signals)) # type: ignore

        signals = {ch: (rec.p_signal[:, i], rec.fs, rec.units[i], 0.0)  # type: ignore
                   for i, ch in enumerate(rec.sig_name)}                 # type: ignore
        write_wave(out_dir / row.wave_path, signals) # type: ignore


def convert(raw_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    build_subjects(raw_dir).to_parquet(out_dir / "subjects.parquet", index=False)

    records = build_records(raw_dir)
    records.to_parquet(out_dir / "records.parquet", index=False)

    build_bp(raw_dir, records).to_parquet(out_dir / "bp.parquet", index=False)

    write_waves(raw_dir, out_dir, records)
