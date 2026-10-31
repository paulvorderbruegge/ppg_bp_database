"""PTT-PPG dataset (Mehrgardt et al., PhysioNet v1.1.0): conversion to the shared format.

Raw:     data/raw/pttppg/{<record>.hea/.dat (WFDB), csv/subjects_info.csv}
Records: one per activity (sit / walk / run) per subject, 500 Hz
Signals: ecg (a.u.) and the six PPG channels pleth_1..6, renamed to ppg_1..6
         (1-3: distal phalanx, 4-6: proximal phalanx; one channel per wavelength each. The README is contradictory
         about which number is red / infrared / green, so the original numbering is kept.)
BP:      cuff readings at the start and at the end of each record
"""
from pathlib import Path

import pandas as pd
import wfdb

from ppg_bp_database.utils import mean_arterial_pressure, normalize_sex, prepare_out_dir, save_tables, write_wave

DATASET = "pttppg"
SYNTHETIC = False  # measured data
SUBJECTS_INFO = "csv/subjects_info.csv"

# WFDB channel -> our signal name
CHANNELS = {"ecg": "ecg", **{f"pleth_{i}": f"ppg_{i}" for i in range(1, 7)}}
# WFDB unit -> our unit. The headers call the ECG "mV", but the values (~ 1e4) are raw ADC counts.
UNITS = {"NU": "a.u.", "mV": "a.u."}


def _read_info(raw_dir: Path) -> pd.DataFrame:
    """One row per record: activity, demographics and the BP readings."""
    return pd.read_csv(raw_dir / SUBJECTS_INFO)


def _subject(record: pd.Series) -> pd.Series:
    """'s1_walk' -> 'pttppg_s1'"""
    return DATASET + "_" + record.str.split("_").str[0]


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    info = _read_info(raw_dir)
    info["subject_id"] = _subject(info["record"])

    # demographics are repeated for every record of a subject: take the first
    info = info.drop_duplicates("subject_id")
    return pd.DataFrame({
        "subject_id": info["subject_id"],
        "dataset": DATASET,
        "age": info["age"],
        "sex": normalize_sex(info["gender"]),
        "height_cm": info["height"],
        "weight_kg": info["weight"],
    })


def write_records(raw_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Convert every WFDB record (ECG + PPG channels) to one HDF5 file and return the records table."""
    rows = []
    for row in _read_info(raw_dir).itertuples():
        rec = wfdb.rdrecord(str(raw_dir / row.record), channel_names=list(CHANNELS))
        assert rec.sig_name == list(CHANNELS), f"{row.record}: channels {rec.sig_name}"

        signals = {CHANNELS[ch]: (rec.p_signal[:, i], rec.fs, UNITS.get(rec.units[i], rec.units[i]), 0.0)  # type: ignore
                   for i, ch in enumerate(rec.sig_name)}

        record_id = f"{DATASET}_{row.record}"
        wave_path = f"waves/{record_id}.h5"
        write_wave(out_dir / wave_path, signals)

        rows.append({
            "record_id": record_id,
            "subject_id": f"{DATASET}_{str(row.record).split('_')[0]}",
            "dataset": DATASET,
            "condition": row.activity,
            "duration_s": rec.sig_len / rec.fs,  # type: ignore
            "signals": list(signals),
            "sampling_rates_hz": [float(rec.fs)] * len(signals),  # type: ignore
            "wave_path": wave_path,
        })
    return pd.DataFrame(rows)


def build_bp(raw_dir: Path, records: pd.DataFrame) -> pd.DataFrame:
    """Two cuff readings per record: at the start (t = 0) and at the end (t = duration)."""
    info = _read_info(raw_dir)
    record_id = DATASET + "_" + info["record"]
    durations = records.set_index("record_id")["duration_s"]

    start = pd.DataFrame({"record_id": record_id, "time_s": 0.0,
                          "sbp": info["bp_sys_start"], "dbp": info["bp_dia_start"]})
    end = pd.DataFrame({"record_id": record_id, "time_s": record_id.map(durations),
                        "sbp": info["bp_sys_end"], "dbp": info["bp_dia_end"]})

    bp = pd.concat([start, end], ignore_index=True)
    bp["map"] = mean_arterial_pressure(bp["sbp"], bp["dbp"])
    bp["method"] = "cuff"
    return bp


def convert(raw_dir: Path, out_dir: Path) -> None:
    prepare_out_dir(out_dir)
    records = write_records(raw_dir, out_dir)
    save_tables(out_dir, build_subjects(raw_dir), records, build_bp(raw_dir, records), synthetic=SYNTHETIC)
