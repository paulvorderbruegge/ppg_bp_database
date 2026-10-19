"""PTT-PPG dataset (PhysioNet, v1.1.0): download and conversion to the shared format."""
from pathlib import Path

import pandas as pd
import wfdb
import h5py

DATASET = "pttppg"
BASE_URL = "https://physionet.org/files/pulse-transit-time-ppg/1.1.0/csv/"
SUBJECTS_INFO = "csv/subjects_info.csv"

#-----------------------------------------------------------------------------------------
# TEMPORARY: only use subjects 1..MAX_SUBJECT while the download is incomplete (None = all)
MAX_SUBJECT = 5

def _read_info(raw_dir: Path) -> pd.DataFrame:
    """Read subjects_info.csv, optionally restricted to the first MAX_SUBJECT subjects."""
    info = pd.read_csv(raw_dir / SUBJECTS_INFO)
    if MAX_SUBJECT is not None:
        subject_nr = info["record"].str.extract(r"^s(\d+)_")[0].astype(int)
        info = info[subject_nr <= MAX_SUBJECT]
    return info.reset_index(drop=True)
#-----------------------------------------------------------------------------------------

def build_subjects(raw_dir: Path) -> pd.DataFrame:
    #info = pd.read_csv(raw_dir / SUBJECTS_INFO)
    info = _read_info(raw_dir)
    
    # One row per subject: use the "walk" record (demographics are the same in all three)
    info = info[info["record"].str.endswith("_walk")].copy()

    subjects = pd.DataFrame({
        "subject_id": DATASET + "_" + info["record"].str.replace("_walk", "", regex=False),
        "dataset": DATASET,
        "age": info["age"],
        "sex": info["gender"],
        "height_cm": info["height"],
        "weight_kg": info["weight"],
    })
    return subjects.reset_index(drop=True)
    
def _filter_ppg_channels(h) -> list[str]:
    return [s for s in h.sig_name if s.startswith("pleth")]    

def build_records(raw_dir: Path) -> pd.DataFrame:
    # info = pd.read_csv(raw_dir / SUBJECTS_INFO)
    info = _read_info(raw_dir)

    # Read the headers to get length, sampling rate and channel names
    headers = [wfdb.rdheader(str(raw_dir / name)) for name in info["record"]]

    records = pd.DataFrame({
        "record_id": DATASET + "_" + info["record"],
        "subject_id": DATASET + "_" + info["record"].str.split("_").str[0],
        "dataset": DATASET,
        "condition": info["activity"],
        "duration_s": [h.sig_len / h.fs for h in headers], # type: ignore
        "sampling_rate_hz": [h.fs for h in headers],
        "signals": [_filter_ppg_channels(h) for h in headers],
        "wave_path": "waves/" + DATASET + "_" + info["record"] + ".h5",
    })
    return records.reset_index(drop=True)

def build_bp(raw_dir: Path) -> pd.DataFrame:
    # info = pd.read_csv(raw_dir / SUBJECTS_INFO)
    info = _read_info(raw_dir)
    record_id = DATASET + "_" + info["record"]

    # Duration of each record, needed as the time of the "end" measurement
    durations = build_records(raw_dir).set_index("record_id")["duration_s"]

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
    bp["map"] = bp["dbp"] + (bp["sbp"] - bp["dbp"]) / 3
    bp["method"] = "cuff"

    return bp.sort_values(["record_id", "time_s"]).reset_index(drop=True)

def write_waves(raw_dir: Path, out_dir: Path) -> None:
    """Write one HDF5 file per record containing only the PPG channels (float32)."""
    info = _read_info(raw_dir) 
    # info = pd.read_csv(raw_dir / SUBJECTS_INFO)
    records = build_records(raw_dir)

    (out_dir / "waves").mkdir(parents=True, exist_ok=True)

    for name, signals, wave_path in zip(info["record"], records["signals"], records["wave_path"]):
        print(name)
        rec = wfdb.rdrecord(str(raw_dir / name), channel_names=signals)

        with h5py.File(out_dir / wave_path, "w") as f:
            f.attrs["fs"] = rec.fs
            for i, channel in enumerate(rec.sig_name): # type: ignore
                ds = f.create_dataset(
                    channel,
                    data=rec.p_signal[:, i].astype("float32"), # type: ignore
                    compression="gzip",
                    chunks=True,
                )
                ds.attrs["unit"] = rec.units[i] # type: ignore

def _debug():
    raw = Path("data/raw/pttppg")

    for hea in sorted(raw.glob("*.hea")):
        name = hea.stem
        h = wfdb.rdheader(str(raw / name))
        dat = raw / f"{name}.dat"
        expected = h.sig_len * h.n_sig * 1.5 if h.fmt[0] == "212" else None
        actual = dat.stat().st_size if dat.exists() else None
        ok = expected is not None and actual == expected
        print(f"{name:10} fmt={h.fmt[0]} expected={expected} actual={actual} {'OK' if ok else '<-- CHECK'}")

def convert(raw_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    subjects = build_subjects(raw_dir)
    subjects.to_parquet(out_dir / "subjects.parquet", index=False)

    records = build_records(raw_dir)
    records.to_parquet(out_dir / "records.parquet", index=False)

    bp = build_bp(raw_dir)
    bp.to_parquet(out_dir / "bp.parquet", index=False)

    _debug()

    write_waves(raw_dir, out_dir)

    # TODO: wave files (one parquet per record)