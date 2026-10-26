"""PPG-BP dataset (Liang et al. 2018): Conversion to the shared format."""
from pathlib import Path

import numpy as np
import pandas as pd

from ppg_bp_database.utils import mean_arterial_pressure, normalize_sex, write_wave

DATASET = "ppgbp"
SAMPLE_RATE_HZ = 1000.0  # 2100 samples per segment
EXCEL_FILE = "PPG-BP Database/Data File/PPG-BP dataset.xlsx"
WAVE_DIR = "PPG-BP Database/Data File/0_subject"

# Excel column -> our name
COLUMNS = {
    "subject_ID": "subject",
    "Sex(M/F)": "sex",
    "Age(year)": "age",
    "Height(cm)": "height_cm",
    "Weight(kg)": "weight_kg",
    "Systolic Blood Pressure(mmHg)": "sbp",
    "Diastolic Blood Pressure(mmHg)": "dbp",
    "Heart Rate(b/m)": "heart_rate_bpm",
    "BMI(kg/m^2)": "bmi_kg_m2",
    "Hypertension": "hypertension_stage",
    "Diabetes": "diabetes",
    "cerebral infarction": "cerebral_infarction",
    "cerebrovascular disease": "cerebrovascular_disease",
}


def _read_info(raw_dir: Path) -> pd.DataFrame:
    """Read the subject table (one row per subject)."""
    info = pd.read_excel(raw_dir / EXCEL_FILE, header=1)  # first row is a title
    info.columns = info.columns.str.strip()
    info = info.rename(columns=COLUMNS)
    info["subject_id"] = DATASET + "_" + info["subject"].astype(int).astype(str)
    return info


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    info = _read_info(raw_dir)
    return pd.DataFrame({
        "subject_id": info["subject_id"],
        "dataset": DATASET,
        "age": info["age"],
        "sex": normalize_sex(info["sex"]),
        "height_cm": info["height_cm"],
        "weight_kg": info["weight_kg"],
        "heart_rate_bpm": info["heart_rate_bpm"],
        "bmi_kg_m2": info["bmi_kg_m2"],
        "hypertension_stage": info["hypertension_stage"],
        "diabetes": info["diabetes"],
        "cerebral_infarction": info["cerebral_infarction"],
        "cerebrovascular_disease": info["cerebrovascular_disease"],
    })


def write_records(raw_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Convert every 2.1 s segment to one HDF5 file and return the records table."""
    rows = []
    for txt in sorted((raw_dir / WAVE_DIR).glob("*.txt")):
        subject, segment = txt.stem.split("_")                 # "2_1.txt"
        data = np.array(txt.read_text().split(), dtype=float)  # one line, tab-separated, trailing tab

        record_id = f"{DATASET}_{subject}_{segment}"
        wave_path = f"waves/{record_id}.h5"
        write_wave(out_dir / wave_path, {"ppg": (data, SAMPLE_RATE_HZ, "a.u.", 0.0)})

        rows.append({
            "record_id": record_id,
            "subject_id": f"{DATASET}_{subject}",
            "dataset": DATASET,
            "condition": "rest",
            "duration_s": len(data) / SAMPLE_RATE_HZ,
            "signals": ["ppg"],
            "sampling_rates_hz": [SAMPLE_RATE_HZ],
            "wave_path": wave_path,
        })
    return pd.DataFrame(rows)


def build_bp(raw_dir: Path, records: pd.DataFrame) -> pd.DataFrame:
    """One cuff reading per subject, copied to each of its 3 records."""
    info = _read_info(raw_dir)

    bp = records[["record_id", "subject_id"]].merge(info[["subject_id", "sbp", "dbp"]], on="subject_id")
    bp["time_s"] = 0.0  # the paper: BP was measured at the start
    bp["map"] = mean_arterial_pressure(bp["sbp"], bp["dbp"])
    bp["method"] = "cuff"
    return bp[["record_id", "time_s", "sbp", "dbp", "map", "method"]]


def convert(raw_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    build_subjects(raw_dir).to_parquet(out_dir / "subjects.parquet", index=False)

    records = write_records(raw_dir, out_dir)
    records.to_parquet(out_dir / "records.parquet", index=False)

    build_bp(raw_dir, records).to_parquet(out_dir / "bp.parquet", index=False)