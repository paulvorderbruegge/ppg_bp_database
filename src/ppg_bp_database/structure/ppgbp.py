"""PPG-BP dataset (Liang et al. 2018, figshare): conversion to the shared format.

Raw:     data/raw/ppgbp/PPG-BP Database/Data File/{PPG-BP dataset.xlsx, 0_subject/<subject>_<segment>.txt}
Records: one per 2.1 s PPG segment (3 per subject), 1000 Hz
BP:      one cuff reading per subject, copied to each of its records
"""
from pathlib import Path

import numpy as np
import pandas as pd

from ppg_bp_database.utils import mean_arterial_pressure, normalize_sex, prepare_out_dir, save_tables, write_wave

DATASET = "ppgbp"
SYNTHETIC = False  # measured data
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
EXTRA_SUBJECT_COLUMNS = ["heart_rate_bpm", "bmi_kg_m2", "hypertension_stage", "diabetes",
                         "cerebral_infarction", "cerebrovascular_disease"]


def _read_info(raw_dir: Path) -> pd.DataFrame:
    """Read the subject table (one row per subject)."""
    info = pd.read_excel(raw_dir / EXCEL_FILE, header=1)  # first row is a title
    info.columns = info.columns.str.strip()
    missing = set(COLUMNS) - set(info.columns)
    assert not missing, f"columns not found in the Excel file: {missing}"
    info = info.rename(columns=COLUMNS)
    info["subject_id"] = DATASET + "_" + info["subject"].astype(int).astype(str)
    return info


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    info = _read_info(raw_dir)
    subjects = info[["subject_id", "age", "height_cm", "weight_kg", *EXTRA_SUBJECT_COLUMNS]].copy()
    subjects["dataset"] = DATASET
    subjects["sex"] = normalize_sex(info["sex"])
    return subjects


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
    assert rows, f"no PPG-BP .txt files found in {raw_dir / WAVE_DIR}"
    return pd.DataFrame(rows)


def build_bp(raw_dir: Path, records: pd.DataFrame) -> pd.DataFrame:
    """One cuff reading per subject, copied to each of its 3 records."""
    info = _read_info(raw_dir)

    bp = records[["record_id", "subject_id"]].merge(info[["subject_id", "sbp", "dbp"]], on="subject_id")
    bp["time_s"] = 0.0  # the paper: BP was measured at the start
    bp["map"] = mean_arterial_pressure(bp["sbp"], bp["dbp"])
    bp["method"] = "cuff"
    bp["source"] = "dataset"
    return bp.drop(columns="subject_id")


def convert(raw_dir: Path, out_dir: Path) -> None:
    prepare_out_dir(out_dir)
    records = write_records(raw_dir, out_dir)
    save_tables(out_dir, build_subjects(raw_dir), records, build_bp(raw_dir, records), synthetic=SYNTHETIC)
