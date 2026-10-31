"""MIMIC-BP dataset (Sanches et al. 2024, Harvard Dataverse): conversion to the shared format.

Raw:     data/raw/mimicbp/**/<source_id>_{ppg,abp,labels}.npy (found wherever the zips were unpacked)
Records: one per 30 s segment (30 per subject), ppg + abp at 125 Hz
BP:      median SBP/DBP per segment as provided by the dataset (derived from the ABP)
         TODO: later many single values per record, derived from the abp wave; the labels can validate them
"""
from pathlib import Path

import numpy as np
import pandas as pd

from ppg_bp_database.utils import mean_arterial_pressure, prepare_out_dir, save_tables, write_wave

DATASET = "mimicbp"
SYNTHETIC = False  # measured data
SAMPLE_RATE_HZ = 125.0
N_SEGMENTS = 30         # per subject
SEGMENT_SAMPLES = 3750  # 30 s at 125 Hz


def _find_files(raw_dir: Path) -> dict[str, Path]:
    """All .npy files by name, wherever the zips were unpacked."""
    files = {p.name: p for p in raw_dir.rglob("*.npy")}
    assert files, f"no MIMIC-BP .npy files found below {raw_dir}"
    return files


def _source_ids(files: dict[str, Path]) -> list[str]:
    """Source subject IDs like 'p093833'."""
    return sorted(name.removesuffix("_ppg.npy") for name in files if name.endswith("_ppg.npy"))


def _record_id(source_id: str, segment: int) -> str:
    return f"{DATASET}_{source_id}_{segment:02d}"  # segment 1..30


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    """One row per subject (1524).
    TODO: age, sex, ... from the MIMIC-III clinical database (source id "p093833" = SUBJECT_ID 93833)."""
    ids = _source_ids(_find_files(raw_dir))
    return pd.DataFrame({"subject_id": [f"{DATASET}_{i}" for i in ids], "dataset": DATASET})


def write_records(raw_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Convert every 30 s segment to one HDF5 file (ppg + abp) and return the records table."""
    files = _find_files(raw_dir)

    rows = []
    for source_id in _source_ids(files):
        ppg = np.load(files[f"{source_id}_ppg.npy"])
        abp = np.load(files[f"{source_id}_abp.npy"])
        assert ppg.shape == abp.shape == (N_SEGMENTS, SEGMENT_SAMPLES), \
            f"{source_id}: ppg {ppg.shape}, abp {abp.shape}"  # one row per segment

        for i in range(N_SEGMENTS):
            record_id = _record_id(source_id, i + 1)
            wave_path = f"waves/{record_id}.h5"
            write_wave(out_dir / wave_path, {
                "ppg": (ppg[i], SAMPLE_RATE_HZ, "a.u.", 0.0),
                "abp": (abp[i], SAMPLE_RATE_HZ, "mmHg", 0.0),
            })

            rows.append({
                "record_id": record_id,
                "subject_id": f"{DATASET}_{source_id}",
                "dataset": DATASET,
                "condition": "icu",
                "duration_s": SEGMENT_SAMPLES / SAMPLE_RATE_HZ,
                "signals": ["ppg", "abp"],
                "sampling_rates_hz": [SAMPLE_RATE_HZ, SAMPLE_RATE_HZ],
                "wave_path": wave_path,
            })
    return pd.DataFrame(rows)


def build_bp(raw_dir: Path, records: pd.DataFrame) -> pd.DataFrame:
    """Median SBP/DBP per segment, as provided by the dataset."""
    files = _find_files(raw_dir)

    parts = []
    for source_id in _source_ids(files):
        labels = np.load(files[f"{source_id}_labels.npy"])
        assert labels.shape == (N_SEGMENTS, 2), f"{source_id}: labels {labels.shape}, expected (30, 2): SBP, DBP"
        parts.append(pd.DataFrame({
            "record_id": [_record_id(source_id, i + 1) for i in range(N_SEGMENTS)],
            "sbp": labels[:, 0],
            "dbp": labels[:, 1],
        }))

    bp = pd.concat(parts, ignore_index=True)
    bp["time_s"] = np.nan  # median over the whole 30 s segment, not a single point in time
    bp["map"] = mean_arterial_pressure(bp["sbp"], bp["dbp"])
    bp["method"] = "derived_from_abp"
    return bp


def convert(raw_dir: Path, out_dir: Path) -> None:
    prepare_out_dir(out_dir)
    records = write_records(raw_dir, out_dir)
    save_tables(out_dir, build_subjects(raw_dir), records, build_bp(raw_dir, records), synthetic=SYNTHETIC)
