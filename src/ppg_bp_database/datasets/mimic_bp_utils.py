"""MIMIC-BP dataset (Sanches et al. 2024, Harvard Dataverse): Conversion to the shared format."""
from pathlib import Path

import numpy as np
import pandas as pd

from ppg_bp_database.utils import mean_arterial_pressure, write_wave

DATASET = "mimicbp"
SAMPLE_RATE_HZ = 125.0
N_SEGMENTS = 30         # per subject
SEGMENT_SAMPLES = 3750  # 30 s at 125 Hz


def _find_files(raw_dir: Path) -> dict[str, Path]:
    """All .npy files by name, wherever the zips were unpacked."""
    return {p.name: p for p in raw_dir.rglob("*.npy")}


def _source_ids(files: dict[str, Path]) -> list[str]:
    """Source subject IDs like 'p093833'."""
    return sorted(name.removesuffix("_ppg.npy") for name in files if name.endswith("_ppg.npy"))


def _record_id(source_id: str, segment: int) -> str:
    return f"{DATASET}_{source_id}_{segment:02d}"  # segment 1..30


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    """One row per subject (1524). TODO: Get data from original dataset."""
    ids = _source_ids(_find_files(raw_dir))
    return pd.DataFrame({
        "subject_id": [f"{DATASET}_{i}" for i in ids],
        "dataset": DATASET,
        "age": np.nan,
        "sex": None,
        "height_cm": np.nan,
        "weight_kg": np.nan,
    })


def write_records(raw_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Convert every 30 s segment to one HDF5 file (ppg + abp) and return the records table."""
    files = _find_files(raw_dir)

    rows = []
    for source_id in _source_ids(files):
        ppg = np.load(files[f"{source_id}_ppg.npy"])
        abp = np.load(files[f"{source_id}_abp.npy"])
        assert ppg.shape == abp.shape == (N_SEGMENTS, SEGMENT_SAMPLES)  # one row per segment

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


def build_bp(raw_dir: Path) -> pd.DataFrame:
    """Median SBP/DBP per segment, as provided by the dataset (derived from the ABP)."""
    files = _find_files(raw_dir)

    parts = []
    for source_id in _source_ids(files):
        labels = np.load(files[f"{source_id}_labels.npy"])  # expected (30, 2): SBP, DBP
        parts.append(pd.DataFrame({
            "record_id": [_record_id(source_id, i + 1) for i in range(N_SEGMENTS)],
            "sbp": labels[:, 0],
            "dbp": labels[:, 1],
        }))

    bp = pd.concat(parts, ignore_index=True)
    bp["time_s"] = np.nan  # median over the whole 30 s segment, not a single point in time
    bp["map"] = mean_arterial_pressure(bp["sbp"], bp["dbp"])
    bp["method"] = "derived_from_abp"
    return bp[["record_id", "time_s", "sbp", "dbp", "map", "method"]]


def convert(raw_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    build_subjects(raw_dir).to_parquet(out_dir / "subjects.parquet", index=False)
    write_records(raw_dir, out_dir).to_parquet(out_dir / "records.parquet", index=False)
    build_bp(raw_dir).to_parquet(out_dir / "bp.parquet", index=False)