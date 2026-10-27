"""MIMIC-BP dataset (Sanches et al. 2024, Harvard Dataverse): Conversion to the shared format."""
from pathlib import Path

import numpy as np
import pandas as pd

from ppg_bp_database.utils import write_wave

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
    """One row per subject (1524)."""
    ids = _source_ids(_find_files(raw_dir))
    # TODO: age, sex, ... from the MIMIC-III clinical database (source id "p093833" = SUBJECT_ID 93833)
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


def convert(raw_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    build_subjects(raw_dir).to_parquet(out_dir / "subjects.parquet", index=False)
    write_records(raw_dir, out_dir).to_parquet(out_dir / "records.parquet", index=False)

    # TODO: derive BP values (many per record) from the abp waves and write bp.parquet;
    #       the label files (median SBP/DBP per segment) can serve as a reference to validate them