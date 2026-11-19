"""MIMIC-BP dataset (Sanches et al. 2024, Harvard Dataverse): conversion to the shared format.

Raw:     data/raw/mimicbp/**/<source_id>_{ppg,abp,labels}.npy (found wherever the zips were unpacked)
Records: one per 30 s segment (30 per subject), ppg + abp at 125 Hz
BP:      source "dataset":    median SBP/DBP per segment as provided by the dataset (derived from the abp), time_s NaN
         source "calculated": one value per beat, calculated from the abp wave (structure/beats.py);
                              the dataset's medians can be used to validate them
"""
from pathlib import Path

import numpy as np
import pandas as pd

from ppg_bp_database.structure import mimic3_clinical
from ppg_bp_database.structure.beats import bp_from_wave
from ppg_bp_database.utils import mean_arterial_pressure, prepare_out_dir, progress, save_tables, write_wave

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


def _subject_id(source_id: str) -> str:
    """'p093833' -> 'mimic3_p093833': shared with the MIMIC part of PulseDB (same person, same subject_id)."""
    return mimic3_clinical.subject_id(source_id)


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    """One row per subject (1524), subject_id "mimic3_<source id>" (see mimic3_clinical, "Subject IDs"). Age, sex, height and weight from the MIMIC-III clinical database
    (source id "p093833" = SUBJECT_ID 93833; see structure/mimic3_clinical.py). MIMIC-BP does not say which
    waveform record a subject's segments come from, so all of the subject's records in the Matched Subset are used
    to find the admission (extra columns mimic3_hadm_id, mimic3_match, age_capped).
    Without data/raw/mimic3 the demographics stay empty."""
    ids = pd.Series(_source_ids(_find_files(raw_dir)), dtype="string")
    subjects = pd.DataFrame({"subject_id": mimic3_clinical.subject_id(ids), "dataset": DATASET})
    if not mimic3_clinical.available():
        print(f"  {mimic3_clinical.missing_message()}")
        return subjects

    query = pd.DataFrame({"key": subjects["subject_id"], "SUBJECT_ID": mimic3_clinical.subject_number(ids)})
    query = query.merge(mimic3_clinical.record_times(), on="SUBJECT_ID", how="left")  # time NaT if no record
    demo = mimic3_clinical.demographics(query)
    print(f"  MIMIC-III admission match: {demo['mimic3_match'].value_counts().to_dict()}")
    return subjects.join(demo, on="subject_id")


def write_records(raw_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Convert every 30 s segment to one HDF5 file (ppg + abp) and return the records table."""
    files = _find_files(raw_dir)

    rows = []
    for source_id in progress(_source_ids(files), "files converted"):
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
                "subject_id": _subject_id(source_id),
                "dataset": DATASET,
                "condition": "icu",
                "duration_s": SEGMENT_SAMPLES / SAMPLE_RATE_HZ,
                "signals": ["ppg", "abp"],
                "sampling_rates_hz": [SAMPLE_RATE_HZ, SAMPLE_RATE_HZ],
                "wave_path": wave_path,
            })
    return pd.DataFrame(rows)


def build_bp(raw_dir: Path, records: pd.DataFrame) -> pd.DataFrame:
    """Median SBP/DBP per segment as provided by the dataset, plus one calculated value per beat of the abp wave."""
    files = _find_files(raw_dir)

    given, beats = [], []
    for source_id in progress(_source_ids(files), "files: BP extracted"):
        record_ids = [_record_id(source_id, i + 1) for i in range(N_SEGMENTS)]
        labels = np.load(files[f"{source_id}_labels.npy"])
        assert labels.shape == (N_SEGMENTS, 2), f"{source_id}: labels {labels.shape}, expected (30, 2): SBP, DBP"
        given.append(pd.DataFrame({"record_id": record_ids, "sbp": labels[:, 0], "dbp": labels[:, 1]}))

        abp = np.load(files[f"{source_id}_abp.npy"])
        beats += [bp_from_wave(rid, abp[i], SAMPLE_RATE_HZ, 0.0, "abp") for i, rid in enumerate(record_ids)]

    given = pd.concat(given, ignore_index=True)
    given["time_s"] = np.nan  # median over the whole 30 s segment, not a single point in time
    given["map"] = mean_arterial_pressure(given["sbp"], given["dbp"])
    given["method"] = "abp"
    given["source"] = "dataset"
    return pd.concat([given, *beats], ignore_index=True)


def convert(raw_dir: Path, out_dir: Path) -> None:
    prepare_out_dir(out_dir)
    records = write_records(raw_dir, out_dir)
    save_tables(out_dir, build_subjects(raw_dir), records, build_bp(raw_dir, records), synthetic=SYNTHETIC)
