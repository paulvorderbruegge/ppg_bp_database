"""PulseDB dataset (Wang et al. 2023, v2_0, MIMIC-III + VitalDB): Conversion to the shared format.
Notable: One record = one consecutive clip of one .mat file (NOT one 10 s segment: that would be ~5 million HDF5 files).
         The 10 s segments that PulseDB removed for bad quality leave gaps inside a clip; they are filled with NaN, so
         the time axis (t0 + i / fs) stays valid and equals PulseDB's own clip time `T`.
"""
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
from mat73 import loadmat

from ppg_bp_database.utils import normalize_sex, write_wave

DATASET = "pulsedb"
SAMPLE_RATE_HZ = 125.0
SEGMENT_SAMPLES = 1250  # 10 s at 125 Hz
SEGMENT_S = SEGMENT_SAMPLES / SAMPLE_RATE_HZ

# folder name -> (source tag used in the IDs, condition)
SOURCES = {
    "PulseDB_MIMIC": ("mimic", "icu"),
    "PulseDB_Vital": ("vital", "surgery"),
}

# our signal name -> (field in the .mat file, unit)
SIGNALS = {
    "ecg": ("ECG_Record", "a.u."),  # TODO: check unit (mV?)
    "ppg": ("PPG_Record", "a.u."),
    "abp": ("ABP_Raw", "mmHg"),
}

# other fields we read
FIELDS = {
    "time": "T",           # (n_segments, 1250), absolute time in s within the clip
    "age": "Subj_Age",
    "sex": "Subj_Gender",  # 'M' / 'F'
    "height": "Subj_Height",
    "weight": "Subj_Weight",
}


def _load(path: Path) -> dict:
    """Load a v7.3 .mat file; if everything is wrapped in one struct, return that struct."""
    data = loadmat(str(path))
    if len(data) == 1:
        inner = next(iter(data.values()))
        if isinstance(inner, dict):
            return inner
    return data


def _matrix(x) -> np.ndarray:
    """(n_segments, 1250) float array, whether the file gave a matrix or a list of rows."""
    return np.asarray(x, dtype=float).reshape(-1, SEGMENT_SAMPLES)


def _first(x):
    """Subject-level value, whether stored once or repeated for every segment."""
    return np.asarray(x).reshape(-1)[0]


def _split_clips(starts: np.ndarray) -> list[np.ndarray]:
    """Segment indices per clip. T restarts at 0 for every clip, so a new clip begins where the start time
    stops increasing."""
    breaks = np.flatnonzero(np.diff(starts) <= 0) + 1
    return np.split(np.arange(len(starts)), breaks)


def convert_file(mat_path: Path, out_dir: Path, source: str, condition: str
                 ) -> tuple[dict, list[dict]]:
    """Write the HDF5 file(s) of one .mat and return (subject row, record rows)."""
    d = _load(mat_path)

    time = _matrix(d[FIELDS["time"]])
    n = len(time)
    starts = time[:, 0]
    assert np.allclose(np.diff(time, axis=1), 1 / SAMPLE_RATE_HZ, atol=1e-3), "T is not 125 Hz in seconds"

    data = {name: _matrix(d[field]) for name, (field, _) in SIGNALS.items()}
    for name, m in data.items():
        assert m.shape == (n, SEGMENT_SAMPLES), f"{mat_path.name}: {name} has shape {m.shape}"

    subject_id = f"{DATASET}_{source}_{mat_path.stem}"
    subject = {
        "subject_id": subject_id,
        "dataset": DATASET,
        "age": float(_first(d[FIELDS["age"]])),
        "sex": normalize_sex(pd.Series([str(_first(d[FIELDS["sex"]]))])).iloc[0],
        "height_cm": float(_first(d[FIELDS["height"]])),  # NaN for MIMIC subjects
        "weight_kg": float(_first(d[FIELDS["weight"]])),
    }

    records = []
    for clip, idx in enumerate(_split_clips(starts), start=1):
        t0 = float(starts[idx[0]])
        offsets = np.rint((starts[idx] - t0) * SAMPLE_RATE_HZ).astype(int)
        assert np.all(np.diff(offsets) >= SEGMENT_SAMPLES), "overlapping segments"
        total = int(offsets[-1]) + SEGMENT_SAMPLES

        signals = {}
        for name, (_, unit) in SIGNALS.items():
            full = np.full(total, np.nan)  # NaN where PulseDB removed a bad segment
            for off, row in zip(offsets, data[name][idx]):
                full[off:off + SEGMENT_SAMPLES] = row
            signals[name] = (full, SAMPLE_RATE_HZ, unit, t0)

        record_id = f"{subject_id}_c{clip:02d}"
        wave_path = f"waves/{record_id}.h5"
        write_wave(out_dir / wave_path, signals)

        records.append({
            "record_id": record_id,
            "subject_id": subject_id,
            "dataset": DATASET,
            "condition": condition,
            "duration_s": total / SAMPLE_RATE_HZ,
            "signals": list(signals),
            "sampling_rates_hz": [SAMPLE_RATE_HZ] * len(signals),
            "wave_path": wave_path,
        })

    return subject, records


def _mat_files(raw_dir: Path) -> Iterator[tuple[Path, str, str]]:
    """(file, source tag, condition) for all subject files, wherever the archives were unpacked."""
    for folder, (source, condition) in SOURCES.items():
        for directory in sorted(raw_dir.rglob(folder)):
            if directory.is_dir():
                for path in sorted(directory.glob("*.mat")):
                    yield path, source, condition


def convert(raw_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    files = list(_mat_files(raw_dir))
    assert files, f"no PulseDB .mat files found below {raw_dir}"

    subjects, records = [], []
    for i, (path, source, condition) in enumerate(files, start=1):
        subject, recs = convert_file(path, out_dir, source, condition)
        subjects.append(subject)
        records.extend(recs)
        if i % 100 == 0:
            print(f"{i}/{len(files)} files converted")

    pd.DataFrame(subjects).to_parquet(out_dir / "subjects.parquet", index=False)
    pd.DataFrame(records).to_parquet(out_dir / "records.parquet", index=False)

    # TODO: derive BP values (beat-wise SBP/DBP/MAP) from the abp wave and write bp.parquet
