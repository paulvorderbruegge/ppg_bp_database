"""Graphene bioimpedance BP dataset (Kireev et al., PhysioNet v1.0.0): conversion to the shared format.

Raw:     data/raw/bpgraphene/**/subject<N>_day<D>/setup<NN>_<condition>/data_trial<TT>_<signal>.csv
         (found wherever the PhysioNet zip was unpacked below data/raw/bpgraphene)
Records: one per trial; signals ppg, finapres_ppg, finapres_bp (bioimpedance is deliberately skipped)
BP:      one value per beat, calculated from the finapres_bp wave (ppg_bp_database.beats), source "calculated"
"""
from pathlib import Path

import pandas as pd

from ppg_bp_database.beats import bp_from_wave
from ppg_bp_database.utils import prepare_out_dir, save_tables, write_wave

DATASET = "bpgraphene"
SYNTHETIC = False  # measured data

# file suffix -> (signal name, documented sampling rate in Hz, unit)
SIGNALS = {
    "ppg": ("ppg", 1250.0, "a.u."),
    "finapresPPG": ("finapres_ppg", 75.0, "a.u."),
    "finapresBP": ("finapres_bp", 200.0, "mmHg"),
}
FS_TOLERANCE = 0.01  # measured sampling rate may differ by 1 % from the documented one


def _trial_files(raw_dir: Path) -> list[Path]:
    """The finapresBP file of every trial; the other signals are next to it."""
    files = sorted(raw_dir.rglob("subject*_day*/setup*/data_trial*_finapresBP.csv"))
    assert files, f"no bp-graphene trial files found below {raw_dir} (zip not unpacked?)"
    return files


def _subject_day(bp_file: Path) -> tuple[str, str]:
    """.../subject1_day2/setup02_hgcp/data_trial01_finapresBP.csv -> ('1', '2')"""
    subject, day = bp_file.parent.parent.name.removeprefix("subject").split("_day")
    return subject, day


def _read_signal(path: Path, fs: float) -> tuple[pd.Series, float]:
    """Signal values and start time. First column = time in s, second column = the signal."""
    df = pd.read_csv(path)
    time = df.iloc[:, 0]
    fs_measured = 1 / time.diff().median()
    assert abs(fs_measured / fs - 1) < FS_TOLERANCE, \
        f"{path}: measured {fs_measured:.2f} Hz, documented {fs} Hz"
    return df.iloc[:, 1], float(time.iloc[0])


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    """One row per subject. The dataset gives no demographics."""
    numbers = sorted({int(_subject_day(f)[0]) for f in _trial_files(raw_dir)})
    return pd.DataFrame({"subject_id": [f"{DATASET}_s{n}" for n in numbers], "dataset": DATASET})


def _record_id(bp_file: Path) -> str:
    subject, day = _subject_day(bp_file)
    setup = bp_file.parent.name.split("_")[0]                     # "setup02"
    trial = bp_file.stem.split("_")[1].removeprefix("trial")      # "01"
    return f"{DATASET}_s{subject}_day{day}_{setup}_trial{trial}"


def write_records(raw_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Convert every trial to one HDF5 file and return the records table."""
    rows = []
    for bp_file in _trial_files(raw_dir):
        subject, day = _subject_day(bp_file)
        setup = bp_file.parent.name                                   # "setup02_hgcp"
        trial = bp_file.stem.split("_")[1].removeprefix("trial")      # "01"

        record_id = _record_id(bp_file)
        wave_path = f"waves/{record_id}.h5"

        signals = {}
        for suffix, (name, fs, unit) in SIGNALS.items():
            values, t0 = _read_signal(bp_file.with_name(f"data_trial{trial}_{suffix}.csv"), fs)
            signals[name] = (values.to_numpy(), fs, unit, t0)
        write_wave(out_dir / wave_path, signals)

        rows.append({
            "record_id": record_id,
            "subject_id": f"{DATASET}_s{subject}",
            "dataset": DATASET,
            "session": f"day{day}",
            "condition": setup.split("_", 1)[1],                         # "baseline", "hgcp", "valsalva", ...
            "duration_s": max(len(d) / fs for d, fs, _, _ in signals.values()),
            "signals": list(signals),
            "sampling_rates_hz": [fs for _, fs, _, _ in signals.values()],
            "wave_path": wave_path,
        })
    return pd.DataFrame(rows)


def build_bp(raw_dir: Path, records: pd.DataFrame) -> pd.DataFrame:
    """One value per beat of the finapres_bp wave (the dataset gives no single BP values)."""
    _, fs, _ = SIGNALS["finapresBP"]
    parts = []
    for bp_file in _trial_files(raw_dir):
        values, t0 = _read_signal(bp_file, fs)
        parts.append(bp_from_wave(_record_id(bp_file), values.to_numpy(), fs, t0, "finapres"))
    return pd.concat(parts, ignore_index=True)


def convert(raw_dir: Path, out_dir: Path) -> None:
    prepare_out_dir(out_dir)
    records = write_records(raw_dir, out_dir)
    save_tables(out_dir, build_subjects(raw_dir), records, build_bp(raw_dir, records), synthetic=SYNTHETIC)
