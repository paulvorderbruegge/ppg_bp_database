"""Pulse Wave Database (Charlton et al. 2019, Zenodo record 3275625, v0.2.0): conversion to the shared format.

SYNTHETIC: the waves are not measured but simulated with a 1D model of the arterial tree, for 4374 virtual subjects
aged 25-75 (6 baseline subjects, one per age decade, plus variations of HR, SV, LVET, DBP, MBP and PWV by -1/0/+1 SD).

Raw:     data/raw/pwdb/{**/pwdb<NNNN>.hea/.dat (WFDB, from PWs_wfdb.zip), pwdb_haemod_params.csv,
         pwdb_model_configs.csv, pwdb_model_variations.csv}   (download: scripts/download_pwdb.sh)
Records: one per virtual subject = one steady-state cardiac cycle, 500 Hz, 13 sites x 4 waves:
         abp_<site> (pressure, mmHg), flow_velocity_<site> (m/s), area_<site> (luminal area, m^2), ppg_<site> (a.u.)
         The last sample of some channels is stored as "invalid" in the WFDB files and therefore NaN.
BP:      brachial SBP/DBP per subject from pwdb_haemod_params.csv (taken from the simulated brachial pressure wave);
         time_s = NaN, the beat is periodic so the value describes the whole record
"""
from pathlib import Path

import pandas as pd
import wfdb

from ppg_bp_database.utils import mean_arterial_pressure, normalize_sex, prepare_out_dir, save_tables, write_wave

DATASET = "pwdb"
SYNTHETIC = True
SAMPLE_RATE_HZ = 500.0
N_SUBJECTS = 4374
HAEMOD_FILE = "pwdb_haemod_params.csv"
CONFIG_FILE = "pwdb_model_configs.csv"
VARIATION_FILE = "pwdb_model_variations.csv"

# WFDB site name -> our site name
SITES = {
    "AorticRoot": "aortic_root",
    "Carotid": "carotid",
    "ThorAorta": "thoracic_aorta",
    "Brachial": "brachial",
    "Radial": "radial",
    "AbdAorta": "abdominal_aorta",
    "IliacBif": "iliac_bifurcation",
    "CommonIliac": "common_iliac",
    "Femoral": "femoral",
    "AntTibial": "anterior_tibial",
    "SupMidCerebral": "sup_middle_cerebral",
    "SupTemporal": "sup_temporal",
    "Digital": "digital",
}
# WFDB wave suffix -> (our signal prefix, expected WFDB unit, our unit)
WAVES = {
    "P": ("abp", "mmHg", "mmHg"),
    "U": ("flow_velocity", "m_per_s", "m/s"),
    "A": ("area", "m2", "m^2"),
    "PPG": ("ppg", "au", "a.u."),
}

# pwdb_haemod_params.csv column -> our subject column
HAEMOD_COLUMNS = {
    "HR [bpm]": "heart_rate_bpm",
    "SV [ml]": "stroke_volume_ml",
    "CO [l/min]": "cardiac_output_l_min",
    "LVET [ms]": "lvet_ms",
    "SBP_a [mmHg]": "sbp_aortic_mmhg",
    "DBP_a [mmHg]": "dbp_aortic_mmhg",
    "PWV_cf [m/s]": "pwv_cf_m_s",
    "PWV_br [m/s]": "pwv_br_m_s",
    "AIx [%]": "aix_percent",
    "SVR [10^6 Pa s / m3]": "svr_1e6_pa_s_m3",
}
# pwdb_model_variations.csv column -> our subject column; values are the prescribed deviation in SD (-1, 0, +1)
VARIATION_COLUMNS = {p: f"variation_{p.lower()}_sd" for p in ["DIA", "HR", "LVET", "MBP", "PWV", "SV"]}


def _read_csv(raw_dir: Path, name: str) -> pd.DataFrame:
    """One row per subject; the headers have spaces after the commas."""
    df = pd.read_csv(raw_dir / name, skipinitialspace=True)
    df = df.rename(columns={df.columns[0]: "subject"})  # "Subject Number" / "SUBJECT NUMBER"
    assert len(df) == N_SUBJECTS and df["subject"].is_unique, f"{name}: {len(df)} rows, expected {N_SUBJECTS}"
    return df.set_index("subject")


def _subject_id(number: int) -> str:
    return f"{DATASET}_{number:04d}"


def _headers(raw_dir: Path) -> dict[int, Path]:
    """Subject number -> path of its WFDB record (without extension), wherever PWs_wfdb.zip was unpacked."""
    files = {int(p.stem.removeprefix("pwdb")): p.with_suffix("") for p in raw_dir.rglob("pwdb[0-9]*.hea")}
    assert len(files) == N_SUBJECTS, \
        f"found {len(files)} PWDB .hea files below {raw_dir}, expected {N_SUBJECTS} (PWs_wfdb.zip not unpacked?)"
    return dict(sorted(files.items()))


def _sex(header: wfdb.Record) -> str:
    """The header comment looks like '<age>: 55 <sex>: male'."""
    comment = " ".join(header.comments)  # type: ignore
    assert "<sex>:" in comment, f"{header.record_name}: no sex in comments {header.comments}"
    return comment.split("<sex>:")[1].split()[0]


def build_subjects(raw_dir: Path) -> pd.DataFrame:
    """One row per virtual subject: age from the model config, sex from the WFDB header (the model is male),
    haemodynamic parameters and the prescribed variations as extra columns."""
    haemod = _read_csv(raw_dir, HAEMOD_FILE)
    config = _read_csv(raw_dir, CONFIG_FILE)
    variations = _read_csv(raw_dir, VARIATION_FILE)
    sex = pd.Series({n: _sex(wfdb.rdheader(str(p))) for n, p in _headers(raw_dir).items()})

    subjects = pd.DataFrame({
        "subject_id": [_subject_id(n) for n in haemod.index],
        "dataset": DATASET,
        "age": config.loc[haemod.index, "age [years]"].to_numpy(),
        "sex": normalize_sex(sex.reindex(haemod.index)).to_numpy(),
        "baseline": config.loc[haemod.index, "base_age"].eq(1).to_numpy(),  # baseline subject of its age group
    })
    for src, dst in HAEMOD_COLUMNS.items():
        subjects[dst] = haemod[src].to_numpy()
    for src, dst in VARIATION_COLUMNS.items():
        subjects[dst] = variations.loc[haemod.index, src].to_numpy()
    return subjects


def write_records(raw_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Convert every subject's WFDB record (one beat, 52 channels) to one HDF5 file and return the records table."""
    rows = []
    for number, path in _headers(raw_dir).items():
        rec = wfdb.rdrecord(str(path))
        assert rec.fs == SAMPLE_RATE_HZ, f"{path.name}: {rec.fs} Hz"

        signals = {}
        for i, (raw_name, unit) in enumerate(zip(rec.sig_name, rec.units)):  # type: ignore
            site, wave = raw_name.rstrip(", ").rsplit("_", 1)                # "AorticRoot_P," -> AorticRoot, P
            prefix, expected_unit, our_unit = WAVES[wave]
            assert unit == expected_unit, f"{path.name} {raw_name}: unit {unit}, expected {expected_unit}"
            signals[f"{prefix}_{SITES[site]}"] = (rec.p_signal[:, i], SAMPLE_RATE_HZ, our_unit, 0.0)  # type: ignore
        assert len(signals) == len(SITES) * len(WAVES), f"{path.name}: {len(signals)} signals"

        record_id = _subject_id(number)  # one record per subject
        wave_path = f"waves/{record_id}.h5"
        write_wave(out_dir / wave_path, signals)

        rows.append({
            "record_id": record_id,
            "subject_id": _subject_id(number),
            "dataset": DATASET,
            "condition": "rest",
            "duration_s": rec.sig_len / SAMPLE_RATE_HZ,  # type: ignore
            "signals": list(signals),
            "sampling_rates_hz": [SAMPLE_RATE_HZ] * len(signals),
            "wave_path": wave_path,
        })
        if len(rows) % 500 == 0:
            print(f"{len(rows)}/{N_SUBJECTS} records converted")
    return pd.DataFrame(rows)


def build_bp(raw_dir: Path, records: pd.DataFrame) -> pd.DataFrame:
    """Brachial SBP/DBP of the simulated pressure wave, one value per record (= per subject)."""
    haemod = _read_csv(raw_dir, HAEMOD_FILE)
    bp = pd.DataFrame({
        "record_id": [_subject_id(n) for n in haemod.index],
        "time_s": float("nan"),  # periodic beat: the value describes the whole record
        "sbp": haemod["SBP_b [mmHg]"].to_numpy(),
        "dbp": haemod["DBP_b [mmHg]"].to_numpy(),
    })
    bp["map"] = mean_arterial_pressure(bp["sbp"], bp["dbp"])
    bp["method"] = "derived_from_abp"
    return bp[bp["record_id"].isin(records["record_id"])]


def convert(raw_dir: Path, out_dir: Path) -> None:
    prepare_out_dir(out_dir)
    records = write_records(raw_dir, out_dir)
    save_tables(out_dir, build_subjects(raw_dir), records, build_bp(raw_dir, records), synthetic=SYNTHETIC)
