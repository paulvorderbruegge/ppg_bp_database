# ppg-bp-database

Converts several PPG / blood-pressure datasets into one shared format, so that notebooks can use them all the same way.

## Layout

```
data/raw/<dataset>/            downloaded data, untouched
data/structured/<dataset>/     converted data (written by scripts/convert.py)
    subjects.parquet           one row per subject
    records.parquet            one row per recording
    bp.parquet                 one row per BP value (missing if the dataset has none yet)
    waves/<record_id>.h5       one HDF5 file per record, one HDF5 dataset per signal
src/ppg_bp_database/
    utils.py                   shared schema, paths, save_tables, write_wave / read_wave
    access.py                  loading functions for notebooks
    datasets/<dataset>.py      one converter per dataset, all with the same functions
scripts/
    convert.py                 raw -> structured, for one, several or all datasets
    download_all.sh            runs the download_<dataset>.sh scripts below (all, or the ones given)
    download_<dataset>.sh      download + checksum check + unpack of one dataset into data/raw/<dataset>
    physionet_download.sh      helper for the PhysioNet datasets (pttppg, bpgraphene)
    inspect_pulsedb.py         prints the structure of one PulseDB .mat file (to check the field names)
notebooks/                     analyses; 01_example_access.ipynb shows the loading functions
```

## Datasets

| ID           | Source                                     | Synthetic | Records                       | Signals                    | BP                        | Status |
|--------------|--------------------------------------------|-----------|-------------------------------|----------------------------|---------------------------|--------|
| `ppgbp`      | Liang et al. 2018 (figshare)               | no        | 2.1 s segments, 3 per subject | ppg (1000 Hz)              | 1 cuff reading / subject  | raw data not on this machine |
| `pttppg`     | PhysioNet PTT-PPG v1.1.0                   | no        | sit / walk / run per subject  | ecg, ppg_1..6 (500 Hz)     | cuff at start and end     | raw data not on this machine |
| `bpgraphene` | PhysioNet graphene bioimpedance v1.0.0     | no        | one per trial                 | ppg, finapres_ppg, finapres_bp | TODO (from finapres_bp) | raw data not on this machine |
| `mimicbp`    | Sanches et al. 2024 (Harvard Dataverse)    | no        | 30 s segments, 30 per subject | ppg, abp (125 Hz)          | median per segment        | raw data not on this machine |
| `pulsedb`    | Wang et al. 2023, v2_0                     | no        | one per consecutive clip      | ecg, ppg, abp (125 Hz)     | TODO (from abp)           | raw data not on this machine |
| `pwdb`       | Charlton et al. 2019, PWDB v0.2.0 (Zenodo 3275625) | **yes** | one simulated beat per virtual subject | abp, ppg, flow_velocity, area at 13 sites each (500 Hz) | brachial SBP/DBP of the simulated wave | converted once (test); raw data not on this machine |

## Shared schema

Defined in `src/ppg_bp_database/utils.py` and checked by `save_tables` on every conversion.

- **subjects**: `subject_id, dataset, synthetic, age, sex ("male"/"female"), height_cm, weight_kg`, plus dataset-specific extra columns
- **records**: `record_id, subject_id, dataset, synthetic, session, condition, duration_s, signals, sampling_rates_hz, wave_path`
- **bp**: `record_id, time_s, sbp, dbp, map, method`. `time_s` is on the record's time axis; NaN = the value describes the whole record.

`synthetic` (bool) is `True` for simulated data / virtual subjects (currently only `pwdb`) and `False` for measured data.
Every dataset module declares it as `SYNTHETIC` and passes it to `save_tables`. To use only real data:
`records[~records["synthetic"]]`.

IDs are globally unique (`<dataset>_...`), so tables of different datasets can simply be concatenated.

Each HDF5 signal has the attributes `fs` (Hz), `unit` and `t0` (s): sample *i* is at time `t0 + i / fs`.
Gaps (e.g. removed segments in PulseDB) are NaN.

MAP is approximated as `dbp + (sbp - dbp) / 3`.

## Download

All datasets are open access; no login needed. On a new (Linux) machine:

```bash
scripts/download_all.sh                  # everything
scripts/download_all.sh ppgbp pwdb       # only some
```

Needs `curl`, `unzip`, `md5sum` / `sha1sum` / `sha256sum` and, for PulseDB, 7-Zip (Fedora: `sudo dnf install 7zip`).
Every file is checked against the checksums published by the source; a failed download stops the script for that
dataset with an error. Rerunning is safe (PhysioNet and MIMIC-BP skip complete files, PWDB and PulseDB resume).

| ID           | Script                   | Source                         | Raw size (approx.)                 |
|--------------|--------------------------|--------------------------------|------------------------------------|
| `ppgbp`      | `download_ppgbp.sh`      | figshare 10.6084/m9.figshare.5459299.v5 | 10 MB                     |
| `pttppg`     | `download_pttppg.sh`     | PhysioNet `pulse-transit-time-ppg` 1.1.0 | 2.9 GB                   |
| `bpgraphene` | `download_bpgraphene.sh` | PhysioNet `bp-graphene-bioimpedance` 1.0.0 | 2.3 GB (zip size; incl. unused bioimpedance) |
| `mimicbp`    | `download_mimicbp.sh`    | Harvard Dataverse doi:10.7910/DVN/DBM1NF (v2.2) | 0.5 GB zips + unpacked (`--all` adds ecg, resp: +0.4 GB) |
| `pwdb`       | `download_pwdb.sh`       | Zenodo 3275625 (v0.2.0)        | 0.4 GB (zip + unpacked)            |
| `pulsedb`    | `download_pulsedb.sh`    | Box links from the PulseDB GitHub (v2_0) | very large: archives + unpacked need twice the space |

Archives are kept after unpacking (except PPG-BP); delete them with `rm` once the conversion works.
The converted data needs space too: e.g. `pwdb` is ~1.1 GB in `data/structured`.

## Usage

```bash
uv sync
uv run python scripts/convert.py              # all datasets with raw data
uv run python scripts/convert.py ppgbp        # just one
uv run ppg-bp-database                        # overview of the converted datasets
```

In a notebook:

```python
import ppg_bp_database as db
records = db.load_records(["ppgbp", "pttppg"])        # one, several or all (no argument) datasets
waves = db.load_wave(records.iloc[0], start_s=0, end_s=5)   # {signal: pd.Series indexed by time in s}
```

The data folder can be moved by setting `PPG_BP_DATA_DIR`.

## Checking a download

- All `download_*.sh` scripts check the published checksums themselves (SHA-256 for PhysioNet, MD5 for figshare,
  Dataverse and Zenodo, SHA-1 for PulseDB).
- PhysioNet datasets can be rechecked later: `cd data/raw/pttppg && sha256sum -c --quiet SHA256SUMS.txt`
- All converters fail with a clear message if files are missing or have an unexpected shape or sampling rate.

## TODO

- Single SBP/DBP values from the ABP / Finapres waves (bpgraphene, pulsedb, later all datasets with an ABP wave)
- MIMIC-BP demographics from the MIMIC-III clinical database
- PulseDB: check the ECG unit
