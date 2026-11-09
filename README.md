# ppg-bp-database

Converts several PPG / blood-pressure datasets into one shared format, so that notebooks can use them all the same way.

## Layout

```
data/raw/<dataset>/            downloaded data, untouched
data/structured/<dataset>/     structured data (written by scripts/structure.py)
    subjects.parquet           one row per subject
    records.parquet            one row per recording
    bp.parquet                 one row per BP value (missing if the dataset has none yet)
    waves/<record_id>.h5       one HDF5 file per record, one HDF5 dataset per signal
src/ppg_bp_database/
    structure/                 step 1: raw -> data/structured
        <dataset>.py           one converter per dataset, all with the same functions
        beats.py               beat-to-beat SBP / DBP / MAP from a continuous BP wave
    harmonize/                 step 2: structured -> data/harmonized (config, sources, windows, signals, pipeline)
    utils.py                   shared schema, paths, save_tables, write_wave / read_wave
    access.py                  loading functions for notebooks (structured and harmonized data)
scripts/
    structure.py               step 1, for one, several or all datasets
    harmonize.py               step 2 (data/harmonized/<config>/)
    download_all.sh            runs the download_<dataset>.sh scripts below (all, or the ones given)
    download_<dataset>.sh      download + checksum check + unpack of one dataset into data/raw/<dataset>
    physionet_download.sh      helper for the PhysioNet datasets (pttppg, bpgraphene)
    download_sample.py         a few records of every dataset (~22 MB) to data/raw, for a quick end-to-end test
    inspect_pulsedb.py         prints the structure of one PulseDB .mat file (to check the field names)
    resampling_error.py        round-trip error native -> harmonized rate -> native, per dataset
notebooks/                     analyses; 01_example_access.ipynb shows the loading functions,
                               02_harmonized_check.ipynb checks the harmonized data
tests/                         uv run pytest (harmonization on synthetic data)
```

## Datasets

| ID           | Source                                     | Synthetic | Records                       | Signals                    | BP                        | Status |
|--------------|--------------------------------------------|-----------|-------------------------------|----------------------------|---------------------------|--------|
| `ppgbp`      | Liang et al. 2018 (figshare)               | no        | 2.1 s segments, 3 per subject | ppg (1000 Hz)              | 1 cuff reading / subject  | sample only (`download_sample.py`) |
| `pttppg`     | PhysioNet PTT-PPG v1.1.0                   | no        | sit / walk / run per subject  | ecg, ppg_1..6 (500 Hz)     | cuff at start and end     | sample only (`download_sample.py`) |
| `bpgraphene` | PhysioNet graphene bioimpedance v1.0.0     | no        | one per trial                 | ppg (1250 Hz), finapres_ppg (75 Hz), finapres_bp (200 Hz) | per beat, calculated from finapres_bp | sample only (`download_sample.py`) |
| `mimicbp`    | Sanches et al. 2024 (Harvard Dataverse)    | no        | 30 s segments, 30 per subject | ppg, abp (125 Hz)          | median per segment (dataset) + per beat, calculated from abp | sample only (`download_sample.py`) |
| `pulsedb`    | Wang et al. 2023, v2_0                     | no        | one per consecutive clip      | ecg, ppg, abp (125 Hz)     | per beat, calculated from abp | sample only (`download_sample.py`) |
| `pwdb`       | Charlton et al. 2019, PWDB v0.2.0 (Zenodo 3275625) | **yes** | one simulated beat per virtual subject | abp, ppg, flow_velocity, area at 13 sites each (500 Hz) | brachial SBP/DBP of the simulated wave | sample only (`download_sample.py`) |

## Shared schema

Defined in `src/ppg_bp_database/utils.py` and checked by `save_tables` on every conversion.

- **subjects**: `subject_id, dataset, synthetic, age, sex ("male"/"female"), height_cm, weight_kg`, plus dataset-specific extra columns
- **records**: `record_id, subject_id, dataset, synthetic, session, condition, duration_s, signals, sampling_rates_hz, wave_path`
- **bp**: `record_id, time_s, sbp, dbp, map, method, source`. `time_s` is on the record's time axis; NaN = the value
  describes the whole record. `method` = what measured the pressure (`cuff`, `abp` = arterial line wave, `finapres`).
  `source` = how we got the value: `dataset` (given by the dataset) or `calculated` (calculated by us from a wave).

Calculated values are one row per heartbeat (`ppg_bp_database.structure.beats`), computed during structuring from
the wave at its original sampling rate (before any resampling): `time_s` = beat onset, `sbp` = systolic peak, `dbp` = pressure
at the onset, `map` = mean of the wave over the beat. Only checks that separate beats from non-beats are applied
(beat length 0.25-2 s, pulse pressure >= 5 mmHg, no NaN gap inside the beat); quality filtering belongs to the
harmonization step. To compare, e.g., MIMIC-BP's own segment medians with ours:
`bp.groupby(["record_id", "source"])["sbp"].median().unstack()`.

`synthetic` (bool) is `True` for simulated data / virtual subjects (currently only `pwdb`) and `False` for measured data.
Every dataset module declares it as `SYNTHETIC` and passes it to `save_tables`. To use only real data:
`records[~records["synthetic"]]`.

IDs are globally unique (`<dataset>_...`), so tables of different datasets can simply be concatenated.

Each HDF5 signal has the attributes `fs` (Hz), `unit` and `t0` (s): sample *i* is at time `t0 + i / fs`.
Gaps (e.g. removed segments in PulseDB) are NaN.

MAP is approximated as `dbp + (sbp - dbp) / 3` for values given by a dataset without MAP; calculated values use the
exact mean of the wave.

## Harmonization

`data/structured` keeps every dataset as recorded. `scripts/harmonize.py` turns them into one training-ready
dataset with the same sampling rate, fixed-length windows and the same PPG processing:

```bash
uv run python scripts/harmonize.py                  # all converted datasets, config "default"
uv run python scripts/harmonize.py ppgbp pttppg     # only these (the output then contains only these)
```

```
data/harmonized/<config>/
    config.json          parameters used (harmonize/config.py: CONFIGS)
    subjects.parquet     subjects with at least one window
    windows.parquet      one row per window, no BP values
    bp.parquet           BP rows per window: window_id, dataset, time_s, sbp, dbp, map, method, source
    waves/<dataset>.h5   "ppg" (z-scored), "abp" (mmHg): float32 (n_windows, window_s * fs); row = windows.row
    report.csv           windows kept / rejected per dataset and reason
```

```python
windows = db.load_windows()                       # filter e.g. windows[windows.full_length & windows.has_abp]
bp = db.load_window_bp()                          # join on window_id
arrays = db.load_window_waves(windows.iloc[:64])  # {"ppg": (64, 1250), "abp": (64, 1250)}
```

**Rate: 125 Hz.** PPG content, including the dicrotic notch and second-derivative features, lies below ~15 Hz.
125 Hz is the native rate of the two largest datasets (PulseDB, MIMIC-BP), so a higher rate would only
interpolate them. All other datasets are downsampled (`resample_poly`, anti-aliasing filter), never upsampled.
`scripts/resampling_error.py` shows what the downsampling loses (round trip native -> 125 Hz -> native).

**Per window:**
- **PPG:** cut with 1 s padding, polarity corrected, resampled, band-passed at 0.5-15 Hz (zero-phase), cropped,
  z-scored.
- **ABP:** only resampled, so it stays in mmHg.

**Signals:** `harmonize/sources.py`.
- ppgbp `ppg`
- pttppg `ppg_1` (configurable). It is inverted in the raw data and flipped here.
- bpgraphene `ppg` + `finapres_bp`
- mimicbp / pulsedb `ppg` + `abp`
- pwdb `ppg_digital` + `abp_brachial`

**Windows and BP rows.** These are derived from the kind of values in the structured bp table, with no
per-dataset code:
- **Calculated beats** (mimicbp, pulsedb, bpgraphene): windows tile the record. The BP rows are the valid beats
  whose onset lies in the window; `time_s` is relative to the window start. MIMIC-BP's own 30 s medians are not
  copied, because they don't describe a 10 s window; use them from `data/structured` for validation.
- **Timed readings** (pttppg cuff at start and end): only windows within `cuff_window_s` (60 s) of a reading.
  The BP row is the nearest reading, at its time relative to the window (it can lie outside the window).
- **Whole-record values** (pwdb; `time_s` NaN): windows tile the record, each with that value.
- **Records shorter than one window** (ppgbp 2.1 s, pwdb one beat) keep their own length: `full_length = False`,
  `n_samples` < window, and the arrays are NaN-padded. pwdb's periodic beat is tiled 3x for filtering only.

**Quality** (thresholds in the config; rejected windows are counted in `report.csv`):
- Beats outside SBP 60-250 / DBP 30-150 / pulse pressure >= 10 mmHg are dropped.
- Windows are rejected for:
  - fewer than 0.5 valid beats/s
  - SBP IQR > 15 mmHg
  - NaN (gaps)
  - flat PPG
  - > 5 % clipped PPG samples
  - PPG and ABP peak counts differing by more than 2
- `windows.ppg_dskew` (skewness of the PPG derivative) is > 0 for a correctly oriented pulse. Use it to check
  polarity per dataset.

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
uv run python scripts/structure.py              # all datasets with raw data
uv run python scripts/structure.py ppgbp        # just one
uv run ppg-bp-database                        # overview of the converted datasets
```

In a notebook:

```python
import ppg_bp_database as db
records = db.load_records(["ppgbp", "pttppg"])        # one, several or all (no argument) datasets
waves = db.load_wave(records.iloc[0], start_s=0, end_s=5)   # {signal: pd.Series indexed by time in s}
```

The data folder can be moved by setting `PPG_BP_DATA_DIR`.

## Quick test with a small sample

`scripts/download_sample.py` downloads a few records per dataset (~22 MB in total) to `data/raw/<dataset>` and marks
each folder with a `SAMPLE` file. It only replaces folders that are empty or carry this marker, so a full download is
never deleted. To get the full data, delete the sample folder and run `download_<dataset>.sh`. Files inside zips (including PulseDB's split zip) are read with HTTP range requests, so the large
archives are never downloaded. Each file is capped at 15 MB and the whole run at 100 MB.

```bash
uv run python scripts/download_sample.py
uv run python scripts/structure.py
uv run python scripts/harmonize.py
```

The sample is only for checking that the code runs, not for analysis. pttppg's `subjects_info.csv` is trimmed to the
downloaded records, and the PulseDB files are the smallest ones, each with a single 10 s segment.

## Checking a download

- All `download_*.sh` scripts check the published checksums themselves (SHA-256 for PhysioNet, MD5 for figshare,
  Dataverse and Zenodo, SHA-1 for PulseDB).
- PhysioNet datasets can be rechecked later: `cd data/raw/pttppg && sha256sum -c --quiet SHA256SUMS.txt`
- All converters fail with a clear message if files are missing or have an unexpected shape or sampling rate.

## TODO

- Validate the beat detection on real data (compare with the MIMIC-BP labels; look at Finapres calibration periods)
- Harmonization on real mimicbp / pulsedb / bpgraphene / pwdb: run notebooks/02_harmonized_check.ipynb
- PulseDB also gives its own SBP/DBP per 10 s segment (`SegSBP`, `SegDBP`); could be added as `source = dataset`
- PulseDB: field names checked only on single-segment files; check how mat73 loads a file with many segments
- MIMIC-BP demographics from the MIMIC-III clinical database
- PulseDB: check the ECG unit
