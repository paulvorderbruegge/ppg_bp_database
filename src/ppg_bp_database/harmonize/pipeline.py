"""Harmonization: data/structured/<dataset>/ -> data/harmonized/<config name>/ (see README.md, "Harmonization").

Output:
    config.json        the HarmonizeConfig used
    subjects.parquet   subjects that have at least one window (SUBJECT_COLUMNS + dataset-specific extras)
    windows.parquet    one row per window (WINDOW_COLUMNS), no BP values
    bp.parquet         BP rows per window (WINDOW_BP_COLUMNS); time_s relative to the window start
    waves/<dataset>.h5 "ppg", "abp": float32 (n_windows, window_s * fs); row i = windows.row i of that dataset
    report.csv         windows kept / rejected per dataset and reason

Everything is streamed record by record, so PulseDB (millions of windows, tens of millions of beats) never has to be
in memory as a whole.
"""
import shutil
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ppg_bp_database import access
from ppg_bp_database.harmonize.config import HarmonizeConfig
from ppg_bp_database.harmonize.signals import Rejected, cut, cut_periodic, derivative_skew, peak_count, process_ppg
from ppg_bp_database.harmonize.sources import Source, sources
from ppg_bp_database.harmonize.windows import plan_record
from ppg_bp_database.utils import BP_COLUMNS, HARMONIZED_DIR, STRUCTURED_DIR, _conform

WINDOW_COLUMNS = {
    "window_id": "string",
    "dataset": "string",
    "subject_id": "string",
    "record_id": "string",
    "synthetic": "bool",
    "condition": "string",
    "row": "int64",           # row in waves/<dataset>.h5
    "t_start_s": "float64",   # time of window sample 0 on the record's time axis
    "n_samples": "int64",     # valid samples; the arrays are NaN-padded after them
    "duration_s": "float64",  # n_samples / fs
    "full_length": "bool",    # False = record shorter than window_s, kept at its own length
    "has_abp": "bool",        # False = abp array is all NaN (dataset has no pressure wave)
    "ppg_signal": "string",   # structured signal the PPG was taken from
    "ppg_dskew": "float64",   # skewness of the PPG derivative; > 0 = correctly oriented pulse
}
WINDOW_BP_COLUMNS = {
    "window_id": "string",
    "dataset": "string",
    **{c: t for c, t in BP_COLUMNS.items() if c != "record_id"},  # time_s now relative to the window start
}
BATCH = 512  # windows per write


class _ParquetAppender:
    def __init__(self, path: Path, columns: dict[str, str]):
        self.columns = columns
        empty = _conform(pd.DataFrame(columns=list(columns)), columns, path.stem, extra_allowed=False)
        self.schema = pa.Schema.from_pandas(empty, preserve_index=False)
        self.writer = pq.ParquetWriter(path, self.schema)

    def write(self, rows: list[dict] | pd.DataFrame) -> None:
        df = pd.DataFrame(rows) if isinstance(rows, list) else rows
        if len(df):
            df = _conform(df, self.columns, "harmonized", extra_allowed=False)
            self.writer.write_table(pa.Table.from_pandas(df, schema=self.schema, preserve_index=False))

    def close(self) -> None:
        self.writer.close()


class _WaveAppender:
    def __init__(self, path: Path, n: int, fs: float):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.file = h5py.File(path, "w")
        self.rows = 0
        for name in ("ppg", "abp"):
            ds = self.file.create_dataset(name, shape=(0, n), maxshape=(None, n), dtype="float32",
                                          chunks=(64, n), compression="gzip")
            ds.attrs["fs"] = fs
        self.file["ppg"].attrs["unit"] = "z-score"
        self.file["abp"].attrs["unit"] = "mmHg"

    def write(self, ppg: list[np.ndarray], abp: list[np.ndarray]) -> None:
        if not ppg:
            return
        for name, arrays in (("ppg", ppg), ("abp", abp)):
            ds = self.file[name]
            ds.resize(self.rows + len(arrays), axis=0)
            ds[self.rows:] = np.stack(arrays)
        self.rows += len(ppg)

    def close(self) -> None:
        self.file.close()


def _record_bp(dataset_dir: Path) -> Iterator[tuple[str, pd.DataFrame]]:
    """(record_id, its BP rows) for every record in bp.parquet, read in batches. The rows of one record are
    contiguous (sorted, or written per file for PulseDB)."""
    path = dataset_dir / "bp.parquet"
    if not path.is_file():
        return
    seen: set[str] = set()
    carry = None
    for batch in pq.ParquetFile(path).iter_batches(batch_size=500_000):
        df = batch.to_pandas()
        if carry is not None:
            df = pd.concat([carry, df], ignore_index=True)
        last = df["record_id"].iloc[-1]
        carry = df[df["record_id"] == last]
        for record_id, rows in df[df["record_id"] != last].groupby("record_id", sort=False):
            assert record_id not in seen, f"{path}: rows of {record_id} are not contiguous"
            seen.add(record_id)
            yield record_id, rows
    if carry is not None and len(carry):
        yield carry["record_id"].iloc[0], carry


def _read_signal(path: Path, name: str) -> tuple[np.ndarray, float, float]:
    with h5py.File(path, "r") as f:
        ds = f[name]
        return ds[:].astype(float), float(ds.attrs["fs"]), float(ds.attrs["t0"])  # type: ignore


def _span(x: np.ndarray, fs: float, t0: float) -> tuple[float, float]:
    """Time span from the first to the last finite sample."""
    finite = np.flatnonzero(np.isfinite(x))
    if not len(finite):
        return np.inf, -np.inf
    return t0 + finite[0] / fs, t0 + (finite[-1] + 1) / fs


def _window(signals: dict, src: Source, start: float, n: int, config: HarmonizeConfig
            ) -> tuple[np.ndarray, np.ndarray | None, float]:
    """(z-scored PPG, ABP in mmHg or None, true start time) of one window. Raises Rejected."""
    fs = config.fs
    x, fs_in, t0 = signals["ppg"]
    if src.periodic:
        y, c = cut_periodic(x, fs_in, fs, n)
        raw, t_start = x[np.isfinite(x)], t0
    else:
        y, c, t_start = cut(x, fs_in, t0, start, n, fs, config.pad_s, "ppg")
        i0 = int(np.ceil((start - t0) * fs_in - 1e-6))
        raw = x[i0:i0 + int(np.ceil(n * fs_in / fs - 1e-6))]
    ppg = process_ppg(y, c, n, raw, fs, config.ppg_band_hz, src.ppg_inverted, config.min_ppg_std,
                      config.max_clip_fraction)

    abp = None
    if "abp" in signals:
        x, fs_in, t0 = signals["abp"]
        y, c = cut_periodic(x, fs_in, fs, n) if src.periodic else cut(x, fs_in, t0, start, n, fs, config.pad_s, "abp")[:2]
        abp = y[c:c + n]
        diff = abs(peak_count(ppg, fs, config.peak_prominence_ratio) - peak_count(abp, fs, config.peak_prominence_ratio))
        if diff > config.max_beat_count_diff:
            raise Rejected("hr_mismatch")
    return ppg, abp, t_start


def harmonize_dataset(name: str, config: HarmonizeConfig, out_dir: Path, windows_out: _ParquetAppender,
                      bp_out: _ParquetAppender, structured_dir: Path = STRUCTURED_DIR) -> tuple[Counter, set[str]]:
    """Harmonize one dataset: writes waves/<name>.h5 and appends to windows / bp. Returns (report counts,
    subject IDs with windows)."""
    src = sources(config)[name]
    records = access.load_records(name, structured_dir).set_index("record_id")
    n_full = int(round(config.window_s * config.fs))
    waves = _WaveAppender(out_dir / "waves" / f"{name}.h5", n_full, config.fs)
    report: Counter = Counter()
    subjects: set[str] = set()
    win_rows, bp_rows, ppg_buf, abp_buf = [], [], [], []

    def flush():
        waves.write(ppg_buf, abp_buf)
        windows_out.write(win_rows)
        bp_out.write(pd.concat(bp_rows, ignore_index=True) if bp_rows else pd.DataFrame())
        for buf in (win_rows, bp_rows, ppg_buf, abp_buf):
            buf.clear()

    done = set()
    for record_id, rec_bp in _record_bp(structured_dir / name):
        rec_bp = rec_bp.reset_index(drop=True)
        done.add(record_id)
        rec = records.loc[record_id]
        wanted = {"ppg": src.ppg, **({"abp": src.abp} if src.abp else {})}
        missing = set(wanted.values()) - set(rec["signals"])
        assert not missing, f"{record_id}: signals {sorted(missing)} not found, has {list(rec['signals'])}"
        path = structured_dir / name / rec["wave_path"]
        signals = {key: _read_signal(path, sig) for key, sig in wanted.items()}

        spans = [_span(*s) for s in signals.values()]
        t_begin, t_end = max(s[0] for s in spans), min(s[1] for s in spans)
        if t_end <= t_begin:
            report["no_signal"] += 1
            continue
        plan = plan_record(t_begin, t_end, rec_bp, config)
        report.update(plan.rejected)
        n = n_full if plan.full_length else int(np.floor(plan.duration_s * config.fs + 1e-6))

        kept_pos, kept_ids, kept_starts = [], [], []
        for k, (start, pos) in enumerate(zip(plan.starts, plan.bp)):
            try:
                ppg, abp, t_start = _window(signals, src, start, n, config)
            except Rejected as e:
                report[str(e)] += 1
                continue
            window_id = f"{record_id}_w{k:04d}"
            win_rows.append({
                "window_id": window_id, "dataset": name, "subject_id": rec["subject_id"], "record_id": record_id,
                "synthetic": rec["synthetic"], "condition": rec["condition"], "row": waves.rows + len(ppg_buf),
                "t_start_s": t_start, "n_samples": n, "duration_s": n / config.fs, "full_length": plan.full_length,
                "has_abp": abp is not None, "ppg_signal": src.ppg, "ppg_dskew": derivative_skew(ppg),
            })
            kept_pos.append(pos)
            kept_ids.append(window_id)
            kept_starts.append(t_start)
            ppg_buf.append(np.pad(ppg, (0, n_full - n), constant_values=np.nan))
            abp_buf.append(np.full(n_full, np.nan) if abp is None else np.pad(abp, (0, n_full - n), constant_values=np.nan))
            subjects.add(rec["subject_id"])
            report["kept"] += 1
        if kept_pos:  # BP rows of all kept windows of this record at once (per-window pandas ops are slow)
            counts = [len(p) for p in kept_pos]
            rows = rec_bp.iloc[np.concatenate(kept_pos)].drop(columns="record_id").reset_index(drop=True)
            rows["window_id"] = np.repeat(kept_ids, counts)
            rows["dataset"] = name
            rows["time_s"] = rows["time_s"].to_numpy() - np.repeat(kept_starts, counts)
            bp_rows.append(rows)
        if len(ppg_buf) >= BATCH:
            flush()
    flush()
    waves.close()
    if without_bp := len(set(records.index) - done):
        report["record_without_bp"] += without_bp
    return report, subjects


def harmonize(datasets: list[str] | None = None, config: HarmonizeConfig = HarmonizeConfig(),
              structured_dir: Path = STRUCTURED_DIR, harmonized_dir: Path = HARMONIZED_DIR) -> Path:
    """Harmonize the given datasets (default: all converted ones) into harmonized_dir/<config.name>/.
    The folder is rebuilt from scratch and then contains exactly these datasets."""
    names = datasets or access.available_datasets(structured_dir)
    unknown = set(names) - set(sources(config))
    assert not unknown, f"no harmonization source defined for {sorted(unknown)} (see harmonize/sources.py)"
    out_dir = harmonized_dir / config.name
    shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.mkdir(parents=True)
    config.save(out_dir / "config.json")

    windows_out = _ParquetAppender(out_dir / "windows.parquet", WINDOW_COLUMNS)
    bp_out = _ParquetAppender(out_dir / "bp.parquet", WINDOW_BP_COLUMNS)
    report_rows, subject_ids = [], set()
    try:
        for name in names:
            report, subjects = harmonize_dataset(name, config, out_dir, windows_out, bp_out, structured_dir)
            subject_ids |= subjects
            report_rows += [{"dataset": name, "reason": r, "count": c} for r, c in sorted(report.items())]
            print(f"{name}: {report['kept']} windows kept, rejected: "
                  f"{ {r: c for r, c in report.items() if r != 'kept'} }")
    finally:
        windows_out.close()
        bp_out.close()

    subjects = access.load_subjects(names, structured_dir)
    subjects[subjects["subject_id"].isin(subject_ids)].to_parquet(out_dir / "subjects.parquet", index=False)
    pd.DataFrame(report_rows).to_csv(out_dir / "report.csv", index=False)
    print(f"-> {out_dir}")
    return out_dir
