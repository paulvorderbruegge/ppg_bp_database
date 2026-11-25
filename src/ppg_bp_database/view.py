"""Notebook helpers for looking at the harmonized data: tables and window waves.

    from ppg_bp_database import view
    windows = view.load_table("windows")
    view.show(windows)                    # any DataFrame, e.g. view.show(windows[windows.has_abp])
    view.plot_window(windows.iloc[0])     # or a window_id string
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ppg_bp_database.access import (_harmonized, load_harmonized_subjects, load_window_bp, load_window_waves,
                                    load_windows)
from ppg_bp_database.utils import HARMONIZED_DIR

TABLES = ("windows", "bp", "subjects")


def load_table(name: str, config: str = "default", harmonized_dir: Path = HARMONIZED_DIR) -> pd.DataFrame:
    """One harmonized table: "windows", "bp" or "subjects"."""
    if name == "windows":
        return load_windows(config, harmonized_dir=harmonized_dir)
    if name == "bp":
        return load_window_bp(config, harmonized_dir=harmonized_dir)
    if name == "subjects":
        return load_harmonized_subjects(config, harmonized_dir=harmonized_dir)
    raise ValueError(f"unknown table {name!r}, expected one of {TABLES}")


def show(df: pd.DataFrame, n: int = 20, drop_empty: bool = True) -> None:
    """Display the first n rows with all columns visible. drop_empty hides columns that are empty in the whole
    table (e.g. dataset-specific subject columns)."""
    from IPython.display import display

    shown = df.dropna(axis=1, how="all") if drop_empty else df
    hidden = df.shape[1] - shown.shape[1]
    print(f"{len(df)} rows x {df.shape[1]} columns" + (f" ({hidden} empty columns hidden)" if hidden else ""))
    with pd.option_context("display.max_columns", None, "display.width", None, "display.max_colwidth", 60):
        display(shown.head(n))


def plot_window(window: pd.Series | str, config: str = "default", start_s: float | None = None,
                end_s: float | None = None, harmonized_dir: Path = HARMONIZED_DIR):
    """PPG and ABP of one window (a row of the windows table or a window_id) with its BP values. Returns the figure.
    start_s / end_s (s, relative to the window start) show only that part; they must lie inside the window.
    Calculated beats: x = beat onset (DBP), dashed line = SBP of that beat. Values given by the dataset (cuff
    readings, record-level values) are in the title, with their time relative to the window start if they have one."""
    import matplotlib.pyplot as plt

    if isinstance(window, str):
        windows = load_windows(config, harmonized_dir=harmonized_dir)
        match = windows[windows["window_id"] == window]
        if match.empty:
            raise KeyError(f"no window {window!r} in config {config!r}")
        window = match.iloc[0]
    fs = json.loads((_harmonized(config, harmonized_dir) / "config.json").read_text())["fs"]
    arrays = load_window_waves(pd.DataFrame([window]).infer_objects(), config, harmonized_dir)
    bp = load_window_bp(config, datasets=window["dataset"], harmonized_dir=harmonized_dir)
    bp = bp[bp["window_id"] == window["window_id"]]
    beats, given = bp[bp["source"] == "calculated"], bp[bp["source"] != "calculated"]

    n = int(window["n_samples"])
    length = n / fs
    start = 0.0 if start_s is None else start_s
    end = length if end_s is None else end_s
    if not 0 <= start < end <= length:
        raise ValueError(f"range {start:g}-{end:g} s is not inside window {window['window_id']} (0-{length:g} s)")
    i0, i1 = round(start * fs), min(round(end * fs) + 1, n)  # +1: include the sample at end_s
    t = np.arange(i0, i1) / fs
    beats = beats[(beats["time_s"] >= start) & (beats["time_s"] <= end)]
    has_abp = bool(window["has_abp"])
    fig, axes = plt.subplots(2 if has_abp else 1, 1, sharex=True, figsize=(12, 5 if has_abp else 2.8), squeeze=False)
    axes = axes[:, 0]
    axes[0].plot(t, arrays["ppg"][0, i0:i1], lw=1)
    axes[0].set_ylabel("PPG (z-score)")
    if has_abp:
        axes[1].plot(t, arrays["abp"][0, i0:i1], lw=1, color="C1")
        axes[1].plot(beats["time_s"], beats["dbp"], "kx", label="onset / DBP")
        axes[1].hlines(beats["sbp"], beats["time_s"], beats["time_s"] + 0.3, colors="k", linestyles="--", label="SBP")
        axes[1].set_ylabel("ABP (mmHg)")
        if not beats.empty:
            axes[1].legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("time in window (s)")
    axes[-1].set_xlim(start, end)

    part = "" if start_s is None and end_s is None else f" [{start:g}-{end:g} s]"
    title = f"{window['window_id']}{part}  |  {window['condition']}, {len(beats)} beats"
    for _, r in given.iterrows():
        at = "" if pd.isna(r["time_s"]) else f" at {r['time_s']:+.0f} s"
        title += f"  |  {r['method']}{at}: {r['sbp']:.0f}/{r['dbp']:.0f} mmHg"
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    return fig
