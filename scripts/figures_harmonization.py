"""Figures for harmonization_description.md, saved to figures/. Uses only the structured data (data/structured), so
the figures show the signals before harmonization.

    uv run python scripts/figures_harmonization.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ppg_bp_database.harmonize.pipeline import _read_signal
from ppg_bp_database.harmonize.signals import bandpass, derivative_skew
from ppg_bp_database.utils import PROJECT_DIR, STRUCTURED_DIR

FIG_DIR = PROJECT_DIR / "figures"
WINDOW_S = 10.0
BAND = (0.5, 15.0)


def _windows(x: np.ndarray, fs: float) -> list[slice]:
    n = int(WINDOW_S * fs)
    return [slice(i, i + n) for i in range(0, len(x) - n + 1, n)]


def dskew_example(show_s: float = 5.0) -> Path:
    """Same 10 s of a seated pttppg record at the native 500 Hz, the used distal channel vs proximal ppg_4. Raw signal,
    only flipped (pttppg is stored upside down, step 2). dskew is calculated on the 0.5-15 Hz band-passed signal.
    Picked: the window where ppg_4 has the most negative band-passed dskew while the used channel is clearly positive (> 1)."""
    from ppg_bp_database.harmonize.config import HarmonizeConfig
    used = HarmonizeConfig().pttppg_ppg
    records = pd.read_parquet(STRUCTURED_DIR / "pttppg" / "records.parquet")
    best = None
    for _, rec in records[records["condition"] == "sit"].iterrows():
        path = STRUCTURED_DIR / "pttppg" / rec["wave_path"]
        (g, fs, _), (b, _, _) = _read_signal(path, used), _read_signal(path, "ppg_4")
        g, b = -g, -b
        gf, bf = bandpass(g, fs, BAND), bandpass(b, fs, BAND)
        for w in _windows(g, fs):
            sg, sb = derivative_skew(gf[w]), derivative_skew(bf[w])
            if sg > 1 and (best is None or sb < best["sb"]):
                best = dict(sb=sb, sg=sg, record_id=rec["record_id"], fs=fs, w=w, g=g[w], b=b[w])
    fs, record_id = best["fs"], best["record_id"]
    print(f"dskew example: {record_id} starting at {best['w'].start / fs:.0f} s")

    n = int(show_s * fs)
    t = np.arange(n) / fs
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 4.2))
    for ax, x, name, skew, color in ((axes[0], best["g"], f"{used} (distal phalanx)", best["sg"], "C0"),
                                     (axes[1], best["b"], "ppg_4 (proximal phalanx)", best["sb"], "C3")):
        ax.plot(t, x[:n] - x[:n].mean(), lw=0.8, color=color)
        ax.set_title(f"{name}: dskew (0.5-15 Hz) = {skew:+.2f}", fontsize=10, loc="left")
        ax.set_ylabel("PPG (a.u., flipped)")
    axes[-1].set_xlabel("time (s)")
    axes[-1].set_xlim(0, show_s)
    fig.suptitle(f"pttppg {record_id} (seated), raw at {fs:.0f} Hz, from {best['w'].start / fs:.0f} s", fontsize=10)
    fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    path = FIG_DIR / "dskew_pttppg_raw.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def finapres_calibration_example(record_id: str = "bpgraphene_s3_day1_setup08_trial33", t1: float = 18.0,
                                  t2: float = 42.0) -> Path:
    """Raw Finapres pressure wave and PPG of one bpgraphene record during a Finapres calibration (Physiocal), with
    the stretches excluded from the usable segments (step 3)."""
    from ppg_bp_database.harmonize.config import HarmonizeConfig
    from ppg_bp_database.harmonize.segments import usable_segments

    config = HarmonizeConfig()
    records = pd.read_parquet(STRUCTURED_DIR / "bpgraphene" / "records.parquet").set_index("record_id")
    path = STRUCTURED_DIR / "bpgraphene" / records.loc[record_id, "wave_path"]
    signals = {"ppg": _read_signal(path, "ppg"), "abp": _read_signal(path, "finapres_bp")}
    segments = usable_segments(signals, {"abp"}, config)

    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 4.2))
    for ax, key, label, color in ((axes[0], "abp", "finapres_bp (mmHg)", "C1"), (axes[1], "ppg", "ppg (a.u.)", "C0")):
        x, fs, t0 = signals[key]
        t = t0 + np.arange(len(x)) / fs
        m = (t >= t1) & (t < t2)
        ax.plot(t[m], x[m], lw=0.8, color=color)
        ax.set_ylabel(label)
        excluded = [(b, a) for (_, b), (a, _) in zip(segments[:-1], segments[1:])]
        for k, (a, b) in enumerate(excluded):
            ax.axvspan(a, b, color="0.6", alpha=0.25, lw=0, label="excluded (flat + 1 s margin)" if k == 0 else None)
    axes[0].legend(loc="upper left", fontsize=8)
    axes[-1].set_xlabel("time in record (s)")
    axes[-1].set_xlim(t1, t2)
    fig.suptitle(f"bpgraphene {record_id}: Finapres calibration, raw signals", fontsize=10)
    fig.tight_layout()
    out = FIG_DIR / "finapres_calibration.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _pttppg_piece(record_id: str, channel: str, t_start: float) -> tuple[np.ndarray, int, int]:
    """Flipped pttppg channel resampled to the target rate: (window with config.pad_s padding, n, pad)."""
    from ppg_bp_database.harmonize.config import HarmonizeConfig
    from ppg_bp_database.harmonize.signals import resample

    config = HarmonizeConfig()
    fs_out, n, pad = config.fs, int(config.window_s * config.fs), int(config.pad_s * config.fs)
    records = pd.read_parquet(STRUCTURED_DIR / "pttppg" / "records.parquet").set_index("record_id")
    x, fs, _ = _read_signal(STRUCTURED_DIR / "pttppg" / records.loc[record_id, "wave_path"], channel)
    i = int(t_start * fs_out)
    return resample(-x, fs, fs_out)[i - pad:i + n + pad], n, pad


def _plot_rows(rows: list[tuple[np.ndarray, str, str]], title: str, name: str) -> Path:
    from ppg_bp_database.harmonize.config import HarmonizeConfig

    config = HarmonizeConfig()
    t = np.arange(len(rows[0][0])) / config.fs
    fig, axes = plt.subplots(len(rows), 1, sharex=True, figsize=(9, 2.1 * len(rows)))
    for ax, (y, label, color) in zip(axes, rows):
        ax.plot(t, y, lw=0.8, color=color)
        ax.set_title(label, fontsize=10, loc="left")
        ax.set_ylabel("PPG (a.u.)")
    axes[-1].set_xlabel("time in window (s)")
    axes[-1].set_xlim(0, config.window_s)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    out = FIG_DIR / name
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def channel_shape_example(record_id: str = "pttppg_s21_sit", t_start: float = 72.0) -> Path:
    """Same 10 s window of pttppg ppg_2 and ppg_3 after the band-pass (step 1): ppg_3 has the higher dskew but a broad top
    with a second peak almost as high as the systolic one; ppg_2 shows the classic pulse shape. The window is typical
    for ppg_3 (its median beat is closer to the channel's median beat than 62 % of the seated windows)."""
    from ppg_bp_database.harmonize.config import HarmonizeConfig

    config = HarmonizeConfig()
    rows = []
    for channel, color in (("ppg_2", "C0"), ("ppg_3", "C1")):
        y, n, pad = _pttppg_piece(record_id, channel, t_start)
        rows.append((bandpass(y, config.fs, config.ppg_band_hz)[pad:pad + n], f"{channel}, band-pass 0.5-15 Hz", color))
    return _plot_rows(rows, f"{record_id} (seated), window from {t_start:.0f} s", "pttppg_channels.png")


def bandpass_example(record_id: str = "pttppg_s12_sit", t_start: float = 312.0) -> Path:
    """One 10 s window of the used pttppg channel before and after the band-pass (step 9): flipped, resampled to
    125 Hz, then filtered with padding like in the pipeline. Picked from seated pttppg windows with 30-70 % of the
    power below 0.5 Hz as one with a very regular pulse after filtering."""
    from ppg_bp_database.harmonize.config import HarmonizeConfig

    config = HarmonizeConfig()
    y, n, pad = _pttppg_piece(record_id, config.pttppg_ppg, t_start)
    rows = [(y[pad:pad + n], "before: flipped and resampled to 125 Hz", "C3"),
            (bandpass(y, config.fs, config.ppg_band_hz)[pad:pad + n], "after: band-pass 0.5-15 Hz", "C0")]
    return _plot_rows(rows, f"{record_id} {config.pttppg_ppg} (seated), window from {t_start:.0f} s",
                      "bandpass_example.png")

if __name__ == "__main__":
    print(dskew_example())
    print(finapres_calibration_example())
    print(channel_shape_example())
    print(bandpass_example())
