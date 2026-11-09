"""Signal processing for one window: cut with padding, resample, filter, quality checks. Pure numpy / scipy."""
from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy.signal import butter, find_peaks, resample_poly, sosfiltfilt


class Rejected(Exception):
    """The window fails a quality check; the message is the reason (counted in report.csv)."""


def resample(x: np.ndarray, fs_in: float, fs_out: float) -> np.ndarray:
    """Polyphase resampling with anti-aliasing filter. padtype="line": the signal is continued linearly beyond its
    ends instead of with zeros, so a pressure wave (~80 mmHg) is not pulled towards 0 at the edges."""
    if fs_in == fs_out:
        return x
    ratio = Fraction(fs_out / fs_in).limit_denominator(1000)
    assert abs(float(ratio) - fs_out / fs_in) < 1e-9, f"no exact rational ratio for {fs_in} -> {fs_out} Hz"
    return resample_poly(x, ratio.numerator, ratio.denominator, padtype="line")


@lru_cache
def _bandpass_sos(fs: float, band: tuple[float, float]) -> np.ndarray:
    return butter(2, band, btype="band", fs=fs, output="sos")


def bandpass(x: np.ndarray, fs: float, band: tuple[float, float]) -> np.ndarray:
    """Zero-phase 2nd-order Butterworth band-pass (no time shift)."""
    return sosfiltfilt(_bandpass_sos(fs, tuple(band)), x)


def zscore(x: np.ndarray) -> np.ndarray:
    return (x - x.mean()) / x.std()


def peak_count(x: np.ndarray, fs: float, ratio: float) -> int:
    """Number of pulse peaks (at most one per 0.25 s = 240 bpm). A peak must reach `ratio` x the prominence of a
    typical peak (75th percentile), so dicrotic waves, which are much smaller than the systolic peak, are not
    counted. A fixed threshold would count them in pressure waves with a pronounced dicrotic wave."""
    peaks, props = find_peaks(x, distance=max(1, int(0.25 * fs)), prominence=0)
    if not len(peaks):
        return 0
    prominence = props["prominences"]
    return int((prominence >= ratio * np.percentile(prominence, 75)).sum())


def derivative_skew(x: np.ndarray) -> float:
    """Skewness of the first derivative. > 0 for a correctly oriented pulse wave (steep systolic upstroke,
    slow diastolic decay); < 0 hints at an inverted PPG."""
    d = np.diff(x)
    d = d - d.mean()
    return float((d ** 3).mean() / (d ** 2).mean() ** 1.5)


def cut(x: np.ndarray, fs_in: float, t0: float, start_s: float, n_out: int, fs_out: float, pad_s: float,
        name: str) -> tuple[np.ndarray, np.ndarray, float]:
    """Cut a window that starts at start_s (record time axis) and is n_out samples long at fs_out.
    Returns (padded segment resampled to fs_out, index of the window inside it, true time of window sample 0).
    The padding stops at NaN gaps and at the record ends. Raises Rejected if the window itself contains NaN."""
    i0 = int(np.ceil((start_s - t0) * fs_in - 1e-6))
    i1 = i0 + int(np.ceil(n_out * fs_in / fs_out - 1e-6))
    if i1 > len(x):
        raise Rejected("incomplete")
    if not np.isfinite(x[i0:i1]).all():
        raise Rejected(f"{name}_nan")

    pad = int(round(pad_s * fs_in))
    a, b = max(0, i0 - pad), min(len(x), i1 + pad)
    nan_before = np.flatnonzero(~np.isfinite(x[a:i0]))
    nan_after = np.flatnonzero(~np.isfinite(x[i1:b]))
    if len(nan_before):
        a += nan_before[-1] + 1
    if len(nan_after):
        b = i1 + nan_after[0]

    y = resample(x[a:b], fs_in, fs_out)
    c = int(round((i0 - a) * fs_out / fs_in))
    if c + n_out > len(y):
        raise Rejected("incomplete")
    return y, c, t0 + a / fs_in + c / fs_out


def cut_periodic(x: np.ndarray, fs_in: float, fs_out: float, n_out: int) -> tuple[np.ndarray, int]:
    """One periodic beat (pwdb): tiled 3x, so resampling and filtering see a continuous signal; the window is the
    middle copy. Returns (tiled resampled signal, index of the window inside it)."""
    x = x[np.isfinite(x)]  # some pwdb channels store their last sample as invalid (NaN)
    y = resample(np.tile(x, 3), fs_in, fs_out)
    c = int(round(len(x) * fs_out / fs_in))
    return y, c


def process_ppg(y: np.ndarray, c: int, n: int, raw: np.ndarray, fs: float, band: tuple[float, float],
                inverted: bool, min_std: float, max_clip: float) -> np.ndarray:
    """Band-pass the padded, resampled PPG, crop the window, check and z-score it.
    raw: the window at its original rate (for the clipping check)."""
    if inverted:
        y = -y
    w = bandpass(y, fs, band)[c:c + n]
    if w.std() < min_std:
        raise Rejected("ppg_flat")
    at_limit = (raw == raw.max()) | (raw == raw.min())
    if at_limit.mean() > max_clip:
        raise Rejected("ppg_clipped")
    return zscore(w)
