"""Signal processing for one window: cut with padding, resample, filter, normalize, quality checks. Pure numpy / scipy."""
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
    """Zero-phase 2nd-order Butterworth band-pass (no time shift). Beyond its ends the signal is extended by 1 s,
    mirrored in time (padtype="even"): where the window has no real padding (step 7), this halves the edge error
    compared to scipy's default (15 samples, padtype="odd")."""
    return sosfiltfilt(_bandpass_sos(fs, tuple(band)), x, padtype="even", padlen=min(int(round(fs)), len(x) - 1))


def zscore(x: np.ndarray) -> np.ndarray:
    sd = x.std()
    if sd < 1e-6:  # flat line, nothing to scale
        raise Rejected("ppg_flat")
    return (x - x.mean()) / sd


def peak_count(x: np.ndarray, fs: float, ratio: float, min_beats_per_s: float) -> int:
    """Number of pulse peaks (at most one per 0.25 s = 240 bpm). A peak must reach `ratio` x the prominence of a
    systolic peak, so dicrotic waves, which are much smaller than the systolic peak, are not counted. A fixed
    threshold would count them in pressure waves with a pronounced dicrotic wave.
    Reference systolic peak: the k-th most prominent peak, k = minimum number of beats in the window
    (min_beats_per_s x duration). As long as the heart rate is above that minimum, at least k peaks are systolic,
    so the reference is one too. A percentile of all candidate peaks fails at a slow heart rate: each beat then
    gives several small candidates (dicrotic wave, wiggles), and the percentile lands on one of them."""
    peaks, props = find_peaks(x, distance=max(1, int(0.25 * fs)), prominence=0)
    if not len(peaks):
        return 0
    prominence = np.sort(props["prominences"])[::-1]
    k = max(1, int(np.ceil(min_beats_per_s * len(x) / fs - 1e-6)))
    return int((prominence >= ratio * prominence[min(k, len(prominence)) - 1]).sum())


def derivative_skew(x: np.ndarray) -> float:
    """Skewness of the first derivative. > 0 for a correctly oriented pulse wave (steep systolic upstroke,
    slow diastolic decay); < 0 hints at an inverted PPG."""
    d = np.diff(x)
    d = d - d.mean()
    return float((d ** 3).mean() / (d ** 2).mean() ** 1.5)


def cut(x: np.ndarray, fs_in: float, t0: float, start_s: float, n_out: int, fs_out: float, pad_s: float
        ) -> tuple[np.ndarray, int, float]:
    """Cut a window that starts at start_s (record time axis) and is n_out samples long at fs_out.
    Returns (padded segment resampled to fs_out, index of the window inside it, true time of window sample 0).
    The padding stops at NaN gaps and at the record ends. The window itself lies inside a usable segment
    (segments.py)."""
    i0 = int(np.ceil((start_s - t0) * fs_in - 1e-6))
    i1 = i0 + int(np.ceil(n_out * fs_in / fs_out - 1e-6))

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
    return y, c, t0 + a / fs_in + c / fs_out


def cut_periodic(x: np.ndarray, fs_in: float, fs_out: float, n_out: int) -> tuple[np.ndarray, int]:
    """One periodic beat (pwdb): tiled 3x, so resampling and filtering see a continuous signal; the window is the
    middle copy. Returns (tiled resampled signal, index of the window inside it)."""
    x = x[np.isfinite(x)]  # some pwdb channels store their last sample as invalid (NaN)
    y = resample(np.tile(x, 3), fs_in, fs_out)
    c = int(round(len(x) * fs_out / fs_in))
    return y, c


def process_ppg(y: np.ndarray, c: int, n: int, fs: float, band: tuple[float, float] | None) -> np.ndarray:
    """Band-pass the padded, resampled PPG (skipped if band is None), crop the window and z-score it.
    Orientation is already corrected (pipeline)."""
    return zscore((y if band is None else bandpass(y, fs, band))[c:c + n])
