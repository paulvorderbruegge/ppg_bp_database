"""Usable time segments of a record: the stretches where all used signals contain valid signal. Windows are only
placed inside these segments (windows.py), so they never cover a gap or a device calibration.

Invalid are
    NaN                       in any used signal (record ends, PulseDB's removed segments)
    flat pressure wave        only for sources with abp_flat_check (Finapres): its calibration (Physiocal) holds
                              the pressure at constant levels for ~0.5-1 s; widened by flat_margin_s on both sides,
                              because the beats right next to a calibration step are distorted
All times are on the record's time axis (sample i at t0 + i / fs).
"""
import numpy as np
from scipy.ndimage import maximum_filter1d, minimum_filter1d

from ppg_bp_database.harmonize.config import HarmonizeConfig

Interval = tuple[float, float]


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """[start, end) sample indices of the True runs of mask."""
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def flat_mask(x: np.ndarray, fs: float, duration_s: float, max_range: float) -> np.ndarray:
    """True for every sample inside a stretch of at least duration_s whose values stay within max_range."""
    n = max(int(round(duration_s * fs)), 2)
    if len(x) < n:
        return np.zeros(len(x), bool)
    hi = maximum_filter1d(np.where(np.isfinite(x), x, np.inf), n, origin=-(n // 2))   # max of x[i:i + n]
    lo = minimum_filter1d(np.where(np.isfinite(x), x, -np.inf), n, origin=-(n // 2))
    starts = (hi - lo)[: len(x) - n + 1] < max_range
    return np.convolve(starts, np.ones(n), "full")[: len(x)] > 0


def _span(x: np.ndarray, fs: float, t0: float) -> Interval:
    """Time span from the first to the last finite sample."""
    finite = np.flatnonzero(np.isfinite(x))
    if not len(finite):
        return np.inf, -np.inf
    return t0 + finite[0] / fs, t0 + (finite[-1] + 1) / fs


def usable_segments(signals: dict[str, tuple[np.ndarray, float, float]], flat_check: set[str],
                    config: HarmonizeConfig) -> list[Interval]:
    """signals: name -> (data, fs, t0). flat_check: names of the signals checked for flat stretches.
    Returns the sorted, non-overlapping segments (begin, end); empty if the signals do not overlap."""
    spans = [_span(*s) for s in signals.values()]
    begin, end = max(s[0] for s in spans), min(s[1] for s in spans)
    if end <= begin:
        return []
    invalid: list[Interval] = []
    for name, (x, fs, t0) in signals.items():
        invalid += [(t0 + a / fs, t0 + b / fs) for a, b in _runs(~np.isfinite(x))]
        if name in flat_check:
            m = config.flat_margin_s
            invalid += [(t0 + a / fs - m, t0 + b / fs + m)
                        for a, b in _runs(flat_mask(x, fs, config.flat_abp_s, config.flat_abp_range))]

    segments, t = [], begin
    for a, b in sorted(invalid):
        if a > t:
            segments.append((t, min(a, end)))
        t = max(t, b)
        if t >= end:
            break
    if t < end:
        segments.append((t, end))
    return [(a, b) for a, b in segments if b > a]
