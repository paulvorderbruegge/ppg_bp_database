"""Synthetic arterial pressure / PPG waves with known per-beat SBP, DBP and MAP (for tests)."""
import numpy as np
from scipy.signal import butter, sosfiltfilt


def beat_shape(n: int) -> np.ndarray:
    """0 at onset, 1 at the systolic peak, dicrotic wave, decay back to 0."""
    p = np.linspace(0, 1, n, endpoint=False)
    sys = np.where(p < 0.12, np.sin(np.clip(p / 0.12, 0, 1) * np.pi / 2), np.exp(-((p - 0.12) / 0.22) ** 2))
    dic = 0.18 * np.exp(-((p - 0.42) / 0.05) ** 2)
    decay = 0.35 * np.exp(-(p - 0.42) / 0.25) * (p > 0.42) * (1 - p) / 0.58
    return np.maximum(sys, dic + decay * (p > 0.42) + sys * (p <= 0.42))


def abp(fs: float, duration_s: float, rng: np.random.Generator, noise: float = 0.3):
    """(wave in mmHg, truth array with one row per beat: onset time, sbp, dbp, map)."""
    x, truth, t = [], [], 0
    n_total = int(duration_s * fs)
    dbps = 72 + np.cumsum(rng.normal(0, 1, int(duration_s * 4) + 2)).clip(-15, 15)  # > 1 value per beat
    k = 0
    while t < n_total:
        n = int(round(fs * 60 / rng.uniform(55, 100)))
        sbp = dbps[k] + rng.uniform(35, 55)
        w = dbps[k] + (dbps[k + 1] - dbps[k]) * np.linspace(0, 1, n, endpoint=False) + (sbp - dbps[k]) * beat_shape(n)
        truth.append((t / fs, w.max(), w[0], w.mean()))
        x.append(w)
        t += n
        k += 1
    return np.concatenate(x)[:n_total] + rng.normal(0, noise, n_total), np.array(truth)


def ppg_from_abp(wave: np.ndarray, fs: float, rng: np.random.Generator, delay_s: float = 0.2) -> np.ndarray:
    """A PPG-like signal: delayed, smoothed pressure with arbitrary scale / offset, baseline wander and noise."""
    d = int(delay_s * fs)
    x = np.concatenate([np.full(d, wave[0]), wave[:-d]])
    x = sosfiltfilt(butter(2, 6, fs=fs, output="sos"), x)
    t = np.arange(len(x)) / fs
    return 3000 + 25 * x + 200 * np.sin(2 * np.pi * 0.1 * t) + rng.normal(0, 5, len(x))
