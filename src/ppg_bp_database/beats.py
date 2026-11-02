"""Beat-to-beat SBP / DBP / MAP from a continuous blood pressure wave (arterial line or Finapres).

Beat k runs from its onset (the pressure minimum before its systolic peak) to the onset of beat k + 1:

    time_s  onset of the beat on the record's time axis
    sbp     systolic peak of the beat
    dbp     pressure at the onset (end-diastolic pressure)
    map     mean of the wave over the beat (exact, not the dbp + (sbp - dbp) / 3 approximation)

Peaks and onsets are located on a low-pass filtered copy (robust against noise), but the values are read from the
unfiltered wave at its original sampling rate, so peaks are not attenuated by filtering or resampling. Only the checks
needed to tell a beat from a non-beat are applied (beat length, pulse pressure, no gap inside the beat); further
quality filtering is left to the harmonization step.

NaN gaps (e.g. removed PulseDB segments) split the wave into runs; beats never span a gap, and the first and last
(incomplete) beat of every run are dropped.
"""
import numpy as np
import pandas as pd
from scipy.signal import butter, find_peaks, sosfiltfilt

MIN_BEAT_S = 0.25            # 240 bpm
MAX_BEAT_S = 2.0             # 30 bpm
MIN_PULSE_PRESSURE = 5.0     # mmHg; also the minimum prominence of a systolic peak (rejects flat lines, flushes)
PROMINENCE_RATIO = 0.5       # a peak must reach this fraction of the prominence of its neighbours (rejects
PROMINENCE_NEIGHBOURS = 15   # dicrotic waves); local median over this many peaks, so slow BP trends are followed
MIN_RUN_S = 3.0              # shorter runs between gaps are skipped
DETECTION_CUTOFF_HZ = 8.0    # low-pass for locating peaks and onsets only

BEAT_COLUMNS = ["time_s", "sbp", "dbp", "map"]


def _finite_runs(x: np.ndarray) -> list[tuple[int, int]]:
    """(start, stop) of every stretch without NaN."""
    finite = np.concatenate([[False], np.isfinite(x), [False]])
    edges = np.flatnonzero(np.diff(finite.astype(np.int8)))
    return list(zip(edges[::2], edges[1::2]))


def _beats_in_run(seg: np.ndarray, fs: float) -> tuple[np.ndarray, ...]:
    """Onset index, sbp, dbp, map of every complete beat in one gap-free run."""
    empty = (np.array([], dtype=int),) + (np.array([]),) * 3
    smooth = sosfiltfilt(butter(2, DETECTION_CUTOFF_HZ, fs=fs, output="sos"), seg)
    peaks, props = find_peaks(smooth, distance=max(1, int(MIN_BEAT_S * fs)), prominence=MIN_PULSE_PRESSURE)
    if len(peaks) < 3:
        return empty

    prominence = props["prominences"]
    local = pd.Series(prominence).rolling(PROMINENCE_NEIGHBOURS, center=True, min_periods=1).median().to_numpy()
    peaks = peaks[prominence >= PROMINENCE_RATIO * local]
    if len(peaks) < 3:
        return empty

    # onsets[i] = minimum between peaks[i] and peaks[i + 1], i.e. the onset of the beat whose peak is peaks[i + 1]
    onsets = np.array([p0 + np.argmin(smooth[p0:p1]) for p0, p1 in zip(peaks[:-1], peaks[1:])])
    start, end = onsets[:-1], onsets[1:]  # complete beats: from one onset to the next
    sbp = seg[peaks[1:-1]]  # raw value at the peak (the raw maximum of the beat would be biased upwards by noise)
    dbp = seg[start]
    cumsum = np.concatenate([[0.0], np.cumsum(seg, dtype=float)])
    mean = (cumsum[end] - cumsum[start]) / (end - start)

    duration = (end - start) / fs
    ok = (duration >= MIN_BEAT_S) & (duration <= MAX_BEAT_S) & (sbp - dbp >= MIN_PULSE_PRESSURE)
    return start[ok], sbp[ok], dbp[ok], mean[ok]


def detect_beats(wave: np.ndarray, fs: float, t0: float = 0.0) -> pd.DataFrame:
    """One row per beat: time_s (onset, t0 + i / fs), sbp, dbp, map. wave in mmHg; may contain NaN gaps."""
    wave = np.asarray(wave, dtype=float)
    parts = []
    for a, b in _finite_runs(wave):
        if b - a < MIN_RUN_S * fs:
            continue
        onset, sbp, dbp, mean = _beats_in_run(wave[a:b], fs)
        parts.append(pd.DataFrame({"time_s": t0 + (a + onset) / fs, "sbp": sbp, "dbp": dbp, "map": mean}))
    if not parts:
        return pd.DataFrame({c: pd.Series(dtype="float64") for c in BEAT_COLUMNS})
    return pd.concat(parts, ignore_index=True)


def bp_from_wave(record_id: str, wave: np.ndarray, fs: float, t0: float, method: str) -> pd.DataFrame:
    """Rows for the bp table (BP_COLUMNS): one per detected beat, marked source = "calculated"."""
    beats = detect_beats(wave, fs, t0)
    return beats.assign(record_id=record_id, method=method, source="calculated")
