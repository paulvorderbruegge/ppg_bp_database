"""Where the windows of a record lie and which BP rows belong to them. Pure pandas, no wave I/O.

The kind of BP values in the structured bp table decides the windowing, so there is no per-dataset code:

    calculated beats (time_s set, source "calculated")   windows tile the record; a window keeps the valid beats
                                                         whose onset lies inside it (and is rejected if too few)
    timed readings   (time_s set, source "dataset")      e.g. cuff at start / end: windows only within
                                                         cuff_window_s of a reading; each keeps its nearest reading
    whole-record     (time_s NaN)                        windows tile the record; each keeps the value (time_s NaN)

"Tile the record" means: tile each usable segment (segments.py) separately, so no window covers a gap.
Records shorter than one window give a single window of their own length (full_length = False).
If a record has calculated beats, its dataset-given values are not used (e.g. MIMIC-BP's 30 s medians: they do not
describe a 10 s window; they stay in the structured table for validation).
All times here are on the record's time axis; pipeline.py makes the BP times relative to the window start.
"""
from collections import Counter
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ppg_bp_database.harmonize.config import HarmonizeConfig


@dataclass
class RecordPlan:
    starts: np.ndarray                       # window start times (record time axis)
    duration_s: float                        # same for all windows of the record
    full_length: bool                        # False = record shorter than config.window_s
    bp: list[np.ndarray]                     # per window: positions of its rows in the record's bp table
    rejected: Counter                        # reason -> number of windows


def valid_beats(beats: pd.DataFrame, config: HarmonizeConfig) -> pd.Series:
    """Beat-level quality check (physiological limits)."""
    return (beats["sbp"].between(*config.sbp_range) & beats["dbp"].between(*config.dbp_range)
            & (beats["sbp"] - beats["dbp"] >= config.min_pulse_pressure))


def _tile(begin: float, end: float, config: HarmonizeConfig) -> np.ndarray:
    """Window starts in [begin, end) so that each window fits completely."""
    n = int(np.floor((end - begin - config.window_s) / config.step_s + 1e-6)) + 1
    return begin + config.step_s * np.arange(max(n, 0))


def _tile_all(segments: list[tuple[float, float]], config: HarmonizeConfig) -> np.ndarray:
    return np.concatenate([_tile(a, b, config) for a, b in segments])


def plan_record(segments: list[tuple[float, float]], bp: pd.DataFrame, config: HarmonizeConfig) -> RecordPlan:
    """Windows of one record. segments: sorted usable time segments (segments.py, record time axis).
    bp: the record's rows of the structured bp table; the plan refers to them by position (bp.iloc[...])."""
    full_length = segments[-1][1] - segments[0][0] >= config.window_s
    if full_length:
        duration = config.window_s
    else:  # short record: one window over its longest segment
        t_begin, t_end = max(segments, key=lambda s: s[1] - s[0])
        duration = t_end - t_begin
    rejected: Counter = Counter()

    bp = bp.reset_index(drop=True)
    timed = bp["time_s"].notna()
    beats = bp[timed & (bp["source"] == "calculated")]
    readings = bp[timed & (bp["source"] != "calculated")]
    whole = bp[~timed]

    if len(beats):
        starts = _tile_all(segments, config) if full_length else np.array([t_begin])
        beats = beats[valid_beats(beats, config)].sort_values("time_s")
        times, pos = beats["time_s"].to_numpy(), beats.index.to_numpy()
        lo, hi = np.searchsorted(times, starts, "left"), np.searchsorted(times, starts + duration, "left")
        min_beats = config.min_beats_per_s * duration
        keep, rows = [], []
        for s, a, b in zip(starts, lo, hi):
            if b - a < min_beats:
                rejected["too_few_beats"] += 1
            else:
                keep.append(s)
                rows.append(pos[a:b])
        return RecordPlan(np.array(keep), duration, full_length, rows, rejected)

    if len(readings):
        if not full_length:
            starts = np.array([t_begin])
        else:
            starts = np.unique(np.concatenate([
                _tile(max(a, t - config.cuff_window_s), min(b, t + config.cuff_window_s), config)
                for t in readings["time_s"] for a, b in segments]))
        centers = starts + duration / 2
        nearest = np.abs(centers[:, None] - readings["time_s"].to_numpy()[None, :]).argmin(axis=1)
        pos = readings.index.to_numpy()
        return RecordPlan(starts, duration, full_length, [pos[[i]] for i in nearest], rejected)

    if len(whole):
        starts = _tile_all(segments, config) if full_length else np.array([t_begin])
        return RecordPlan(starts, duration, full_length, [whole.index.to_numpy()] * len(starts), rejected)

    return RecordPlan(np.array([]), duration, full_length, [], Counter({"record_without_bp": 1}))
