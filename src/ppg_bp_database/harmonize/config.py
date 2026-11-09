"""Parameters of one harmonization run. Every run is saved under data/harmonized/<name>/ together with its config.json,
so differently harmonized versions (e.g. 10 s and 30 s windows) can exist side by side."""
import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class HarmonizeConfig:
    name: str = "default"

    # --- windows
    fs: float = 125.0              # target rate; native rate of pulsedb / mimicbp, PPG content is far below 62.5 Hz
    window_s: float = 10.0         # records shorter than this keep their own length (windows.full_length = False)
    step_s: float = 10.0           # distance between window starts; = window_s -> no overlap
    pad_s: float = 1.0             # extra signal read on both sides, so filter edge effects fall outside the window
    cuff_window_s: float = 60.0    # cuff records: only windows within this distance of a cuff reading are used

    # --- PPG processing (ABP is only resampled: it stays in mmHg)
    ppg_band_hz: tuple[float, float] = (0.5, 15.0)   # zero-phase band-pass, then z-score per window
    pttppg_ppg: str = "ppg_1"      # pttppg has 6 PPG channels (1-3 distal phalanx); wavelengths are not documented

    # --- beat quality (calculated beats outside these limits are dropped)
    sbp_range: tuple[float, float] = (60.0, 250.0)
    dbp_range: tuple[float, float] = (30.0, 150.0)
    min_pulse_pressure: float = 10.0

    # --- window quality
    min_beats_per_s: float = 0.5   # windows with calculated beats need >= this many valid beats per second (30 bpm)
    max_sbp_iqr: float = 15.0      # SBP spread (IQR) of the beats inside one window; larger = artifact / flush
    min_ppg_std: float = 1e-6      # std of the band-passed PPG (before z-scoring); below = flat line
    max_clip_fraction: float = 0.05  # fraction of raw PPG samples sitting at the window's min or max (saturation)
    max_beat_count_diff: int = 2   # PPG peaks vs ABP peaks in a window (only where an ABP wave exists)
    peak_prominence_ratio: float = 0.5  # for that count: peaks must reach this x the typical peak prominence

    def save(self, path: Path) -> None:
        path.write_text(json.dumps(asdict(self), indent=2))

    @classmethod
    def load(cls, path: Path) -> "HarmonizeConfig":
        data: dict[str, Any] = json.loads(path.read_text())
        for f in fields(cls):  # JSON has no tuples: turn lists back into tuples (e.g. ppg_band_hz)
            if isinstance(data.get(f.name), list):
                data[f.name] = tuple(data[f.name])
        return cls(**data)


CONFIGS: dict[str, HarmonizeConfig] = {
    "default": HarmonizeConfig(),
}
