"""Round-trip error (Rekonstruktionsfehler) of resampling to the harmonized rate, to justify that rate.

    uv run python scripts/resampling_error.py           # all converted datasets
    uv run python scripts/resampling_error.py ppgbp

Each signal the harmonization uses is resampled native -> 125 Hz -> native (same function as the pipeline) and
compared with the original. The error is what the downsampling lost, in % of the signal's standard deviation:
    error        whole signal
    error_band   only in the PPG band (0.5-15 Hz), i.e. what the model sees
The *_fft columns do the same round trip with Fourier resampling (scipy.signal.resample) for comparison.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.signal import resample as fourier_resample

from ppg_bp_database.access import available_datasets, load_records, wave_file
from ppg_bp_database.harmonize import CONFIGS
from ppg_bp_database.harmonize.pipeline import _read_signal
from ppg_bp_database.harmonize.signals import bandpass, resample
from ppg_bp_database.harmonize.sources import sources

CONFIG = CONFIGS["default"]
RECORDS_PER_DATASET = 200
EDGE_S = 0.5  # left out at both ends: resampling edge effects


def nrmse(a: np.ndarray, b: np.ndarray, fs: float) -> float:
    e = int(EDGE_S * fs)
    a, b = a[e:-e], b[e:-e]
    return 100 * float(np.sqrt(np.mean((a - b) ** 2)) / a.std())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("datasets", nargs="*")
    names = parser.parse_args().datasets or available_datasets()

    rows = []
    for name in names:
        src = sources(CONFIG)[name]
        records = load_records([name])
        records = records.sample(min(len(records), RECORDS_PER_DATASET), random_state=0)
        for _, rec in records.iterrows():
            for signal in [src.ppg, src.abp]:
                if signal is None:
                    continue
                x, fs, _ = _read_signal(wave_file(rec), signal)
                x = x[np.isfinite(x)]  # drop NaN (pwdb's last sample; gaps)
                if src.periodic:
                    x = np.tile(x, 5)  # pwdb: one beat -> continuous signal
                if fs == CONFIG.fs or len(x) < 3 * EDGE_S * fs:
                    continue
                x = x - x.mean()  # only the waveform matters; a large offset would distort the upsampling step
                back = resample(resample(x, fs, CONFIG.fs), CONFIG.fs, fs)[:len(x)]
                back_fft = fourier_resample(fourier_resample(x, round(len(x) * CONFIG.fs / fs)), len(x))
                band = (CONFIG.ppg_band_hz[0], CONFIG.ppg_band_hz[1])
                x_band = bandpass(x, fs, band)
                rows.append({"dataset": name, "signal": signal, "fs_hz": fs,
                             "error": nrmse(x, back, fs),
                             "error_band": nrmse(x_band, bandpass(back, fs, band), fs),
                             "error_fft": nrmse(x, back_fft, fs),
                             "error_band_fft": nrmse(x_band, bandpass(back_fft, fs, band), fs)})

    df = pd.DataFrame(rows)
    print(f"native -> {CONFIG.fs:g} Hz -> native: error in % of the signal's std (median over records)\n")
    print(df.groupby(["dataset", "signal", "fs_hz"])[["error", "error_band", "error_fft", "error_band_fft"]].median().round(2).to_string())


if __name__ == "__main__":
    main()
