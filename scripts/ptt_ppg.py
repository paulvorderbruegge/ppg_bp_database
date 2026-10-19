from ppg_bp_database.datasets import ptt_ppg_utils

RAW = "data/raw/pttppg"
OUT = "data/structured/pttppg"

if __name__ == "__main__":
    from pathlib import Path
    ptt_ppg_utils.convert(Path(RAW), Path(OUT))