from pathlib import Path
from ppg_bp_database.datasets import ptt_ppg_utils


if __name__ == "__main__":
    from pathlib import Path
    ptt_ppg_utils.convert("data/raw/pttppg", "data/structured/pttppg")