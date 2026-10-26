from pathlib import Path
from ppg_bp_database.datasets import ppg_bp_utils


if __name__ == "__main__":
    ppg_bp_utils.convert(Path("data/raw/ppgbp"), Path("data/structured/ppgbp"))
