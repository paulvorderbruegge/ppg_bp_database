from pathlib import Path
from ppg_bp_database.datasets import mimic_bp_utils


if __name__ == "__main__":
    mimic_bp_utils.convert(Path("data/raw/mimicbp"), Path("data/structured/mimicbp"))