from pathlib import Path
from ppg_bp_database.datasets import pulsedb_utils


if __name__ == "__main__":
    pulsedb_utils.convert(Path("data/raw/pulsedb"), Path("data/structured/pulsedb"))
