"""Step 1, structuring: raw datasets (data/raw/<dataset>) -> shared format (data/structured/<dataset>).

    uv run python scripts/structure.py                 # all datasets whose raw folder is not empty
    uv run python scripts/structure.py ppgbp pttppg    # only these
"""
import argparse
import sys

from ppg_bp_database.structure import DATASETS
from ppg_bp_database.utils import RAW_DIR, STRUCTURED_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("datasets", nargs="*", choices=[[], *DATASETS], metavar="dataset",
                        help=f"any of: {', '.join(DATASETS)} (default: all with raw data)")
    names = parser.parse_args().datasets

    if not names:
        names = [n for n in DATASETS if (RAW_DIR / n).is_dir() and any((RAW_DIR / n).iterdir())]
        print(f"Datasets with raw data: {', '.join(names) or 'none'}")

    failed = []
    for name in names:
        print(f"--- {name}")
        try:
            DATASETS[name].convert(RAW_DIR / name, STRUCTURED_DIR / name)
        except Exception as e:
            if len(names) == 1:
                raise
            print(f"FAILED: {type(e).__name__}: {e}")
            failed.append(name)
    if failed:
        sys.exit(f"Failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
