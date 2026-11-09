"""Harmonize converted datasets (data/structured/<dataset>) into data/harmonized/<config>/.

    uv run python scripts/harmonize.py                       # all converted datasets, config "default"
    uv run python scripts/harmonize.py ppgbp pttppg          # only these (the output then contains only these)
    uv run python scripts/harmonize.py --config default

Configs are defined in src/ppg_bp_database/harmonize/config.py (CONFIGS).
"""
import argparse

from ppg_bp_database.access import available_datasets
from ppg_bp_database.harmonize import CONFIGS, harmonize


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("datasets", nargs="*", metavar="dataset",
                        help=f"converted datasets (default: all, currently {', '.join(available_datasets()) or 'none'})")
    parser.add_argument("--config", default="default", choices=list(CONFIGS))
    args = parser.parse_args()
    harmonize(args.datasets or None, CONFIGS[args.config])


if __name__ == "__main__":
    main()
