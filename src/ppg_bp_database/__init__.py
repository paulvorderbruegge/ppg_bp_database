"""PPG / blood-pressure datasets in one shared format. See README.md."""
from ppg_bp_database.access import (available_datasets, load_bp, load_harmonized_subjects, load_records,
                                    load_subjects, load_wave, load_window_bp, load_window_waves, load_windows,
                                    summary, wave_file)


def main() -> None:
    """`ppg-bp-database` command: overview of the converted datasets."""
    import pandas as pd

    from ppg_bp_database.utils import STRUCTURED_DIR

    table = summary()
    if table.empty:
        print(f"No converted datasets in {STRUCTURED_DIR}. Run: uv run python scripts/structure.py")
        return
    with pd.option_context("display.max_colwidth", 60, "display.width", 200):
        print(table.to_string(index=False))


__all__ = ["available_datasets", "load_bp", "load_harmonized_subjects", "load_records", "load_subjects", "load_wave",
           "load_window_bp", "load_window_waves", "load_windows", "main", "summary", "wave_file"]
