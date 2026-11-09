"""Step 1, structuring: raw data (data/raw/<dataset>) -> shared format (data/structured/<dataset>).

beats.py calculates beat-to-beat BP from continuous pressure waves. Besides that, one module per dataset; every
dataset module has the same layout:

    DATASET                          dataset ID; also the folder name below data/raw and data/structured
    SYNTHETIC                        True if the data is simulated (e.g. pwdb), False if measured
    build_subjects(raw_dir)          -> subjects table
    write_records(raw_dir, out_dir)  -> records table (and writes one HDF5 wave file per record)
    build_bp(raw_dir, records)       -> bp table, or None if the dataset has no BP values yet
    convert(raw_dir, out_dir)        runs all of the above and saves the tables via utils.save_tables

(pulsedb builds subjects and records in one pass, see its docstring.)
"""
from types import ModuleType

from ppg_bp_database.structure import bpgraphene, mimicbp, ppgbp, pttppg, pulsedb, pwdb

DATASETS: dict[str, ModuleType] = {m.DATASET: m for m in (ppgbp, pttppg, bpgraphene, mimicbp, pulsedb, pwdb)}
