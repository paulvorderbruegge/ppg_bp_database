"""Harmonization: structured datasets -> one training-ready dataset (fixed rate, fixed-length windows).
See README.md ("Harmonization") and pipeline.py."""
from ppg_bp_database.harmonize.config import CONFIGS, HarmonizeConfig
from ppg_bp_database.harmonize.pipeline import harmonize

__all__ = ["CONFIGS", "HarmonizeConfig", "harmonize"]
