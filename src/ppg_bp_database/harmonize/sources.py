"""Which structured signals each dataset contributes to the harmonized data."""
from dataclasses import dataclass

from ppg_bp_database.harmonize.config import HarmonizeConfig


@dataclass(frozen=True)
class Source:
    ppg: str                    # PPG signal name in the structured HDF5 files
    abp: str | None             # continuous pressure wave, or None if the dataset has none
    ppg_inverted: bool = False  # True if the systolic upstroke points down (checked with windows.ppg_dskew)
    periodic: bool = False      # record is one periodic beat (pwdb): tiled for filtering, then cropped back


def sources(config: HarmonizeConfig) -> dict[str, Source]:
    return {
        "ppgbp": Source(ppg="ppg", abp=None),
        # pttppg: slow rise and sharp drop in every channel, derivative skew < 0 in 63 of 66 records -> inverted
        "pttppg": Source(ppg=config.pttppg_ppg, abp=None, ppg_inverted=True),
        "bpgraphene": Source(ppg="ppg", abp="finapres_bp"),  # finapres_ppg is only 75 Hz
        "mimicbp": Source(ppg="ppg", abp="abp"),
        "pulsedb": Source(ppg="ppg", abp="abp"),
        # finger PPG like the measured datasets; brachial pressure because pwdb's SBP/DBP are brachial
        "pwdb": Source(ppg="ppg_digital", abp="abp_brachial", periodic=True),
    }
