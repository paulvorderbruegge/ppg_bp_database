"""MIMIC-III Clinical Database v1.4: demographics for the MIMIC-based datasets (mimicbp, the MIMIC part of pulsedb).

Not a dataset of its own (not in DATASETS); the dataset modules call `demographics`.
Raw: data/raw/mimic3 (scripts/download_mimic3.sh): PATIENTS, ADMISSIONS, CHARTEVENTS (+ D_ITEMS) as .csv.gz and
     mimic3wdb-matched/RECORDS-waveforms (open access).

Which admission? A subject can have several hospital admissions, years apart. The waveforms of both datasets come
from the MIMIC-III Waveform Database Matched Subset, whose record names carry the record's start time on the same
(date-shifted) time axis as the clinical database, e.g. p000188-2161-12-09-17-50. Every record start is assigned to
the subject's admission whose [ADMITTIME, DISCHTIME] contains it (or lies within ADMISSION_TOLERANCE of it; the
nearest one if several).
    pulsedb:  the record is known (field CaseID = start of the record) -> one record per subject.
    mimicbp:  only the subject is known -> all of the subject's waveform records (RECORDS-waveforms).
If the records lie in more than one admission, the admission of the earliest record is used.

Subject IDs. Both datasets contain patients of the same Matched Subset, partly the same people. Their subjects get the
dataset-independent ID "mimic3_<MIMIC id>" (e.g. mimic3_p000188, see `subject_id`), so a person in both datasets has
the same subject_id and splits by subject_id keep them on one side. Such a person has one subjects row per dataset
(values may differ: other admission, PulseDB's own age / sex): the subjects table is keyed by (subject_id, dataset).
    mimic3_match = "unique" (one admission), "ambiguous" (several; earliest used) or "none" (no record within an
                   admission: age and weight stay NaN, sex and height are still filled)
    mimic3_hadm_id = HADM_ID of the admission used

Values. Reference time = start of the earliest waveform record in the admission used.
    age          (reference time - DOB) in years. MIMIC-III shifts the DOB of patients older than 89 so that they
                 appear ~300 years old; their age is set to AGE_CAP (90, the lower bound) and age_capped = True.
    sex          PATIENTS.GENDER
    height_cm    median of the height measurements of the admission; if there are none, of all the subject's
                 admissions (adult height hardly changes)
    weight_kg    the measurement closest to the reference time within the admission; NaN if there is none
Height and weight are read from CHARTEVENTS (item IDs as in the MIMIC code repository, concepts heightfirstday /
weightfirstday, without their newborn items 3485/3486/4187/4188 lengths and 3580/3581 Present Weight: in MIMIC-III
these are charted almost only for newborns, ~2 kg / ~44 cm, so none survive the adult range) and filtered to plausible adult ranges. CHARTEVENTS is 4 GB gzipped, so
the measurements are extracted once into data/structured/mimic3_clinical/height_weight.parquet (delete it to extract
again; afterwards CHARTEVENTS.csv.gz itself is no longer needed).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pcsv

from ppg_bp_database.utils import RAW_DIR, STRUCTURED_DIR, normalize_sex

MIMIC3_DIR = RAW_DIR / "mimic3"
CACHE_FILE = STRUCTURED_DIR / "mimic3_clinical" / "height_weight.parquet"
RECORDS_FILE = "mimic3wdb-matched/RECORDS-waveforms"

ADMISSION_TOLERANCE = pd.Timedelta(hours=24)  # a record may start shortly before admission / end after discharge
AGE_CAP = 90.0
SHIFTED_AGE = 120.0  # computed ages above this are the shifted ones of patients older than 89
HEIGHT_RANGE_CM = (120.0, 230.0)
WEIGHT_RANGE_KG = (25.0, 350.0)

INCH_CM, POUND_KG = 2.54, 0.45359237
ITEMS = {  # CHARTEVENTS ITEMID -> (quantity, factor to cm / kg). CareVue (< 220000) and MetaVision (>= 220000)
    920: ("height", INCH_CM),     # Admit Ht
    1394: ("height", INCH_CM),    # Height Inches
    226707: ("height", INCH_CM),  # Height
    226730: ("height", 1.0),      # Height (cm)
    762: ("weight", 1.0),         # Admit Wt
    763: ("weight", 1.0),         # Daily Weight
    224639: ("weight", 1.0),      # Daily Weight
    226512: ("weight", 1.0),      # Admission Weight (Kg)
    226531: ("weight", POUND_KG), # Admission Weight (lbs.)
}

OUTPUT_COLUMNS = ["age", "age_capped", "sex", "height_cm", "weight_kg", "mimic3_hadm_id", "mimic3_match"]


def _csv(mimic3_dir: Path, table: str) -> Path:
    return mimic3_dir / f"{table}.csv.gz"


def available(mimic3_dir: Path | None = None, cache_file: Path | None = None) -> bool:
    """True if all files needed by `demographics` are there (CHARTEVENTS only until its extract exists)."""
    mimic3_dir, cache_file = mimic3_dir or MIMIC3_DIR, cache_file or CACHE_FILE
    needed = [_csv(mimic3_dir, "PATIENTS"), _csv(mimic3_dir, "ADMISSIONS"), mimic3_dir / RECORDS_FILE]
    return all(p.is_file() for p in needed) and (cache_file.is_file() or _csv(mimic3_dir, "CHARTEVENTS").is_file())


def missing_message(fields: str = "age / sex / height / weight") -> str:
    return f"MIMIC-III clinical data not found in {MIMIC3_DIR} (scripts/download_mimic3.sh): {fields} stay empty"


SUBJECT_PREFIX = "mimic3_"


def subject_id(source_ids: pd.Series | str) -> pd.Series | str:
    """'p000188' -> 'mimic3_p000188': shared subject_id of a MIMIC-III patient in every MIMIC-based dataset."""
    return SUBJECT_PREFIX + source_ids


def subject_number(ids: pd.Series) -> pd.Series:
    """'p000188' or 'mimic3_p000188' -> 188 (= MIMIC-III SUBJECT_ID)."""
    numbers = ids.str.extract(rf"^(?:{SUBJECT_PREFIX})?p(\d+)$", expand=False)
    assert numbers.notna().all(), f"not a MIMIC subject ID: {ids[numbers.isna()].iloc[0]!r}"
    return numbers.astype("int64")


def record_start(names: pd.Series) -> pd.Series:
    """Start time from a record name / PulseDB CaseID: 'p000188-2161-12-09-17-50' or '2161-12-09-17-50'.
    Unparseable values become NaT."""
    stamp = names.astype("string").str.extract(r"(\d{4}-\d{2}-\d{2}-\d{2}-\d{2})$", expand=False)
    return pd.to_datetime(stamp, format="%Y-%m-%d-%H-%M", errors="coerce").astype("datetime64[s]")


def record_times(mimic3_dir: Path | None = None) -> pd.DataFrame:
    """SUBJECT_ID, time: start of every record of the Waveform Database Matched Subset."""
    mimic3_dir = mimic3_dir or MIMIC3_DIR
    paths = pd.read_csv(mimic3_dir / RECORDS_FILE, header=None, names=["path"], dtype="string")["path"]
    names = paths.str.rsplit("/", n=1).str[-1]
    out = pd.DataFrame({"SUBJECT_ID": subject_number(names.str.split("-").str[0]), "time": record_start(names)})
    assert out["time"].notna().all(), f"unexpected record names in {RECORDS_FILE}"
    return out


def _read_dates(path: Path, columns: list[str], dates: list[str]) -> pd.DataFrame:
    """Read a small table; dates in seconds resolution (nanoseconds overflow for the shifted ~300-year ages)."""
    df = pd.read_csv(path, usecols=columns)
    for c in dates:
        df[c] = pd.to_datetime(df[c]).astype("datetime64[s]")
    return df


def extract_height_weight(mimic3_dir: Path = MIMIC3_DIR, out_file: Path = CACHE_FILE) -> pd.DataFrame:
    """Stream CHARTEVENTS once, keep the height / weight items, convert to cm / kg, drop implausible values and
    save SUBJECT_ID, HADM_ID, CHARTTIME, ITEMID, quantity, value to out_file."""
    path = _csv(mimic3_dir, "CHARTEVENTS")
    print(f"  extracting height / weight from {path.name} (one time only, takes several minutes)", flush=True)
    types = {"SUBJECT_ID": pa.int64(), "HADM_ID": pa.int64(), "ITEMID": pa.int64(), "CHARTTIME": pa.timestamp("s"),
             "VALUENUM": pa.float64(), "ERROR": pa.float64()}
    item_ids = pa.array(list(ITEMS), pa.int64())
    parts, n_rows = [], 0
    reader = pcsv.open_csv(path, read_options=pcsv.ReadOptions(block_size=64 << 20),
                           parse_options=pcsv.ParseOptions(newlines_in_values=True),
                           convert_options=pcsv.ConvertOptions(include_columns=list(types), column_types=types))
    for batch in reader:
        parts.append(batch.filter(pc.is_in(batch["ITEMID"], value_set=item_ids)))
        n_rows += batch.num_rows
        if n_rows // 50_000_000 != (n_rows - batch.num_rows) // 50_000_000:
            print(f"  {n_rows / 1e6:.0f} M rows read", flush=True)
    df = pa.Table.from_batches(parts, schema=reader.schema).to_pandas()

    df = df[(df["ERROR"].fillna(0) != 1) & (df["VALUENUM"] > 0)]  # ERROR = 1: marked as wrong in MetaVision
    df["quantity"] = df["ITEMID"].map(lambda i: ITEMS[i][0]).astype("string")
    df["value"] = df["VALUENUM"] * df["ITEMID"].map(lambda i: ITEMS[i][1])
    lo = df["quantity"].map({"height": HEIGHT_RANGE_CM[0], "weight": WEIGHT_RANGE_KG[0]})
    hi = df["quantity"].map({"height": HEIGHT_RANGE_CM[1], "weight": WEIGHT_RANGE_KG[1]})
    kept = df[df["value"].between(lo, hi)]
    _print_items(mimic3_dir, df, kept, n_rows)

    out = kept[["SUBJECT_ID", "HADM_ID", "CHARTTIME", "ITEMID", "quantity", "value"]].reset_index(drop=True)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_file, index=False)
    return out


def _print_items(mimic3_dir: Path, values: pd.DataFrame, kept: pd.DataFrame, n_rows: int) -> None:
    """Values per item before / after the range filter, with the item labels from D_ITEMS (if downloaded)."""
    counts = pd.DataFrame({"values": values["ITEMID"].value_counts(), "in range": kept["ITEMID"].value_counts()})
    counts = counts.fillna(0).astype(int)
    d_items = _csv(mimic3_dir, "D_ITEMS")
    if d_items.is_file():
        labels = pd.read_csv(d_items, usecols=["ITEMID", "LABEL", "UNITNAME"]).set_index("ITEMID")
        counts = counts.join(labels)
    print(f"  {n_rows:,} CHARTEVENTS rows, height / weight items:\n{counts.sort_index().to_string()}")


def _height_weight(mimic3_dir: Path, cache_file: Path) -> pd.DataFrame:
    if cache_file.is_file():
        return pd.read_parquet(cache_file)
    return extract_height_weight(mimic3_dir, cache_file)


def _assign_admissions(query: pd.DataFrame, admissions: pd.DataFrame) -> pd.DataFrame:
    """query (key, SUBJECT_ID, time) + HADM_ID of the admission containing time (NaN if none within tolerance)."""
    q = query.reset_index(drop=True).rename_axis("row").reset_index()
    m = q.dropna(subset=["time"]).merge(admissions, on="SUBJECT_ID")
    before = (m["ADMITTIME"] - m["time"]).clip(lower=pd.Timedelta(0))  # > 0 if the record starts before admission
    after = (m["time"] - m["DISCHTIME"]).clip(lower=pd.Timedelta(0))   # > 0 if it starts after discharge
    m["distance"] = before + after
    m = m[m["distance"] <= ADMISSION_TOLERANCE].sort_values(["row", "distance"]).drop_duplicates("row")
    return q.merge(m[["row", "HADM_ID"]], on="row", how="left").drop(columns="row")


def demographics(query: pd.DataFrame, mimic3_dir: Path | None = None, cache_file: Path | None = None) -> pd.DataFrame:
    """query: one row per waveform record: key (our subject_id), SUBJECT_ID (MIMIC-III), time (record start, NaT if
    unknown); subjects without a known record need one row with time NaT.
    Returns one row per key, indexed by key, with OUTPUT_COLUMNS (see the module docstring).
    mimic3_dir / cache_file default to MIMIC3_DIR / CACHE_FILE."""
    mimic3_dir, cache_file = mimic3_dir or MIMIC3_DIR, cache_file or CACHE_FILE
    patients = _read_dates(_csv(mimic3_dir, "PATIENTS"), ["SUBJECT_ID", "GENDER", "DOB"], ["DOB"])
    admissions = _read_dates(_csv(mimic3_dir, "ADMISSIONS"), ["SUBJECT_ID", "HADM_ID", "ADMITTIME", "DISCHTIME"],
                             ["ADMITTIME", "DISCHTIME"])
    query = query.assign(time=pd.to_datetime(query["time"]).astype("datetime64[s]"))
    rows = _assign_admissions(query, admissions)

    # per key: admission of the earliest matched record, and how many admissions were matched
    keys = query.groupby("key", sort=False)["SUBJECT_ID"].first().to_frame()
    matched = rows.dropna(subset=["HADM_ID"]).sort_values("time")
    first = matched.drop_duplicates("key").set_index("key")
    out = keys.join(first[["time", "HADM_ID"]].rename(columns={"time": "ref_time"}))
    n_adm = matched.groupby("key")["HADM_ID"].nunique().reindex(out.index, fill_value=0)
    out["mimic3_match"] = pd.Series(np.select([n_adm == 1, n_adm > 1], ["unique", "ambiguous"], "none"),
                                    index=out.index, dtype="string")
    out["mimic3_hadm_id"] = out.pop("HADM_ID").astype("Int64")

    # age and sex
    out = out.join(patients.set_index("SUBJECT_ID"), on="SUBJECT_ID")
    age = (out["ref_time"] - out["DOB"]) / pd.Timedelta(days=365.2425)
    out["age_capped"] = age > SHIFTED_AGE
    out["age"] = age.mask(out["age_capped"], AGE_CAP).astype("float64")
    out["sex"] = normalize_sex(out["GENDER"]).where(out["GENDER"].notna())

    # height: median in the admission, else over all of the subject's admissions
    hw = _height_weight(mimic3_dir, cache_file)
    heights = hw[hw["quantity"] == "height"]
    by_adm = heights.groupby("HADM_ID")["value"].median()
    by_subject = heights.groupby("SUBJECT_ID")["value"].median()
    out["height_cm"] = (out["mimic3_hadm_id"].map(by_adm).astype("float64")
                        .fillna(out["SUBJECT_ID"].map(by_subject)).astype("float64"))

    # weight: closest to the reference time within the admission
    weights = hw[hw["quantity"] == "weight"][["HADM_ID", "CHARTTIME", "value"]]
    w = (out.reset_index()[["key", "mimic3_hadm_id", "ref_time"]].dropna()
         .astype({"mimic3_hadm_id": "int64"})
         .merge(weights, left_on="mimic3_hadm_id", right_on="HADM_ID"))
    w["gap"] = (w["CHARTTIME"].astype("datetime64[s]") - w["ref_time"]).abs()
    closest = w.sort_values("gap").drop_duplicates("key").set_index("key")["value"]
    out["weight_kg"] = closest.reindex(out.index).astype("float64")

    return out[OUTPUT_COLUMNS]
