"""MIMIC-III demographics lookup (structure/mimic3_clinical.py) and its use in pulsedb, on synthetic MIMIC-III tables.

Run: uv run pytest
"""
import numpy as np
import pandas as pd
import pytest

from ppg_bp_database.structure import mimic3_clinical as m3
from ppg_bp_database.structure import pulsedb

CHART_COLUMNS = ["ROW_ID", "SUBJECT_ID", "HADM_ID", "ICUSTAY_ID", "ITEMID", "CHARTTIME", "STORETIME", "CGID", "VALUE",
                 "VALUENUM", "VALUEUOM", "WARNING", "ERROR", "RESULTSTATUS", "STOPPED"]


def _chart(subject, hadm, item, time, value, error=None):
    return {"ROW_ID": 0, "SUBJECT_ID": subject, "HADM_ID": hadm, "ICUSTAY_ID": None, "ITEMID": item,
            "CHARTTIME": time, "STORETIME": time, "CGID": 1, "VALUE": str(value), "VALUENUM": value,
            "VALUEUOM": "x", "WARNING": None, "ERROR": error, "RESULTSTATUS": None, "STOPPED": None}


@pytest.fixture
def mimic3(tmp_path, monkeypatch):
    d = tmp_path / "mimic3"
    (d / "mimic3wdb-matched").mkdir(parents=True)
    pd.DataFrame({
        "ROW_ID": [1, 2, 3, 4], "SUBJECT_ID": [1, 2, 3, 4], "GENDER": ["M", "F", "F", "M"],
        # subject 2: older than 89, DOB shifted ~300 years
        "DOB": ["2100-01-01 00:00:00", "1860-01-01 00:00:00", "2100-06-01 00:00:00", "2090-01-01 00:00:00"],
    }).to_csv(d / "PATIENTS.csv.gz", index=False)
    pd.DataFrame({
        "ROW_ID": [1, 2, 3, 4, 5], "SUBJECT_ID": [1, 2, 2, 3, 4], "HADM_ID": [11, 21, 22, 31, 41],
        "ADMITTIME": ["2150-01-01 00:00:00", "2160-03-01 00:00:00", "2165-06-01 00:00:00", "2170-01-01 00:00:00",
                      "2140-01-01 00:00:00"],
        "DISCHTIME": ["2150-01-10 00:00:00", "2160-03-05 00:00:00", "2165-06-10 00:00:00", "2170-01-05 00:00:00",
                      "2140-01-05 00:00:00"],
    }).to_csv(d / "ADMISSIONS.csv.gz", index=False)
    pd.DataFrame([
        _chart(1, 11, 920, "2150-01-01 12:00:00", 70),              # 70 in = 177.8 cm
        _chart(1, 11, 762, "2150-01-02 09:00:00", 80),              # closest weight to the record
        _chart(1, 11, 763, "2150-01-08 09:00:00", 90),
        _chart(1, 11, 224639, "2150-01-02 10:00:00", 60, error=1),  # marked as error
        _chart(1, 11, 763, "2150-01-02 10:00:00", 1000),            # implausible
        _chart(1, 11, 211, "2150-01-02 10:00:00", 85),              # heart rate, not used
        _chart(2, 22, 226730, "2165-06-02 00:00:00", 165),          # height only in the later admission
        _chart(2, 21, 226531, "2160-03-02 00:00:00", 220),          # lbs
        _chart(3, 31, 3485, "2170-01-02 00:00:00", 160),
    ], columns=CHART_COLUMNS).to_csv(d / "CHARTEVENTS.csv.gz", index=False)
    (d / m3.RECORDS_FILE).write_text("\n".join([
        "p00/p000001/p000001-2150-01-02-10-00",
        "p00/p000002/p000002-2165-06-02-08-00",
        "p00/p000002/p000002-2160-03-01-20-00",  # earlier, other admission -> ambiguous, this one used
        "p00/p000003/p000003-2171-01-01-00-00",  # no admission around it
    ]) + "\n")
    monkeypatch.setattr(m3, "MIMIC3_DIR", d)
    monkeypatch.setattr(m3, "CACHE_FILE", tmp_path / "structured" / "mimic3_clinical" / "height_weight.parquet")
    return d


def test_record_start():
    t = m3.record_start(pd.Series(["p000188-2161-12-09-17-50", "2179-10-26-00-15", "c1487", None]))
    assert t.iloc[0] == pd.Timestamp("2161-12-09 17:50") and t.iloc[1] == pd.Timestamp("2179-10-26 00:15")
    assert t.iloc[2:].isna().all()


def test_demographics(mimic3):
    assert m3.available()
    query = pd.DataFrame({"key": ["a", "b", "c", "d"], "SUBJECT_ID": [1, 2, 3, 4]})
    query = query.merge(m3.record_times(), on="SUBJECT_ID", how="left")
    demo = m3.demographics(query)
    assert list(demo.columns) == m3.OUTPUT_COLUMNS
    assert m3.CACHE_FILE.is_file()

    a, b, c, d = (demo.loc[k] for k in "abcd")
    assert a["mimic3_match"] == "unique" and a["mimic3_hadm_id"] == 11
    assert a["age"] == pytest.approx(50.0, abs=0.01) and not a["age_capped"] and a["sex"] == "male"
    assert a["height_cm"] == pytest.approx(177.8) and a["weight_kg"] == 80

    assert b["mimic3_match"] == "ambiguous" and b["mimic3_hadm_id"] == 21
    assert b["age"] == 90 and b["age_capped"] and b["sex"] == "female"
    assert b["height_cm"] == 165  # from the other admission
    assert b["weight_kg"] == pytest.approx(220 * 0.45359237)

    assert c["mimic3_match"] == "none" and pd.isna(c["mimic3_hadm_id"])
    assert np.isnan(c["age"]) and np.isnan(c["weight_kg"]) and c["height_cm"] == 160 and c["sex"] == "female"

    assert d["mimic3_match"] == "none" and np.isnan(d["height_cm"]) and d["sex"] == "male"

    # second call reads the cache, CHARTEVENTS no longer needed
    (mimic3 / "CHARTEVENTS.csv.gz").unlink()
    assert m3.available()
    pd.testing.assert_frame_equal(m3.demographics(query), demo)


def test_pulsedb_height_weight(mimic3):
    subjects = pd.DataFrame({
        "subject_id": ["mimic3_p000001", "pulsedb_vital_p004409"], "dataset": "pulsedb",
        "age": [50.0, 72.0], "sex": ["male", "male"], "height_cm": [np.nan, 157.4], "weight_kg": [np.nan, 59.7],
        "case_id": ["2150-01-02-10-00", "c1487"],
    })
    out = pulsedb.add_mimic3_height_weight(subjects).set_index("subject_id")
    mimic, vital = out.loc["mimic3_p000001"], out.loc["pulsedb_vital_p004409"]
    assert mimic["height_cm"] == pytest.approx(177.8) and mimic["weight_kg"] == 80
    assert mimic["mimic3_match"] == "unique" and mimic["age"] == 50
    assert vital["height_cm"] == 157.4 and pd.isna(vital["mimic3_match"])


def test_missing_mimic3(tmp_path, monkeypatch):
    monkeypatch.setattr(m3, "MIMIC3_DIR", tmp_path / "nothing")
    monkeypatch.setattr(m3, "CACHE_FILE", tmp_path / "nothing.parquet")
    assert not m3.available()
    subjects = pd.DataFrame({"subject_id": ["mimic3_p000001"], "height_cm": [np.nan], "weight_kg": [np.nan],
                             "case_id": ["2150-01-02-10-00"]})
    assert pulsedb.add_mimic3_height_weight(subjects)["height_cm"].isna().all()
