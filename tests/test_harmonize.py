"""End-to-end test of the harmonization on synthetic structured datasets (one per harmonization source).

Run: uv run pytest
"""
import numpy as np
import pandas as pd
import pytest

from ppg_bp_database import access
from ppg_bp_database.structure.beats import bp_from_wave
from ppg_bp_database.harmonize import HarmonizeConfig, harmonize
from ppg_bp_database.utils import mean_arterial_pressure, save_tables, write_wave
from synth import abp, ppg_from_abp

CONFIG = HarmonizeConfig(name="test")
FS = CONFIG.fs


def _save(root, dataset, records: list[dict], bp: pd.DataFrame, synthetic=False):
    """records: dicts with record_id, subject_id and signals {name: (data, fs, unit, t0)}."""
    out = root / "structured" / dataset
    rows = []
    for r in records:
        path = f"waves/{r['record_id']}.h5"
        write_wave(out / path, r["signals"])
        rows.append({"record_id": r["record_id"], "subject_id": r["subject_id"], "dataset": dataset,
                     "condition": "test", "duration_s": max(len(d) / fs for d, fs, _, _ in r["signals"].values()),
                     "signals": list(r["signals"]), "sampling_rates_hz": [s[1] for s in r["signals"].values()],
                     "wave_path": path})
    subjects = pd.DataFrame({"subject_id": sorted({r["subject_id"] for r in records}), "dataset": dataset})
    save_tables(out, subjects, pd.DataFrame(rows), bp, synthetic=synthetic)


def _given(record_id, time_s, sbp, dbp, method):
    return pd.DataFrame({"record_id": record_id, "time_s": time_s, "sbp": sbp, "dbp": dbp,
                         "map": mean_arterial_pressure(pd.Series(sbp), pd.Series(dbp)), "method": method,
                         "source": "dataset"})


@pytest.fixture(scope="module")
def harmonized(tmp_path_factory):
    root = tmp_path_factory.mktemp("data")
    rng = np.random.default_rng(0)
    truth = {}

    # mimicbp-like: 30 s at 125 Hz, dataset median (time NaN) + calculated beats
    recs, bps = [], []
    for i in range(3):
        rid = f"mimicbp_p{i}_01"
        x, truth[rid] = abp(125.0, 30, rng)
        recs.append({"record_id": rid, "subject_id": f"mimicbp_p{i}",
                     "signals": {"ppg": (ppg_from_abp(x, 125.0, rng), 125.0, "a.u.", 0.0), "abp": (x, 125.0, "mmHg", 0.0)}})
        bps += [_given(rid, [np.nan], [120.0], [80.0], "abp"), bp_from_wave(rid, x, 125.0, 0.0, "abp")]
    x, _ = abp(125.0, 30, rng)  # good ABP, but the PPG is only noise -> every window must be rejected
    recs.append({"record_id": "mimicbp_noise_01", "subject_id": "mimicbp_noise",
                 "signals": {"ppg": (rng.normal(size=len(x)), 125.0, "a.u.", 0.0), "abp": (x, 125.0, "mmHg", 0.0)}})
    bps.append(bp_from_wave("mimicbp_noise_01", x, 125.0, 0.0, "abp"))
    _save(root, "mimicbp", recs, pd.concat(bps))

    # pulsedb-like: 60 s clip starting at t0 = 100 s with a removed segment (NaN) at 120-130 s
    x, truth["pulsedb_c01"] = abp(125.0, 60, rng)
    truth["pulsedb_c01"][:, 0] += 100
    p = ppg_from_abp(x, 125.0, rng)
    x[20 * 125:30 * 125] = np.nan
    p[20 * 125:30 * 125] = np.nan
    _save(root, "pulsedb", [{"record_id": "pulsedb_c01", "subject_id": "pulsedb_s1",
                             "signals": {"ppg": (p, 125.0, "a.u.", 100.0), "abp": (x, 125.0, "mmHg", 100.0)}}],
          bp_from_wave("pulsedb_c01", x, 125.0, 100.0, "abp"))

    # bpgraphene-like: finapres 200 Hz, ppg 1250 Hz (same underlying pulse), t0 = 3.5 s
    x, truth["bpgraphene_t1"] = abp(200.0, 40, rng)
    truth["bpgraphene_t1"][:, 0] += 3.5
    t200 = np.arange(len(x)) / 200
    t1250 = np.arange(int(40 * 1250)) / 1250
    p = ppg_from_abp(np.interp(t1250, t200, x), 1250.0, rng)
    _save(root, "bpgraphene", [{"record_id": "bpgraphene_t1", "subject_id": "bpgraphene_s1",
                                "signals": {"ppg": (p, 1250.0, "a.u.", 3.5), "finapres_bp": (x, 200.0, "mmHg", 3.5)}}],
          bp_from_wave("bpgraphene_t1", x, 200.0, 3.5, "finapres"))

    # pttppg-like: 200 s, 500 Hz, inverted PPG, cuff at start and end; a second channel that must be ignored
    x, _ = abp(500.0, 200, rng)
    p = -ppg_from_abp(x, 500.0, rng)
    _save(root, "pttppg", [{"record_id": "pttppg_s1_sit", "subject_id": "pttppg_s1",
                            "signals": {"ppg_1": (p, 500.0, "a.u.", 0.0), "ppg_2": (rng.normal(size=len(p)), 500.0, "a.u.", 0.0)}}],
          _given("pttppg_s1_sit", [0.0, 199.0], [118.0, 125.0], [76.0, 79.0], "cuff"))

    # ppgbp-like: 2.1 s at 1000 Hz, cuff at t = 0
    x, _ = abp(1000.0, 2.1, rng)
    _save(root, "ppgbp", [{"record_id": "ppgbp_1_1", "subject_id": "ppgbp_1",
                           "signals": {"ppg": (ppg_from_abp(x, 1000.0, rng), 1000.0, "a.u.", 0.0)}}],
          _given("ppgbp_1_1", [0.0], [121.0], [77.0], "cuff"))

    # pwdb-like: one periodic beat at 500 Hz, last PPG sample invalid (NaN), whole-record value
    from synth import beat_shape
    beat = 75 + 45 * beat_shape(400)
    ppg = 1 + beat / 100
    ppg[-1] = np.nan
    _save(root, "pwdb", [{"record_id": "pwdb_0001", "subject_id": "pwdb_0001",
                          "signals": {"ppg_digital": (ppg, 500.0, "a.u.", 0.0), "abp_brachial": (beat, 500.0, "mmHg", 0.0),
                                      "ppg_radial": (ppg, 500.0, "a.u.", 0.0)}}],
          _given("pwdb_0001", [np.nan], [120.0], [75.0], "abp"), synthetic=True)

    harmonize(None, CONFIG, structured_dir=root / "structured", harmonized_dir=root / "harmonized")
    h = root / "harmonized"
    return {
        "dir": h / "test",
        "windows": access.load_windows("test", harmonized_dir=h),
        "bp": access.load_window_bp("test", harmonized_dir=h),
        "waves": lambda w: access.load_window_waves(w, "test", harmonized_dir=h),
        "truth": truth,
    }


def _windows(h, dataset):
    return h["windows"][h["windows"]["dataset"] == dataset].reset_index(drop=True)


def test_shapes_and_rate(harmonized):
    w = harmonized["windows"]
    arrays = harmonized["waves"](w)
    assert arrays["ppg"].shape == arrays["abp"].shape == (len(w), int(CONFIG.window_s * FS))
    full = w["full_length"].to_numpy()
    assert np.isfinite(arrays["ppg"][full]).all()
    assert (w.loc[full, "n_samples"] == 1250).all()
    assert set(w["dataset"]) == {"mimicbp", "pulsedb", "bpgraphene", "pttppg", "ppgbp", "pwdb"}


def test_beats_match_truth_and_wave(harmonized):
    """Each window's BP rows are the true beats inside it, and time_s points to the onset in the ABP array."""
    for dataset, rid in [("mimicbp", "mimicbp_p0_01"), ("bpgraphene", "bpgraphene_t1"), ("pulsedb", "pulsedb_c01")]:
        w = _windows(harmonized, dataset)
        w = w[w["record_id"] == rid]
        assert len(w)
        arrays = harmonized["waves"](w)
        truth = harmonized["truth"][rid]
        for i, win in enumerate(w.itertuples()):
            rows = harmonized["bp"][harmonized["bp"]["window_id"] == win.window_id]
            assert len(rows) >= 0.5 * CONFIG.window_s
            assert rows["source"].eq("calculated").all()          # dataset medians are not copied
            assert rows["time_s"].between(0, CONFIG.window_s).all()
            onsets = win.t_start_s + rows["time_s"].to_numpy()
            j = np.abs(truth[:, 0][None] - onsets[:, None]).argmin(1)
            assert np.abs(rows["sbp"].to_numpy() - truth[j, 1]).max() < 2.5
            assert np.abs(rows["dbp"].to_numpy() - truth[j, 2]).max() < 2.5
            idx = np.round(rows["time_s"].to_numpy() * FS).astype(int)
            idx = idx[idx < win.n_samples]
            assert np.abs(arrays["abp"][i, idx] - rows["dbp"].to_numpy()[: len(idx)]).max() < 4  # resampled onset


def test_gap_rejected_and_reported(harmonized):
    w = _windows(harmonized, "pulsedb")
    assert w["t_start_s"].round(3).tolist() == [100.0, 110.0, 130.0, 140.0, 150.0]   # 120-130 s is the gap
    report = pd.read_csv(harmonized["dir"] / "report.csv")
    rejected = report[(report["dataset"] == "pulsedb") & (report["reason"] != "kept")]["count"].sum()
    assert rejected == 1


def test_noise_ppg_rejected(harmonized):
    assert not (harmonized["windows"]["record_id"] == "mimicbp_noise_01").any()
    report = pd.read_csv(harmonized["dir"] / "report.csv").set_index(["dataset", "reason"])["count"]
    assert report.get(("mimicbp", "hr_mismatch"), 0) == 3


def test_cuff_windows_near_readings(harmonized):
    w = _windows(harmonized, "pttppg")
    starts = w["t_start_s"].to_numpy()
    near_start, near_end = starts < 100, starts >= 100
    assert (starts[near_start] + CONFIG.window_s <= 60 + 1e-6).all()
    assert (starts[near_end] >= 199 - 60 - 1e-6).all()
    assert len(w) == 6 + 6          # 0..50 and 139..189 (record ends at 200)
    bp = harmonized["bp"].merge(w[["window_id", "t_start_s"]], on="window_id")
    assert (bp.loc[bp["t_start_s"] < 100, "sbp"] == 118).all() and (bp.loc[bp["t_start_s"] >= 100, "sbp"] == 125).all()
    assert (w["ppg_signal"] == "ppg_1").all()
    assert (w["ppg_dskew"] > 0).mean() > 0.9   # inverted PPG was flipped
    assert not w["has_abp"].any()


def test_short_records_keep_their_length(harmonized):
    w = _windows(harmonized, "ppgbp")
    assert len(w) == 1 and not w.loc[0, "full_length"] and w.loc[0, "n_samples"] == 262
    ppg = harmonized["waves"](w)["ppg"][0]
    assert np.isfinite(ppg[:262]).all() and np.isnan(ppg[262:]).all()
    bp = harmonized["bp"][harmonized["bp"]["window_id"] == w.loc[0, "window_id"]]
    assert bp["time_s"].tolist() == [0.0] and bp["sbp"].tolist() == [121.0]


def test_periodic_beat(harmonized):
    w = _windows(harmonized, "pwdb")
    assert len(w) == 1 and w.loc[0, "has_abp"] and w.loc[0, "n_samples"] == 99   # 399 valid samples at 500 Hz
    arrays = harmonized["waves"](w)
    a = arrays["abp"][0, :99]
    assert abs(a.max() - 120) < 3 and abs(a.min() - 75) < 3                    # pressure kept in mmHg
    bp = harmonized["bp"][harmonized["bp"]["window_id"] == w.loc[0, "window_id"]]
    assert bp["time_s"].isna().all() and bp["sbp"].tolist() == [120.0]


def test_outputs(harmonized):
    d = harmonized["dir"]
    assert HarmonizeConfig.load(d / "config.json") == CONFIG
    subjects = pd.read_parquet(d / "subjects.parquet")
    assert set(subjects["subject_id"]) == set(harmonized["windows"]["subject_id"])
    assert harmonized["windows"]["window_id"].is_unique
    assert set(harmonized["bp"]["window_id"]) <= set(harmonized["windows"]["window_id"])
