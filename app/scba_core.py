"""Core logic of the SCBA resilient-modulus teaching/planning app (no Streamlit imports -> testable).

Models are fitted at start-up from data/raw/scba-resilient-modulus-rf.xlsx (same screening rules and
hyper-parameters as the notebook), so no pickle/version issues.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from sklearn.ensemble import RandomForestRegressor

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "raw" / "scba-resilient-modulus-rf.xlsx"
SEED = 42

SCBA_LEVELS = (0, 5, 10)
TARGET = "Axial Resilient Modulus-MPa"
S3, SD, SC = "Confining Stress-kPa", "Cyclic Axial Stress-kPa", "Contact Stress-kPa"
RF_FEATURES = ["% CBCA", S3, SD, SC, "theta"]
# Tuned in the notebook (section 7), feature set M4 (inputs known before the test)
RF_PARAMS = dict(n_estimators=50, max_depth=10, min_samples_split=5, max_features=None, random_state=SEED)

# Typical error under leave-one-sequence-out validation (reports/v2_primary_results.csv)
RMSE_RF_V2 = 38.9    # MPa
RMSE_KTH_V2 = 32.5   # MPa

# 15-sequence schedule of the study (matches AASHTO T 307 for base/subbase; verified against the data)
_SIGMA3 = [20.7] * 3 + [34.5] * 3 + [68.9] * 3 + [103.4] * 3 + [137.9] * 3
_SIGMAX = [20.7, 41.4, 62.1, 34.5, 68.9, 103.4, 68.9, 137.9, 206.8, 68.9, 103.4, 206.8, 103.4, 137.9, 275.8]


def protocol() -> pd.DataFrame:
    """Nominal loading schedule. sigma_max = sigma_d + sigma_contact; sigma_contact = 10 % of sigma_max."""
    df = pd.DataFrame({"seq": range(1, 16), "s3": _SIGMA3, "smax": _SIGMAX})
    df["sd"] = 0.9 * df["smax"]
    df["sc"] = 0.1 * df["smax"]
    df["theta"] = 3 * df["s3"] + df["sd"]
    return df


def contact_from_cyclic(sd):
    """sigma_contact for a given cyclic stress (contact = 10 % of sigma_max = sd + contact)."""
    return np.asarray(sd, dtype=float) / 9.0


def load_clean(path: Path = DATA_PATH):
    """Reads the dataset and applies the notebook's screening rules (section 4.0)."""
    d = pd.read_excel(path)
    d["% CBCA"] = d["Experimento"].astype(str).str.extract(r"(\d+)")[0].astype(int)
    d["block"] = d["% CBCA"].astype(str) + "_" + d["Sequence Number"].astype(str)
    zero_strain = d["Axial Resilient def.-mm"].eq(0)
    med = d.groupby("block")[SD].transform("median")
    off_target = (d[SD] / med - 1).abs() > 0.15
    invalid = zero_strain | off_target
    clean = d.loc[~invalid].reset_index(drop=True)
    clean["theta"] = 3 * clean[S3] + clean[SD]
    return clean, int(invalid.sum())


@dataclass
class Bundle:
    data: pd.DataFrame
    n_removed: int
    rf: RandomForestRegressor
    kth: dict            # {scba: (ln k1, k2)}   Mr = exp(a) * theta^b
    sched: pd.DataFrame
    hull: Delaunay


def build_bundle(path: Path = DATA_PATH) -> Bundle:
    data, n_removed = load_clean(path)
    rf = RandomForestRegressor(**RF_PARAMS).fit(data[RF_FEATURES], data[TARGET])
    kth = {}
    for cb, g in data.groupby("% CBCA"):
        b, a = np.polyfit(np.log(g["theta"]), np.log(g[TARGET]), 1)
        kth[int(cb)] = (float(a), float(b))
    sched = protocol()
    hull = Delaunay(sched[["s3", "sd"]].to_numpy())
    return Bundle(data, n_removed, rf, kth, sched, hull)


# ------------------------------------------------------------------ predictions
def predict_rf(b: Bundle, cb, s3, sd):
    s3, sd = np.atleast_1d(np.asarray(s3, float)), np.atleast_1d(np.asarray(sd, float))
    X = pd.DataFrame({"% CBCA": np.full(len(s3), cb), S3: s3, SD: sd,
                      SC: contact_from_cyclic(sd), "theta": 3 * s3 + sd})[RF_FEATURES]
    return b.rf.predict(X)


def predict_kth(b: Bundle, cb, s3, sd):
    a, k2 = b.kth[int(cb)]
    theta = 3 * np.atleast_1d(np.asarray(s3, float)) + np.atleast_1d(np.asarray(sd, float))
    return np.exp(a + k2 * np.log(theta))


# ------------------------------------------------------------------ engineering helpers
def area_m2(diameter_mm: float) -> float:
    return np.pi * (diameter_mm / 1000.0 / 2.0) ** 2


def load_kN(stress_kpa, diameter_mm: float):
    return np.asarray(stress_kpa, float) * area_m2(diameter_mm)          # kPa * m2 = kN


def resilient_strain_microstrain(sd_kpa, mr_mpa):
    return np.asarray(sd_kpa, float) / (np.asarray(mr_mpa, float) * 1000.0) * 1e6


def in_envelope(b: Bundle, s3: float, sd: float) -> bool:
    return bool(b.hull.find_simplex([[s3, sd]])[0] >= 0)


def nearest_sequence(b: Bundle, s3: float, sd: float) -> int:
    sc = b.sched
    dist = np.hypot((sc["s3"] - s3) / sc["s3"].max(), (sc["sd"] - sd) / sc["sd"].max())
    return int(sc.loc[dist.idxmin(), "seq"])


def haversine_pulse(sd: float, t: np.ndarray, t_load: float = 0.1) -> np.ndarray:
    """Haversine cyclic stress: load for t_load seconds, then rest (period = 1 s)."""
    t = np.mod(t, 1.0)
    return np.where(t < t_load, sd * np.sin(np.pi * t / t_load) ** 2, 0.0)


def plan_table(b: Bundle, cb: int, diameter_mm: float, gauge_mm: float) -> pd.DataFrame:
    """Lab-planning table for one SCBA content. Reference = values measured in this study;
    k-theta = smooth fitted law. (The Random Forest is NOT shown here: at the 15 tested
    stress states it is in-sample and simply reproduces the measured values.)"""
    s = b.sched
    kt = predict_kth(b, cb, s["s3"], s["sd"])
    g = b.data[b.data["% CBCA"] == cb].groupby("Sequence Number")[TARGET]
    med, q1, q3 = g.median(), g.quantile(0.25), g.quantile(0.75)
    out = pd.DataFrame({
        "SCBA (%)": cb, "Secuencia": s["seq"],
        "σ3 (kPa)": s["s3"], "σd cíclico (kPa)": s["sd"], "σ contacto (kPa)": s["sc"],
        "σ máx (kPa)": s["smax"], "θ (kPa)": s["theta"],
        "Carga cíclica (kN)": load_kN(s["sd"], diameter_mm),
        "Carga contacto (kN)": load_kN(s["sc"], diameter_mm),
        "Carga máx (kN)": load_kN(s["smax"], diameter_mm),
    })
    out["Mr medido, mediana (MPa)"] = out["Secuencia"].map(med)
    out["P25–P75 entre lecturas (MPa)"] = [
        f"{q1[k]:.0f}–{q3[k]:.0f}" if k in q1.index else "—" for k in out["Secuencia"]]
    out["Mr k–θ (MPa)"] = kt
    ref = out["Mr medido, mediana (MPa)"].fillna(pd.Series(kt, index=out.index))   # fall back to k-theta if no valid data
    out["εr esperada (µε)"] = resilient_strain_microstrain(s["sd"], ref)
    out["ΔL recuperable (mm)"] = out["εr esperada (µε)"] * 1e-6 * gauge_mm
    return out.round({"σ3 (kPa)": 1, "σd cíclico (kPa)": 1, "σ contacto (kPa)": 1, "σ máx (kPa)": 1, "θ (kPa)": 1,
                      "Carga cíclica (kN)": 3, "Carga contacto (kN)": 3, "Carga máx (kN)": 3,
                      "Mr medido, mediana (MPa)": 0, "Mr k–θ (MPa)": 0,
                      "εr esperada (µε)": 0, "ΔL recuperable (mm)": 4})
