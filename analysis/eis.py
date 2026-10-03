"""Shared helpers: loading, error metrics and the linear Kramers-Kronig test.

All impedances are in milliohms, all frequencies in hertz. The suffix `dut`
refers to the analogue front-end reported in the paper and `ref` to the EC301.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
FIGURES = ROOT / "figures"
CELLS = ("Samsung", "LGM50")

BANDS = (("20 mHz -- 1 Hz", 0.0, 1.0),
         ("1 -- 100 Hz", 1.0, 100.0),
         ("100 Hz -- 1 kHz", 100.0, 1000.0),
         ("1 -- 10 kHz", 1000.0, 10001.0))


def load(cell):
    """Return the validation sweep for one cell as a DataFrame."""
    d = pd.read_csv(RAW / f"{cell}_pointwise.csv")
    d = d.rename(columns={"dc_re": "dut_re", "dc_im": "dut_im", "dc_ph": "dut_ph",
                          "ec_re": "ref_re", "ec_im": "ref_im", "ec_ph": "ref_ph"})
    d["dut_mag"] = np.hypot(d.dut_re, d.dut_im)
    d["ref_mag"] = np.hypot(d.ref_re, d.ref_im)
    d["e_re"] = d.dut_re - d.ref_re
    d["e_im"] = d.dut_im - d.ref_im
    return d


def rmse(e):
    """Root-mean-square error over a set of frequencies, Eq. (21)."""
    e = np.asarray(e, float)
    return float(np.sqrt(np.mean(e ** 2)))


def mae(e):
    """Mean absolute error over a set of frequencies, Eq. (21)."""
    return float(np.mean(np.abs(np.asarray(e, float))))


def summary(e, ref=None):
    """RMSE, MAE, maximum absolute deviation, bias and scatter, Eqs. (21)-(22)."""
    e = np.asarray(e, float)
    out = {"N": len(e), "RMSE": rmse(e), "MAE": mae(e),
           "MaxAbs": float(np.max(np.abs(e))), "Bias": float(np.mean(e)),
           "SD": float(np.std(e, ddof=1))}
    if ref is not None:
        out["MAPE_%"] = float(np.mean(np.abs(e / np.asarray(ref, float))) * 100)
    return out


def r_squared(y, yhat):
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    return float(1 - np.sum((y - yhat) ** 2) / np.sum((y - np.mean(y)) ** 2))


def lin_kk(f, Z, M=None):
    """Linear Kramers-Kronig (Voigt) fit, Boukamp (1995).

    Fits Z = R0 + jwL + sum_k Rk / (1 + jw tau_k) with the tau_k fixed and
    logarithmically spaced across the measured range, so the problem is linear
    and is solved in one modulus-weighted least-squares step. Such a model is
    Kramers-Kronig consistent by construction, so the residual scatter about it
    estimates the random error of the measurement.
    """
    w = 2 * np.pi * np.asarray(f, float)
    Z = np.asarray(Z, complex)
    N = len(w)
    M = M or max(6, int(N * 0.6))
    tau = np.logspace(np.log10(1 / w.max()), np.log10(1 / w.min()), M)
    A = np.column_stack([np.ones(N, complex), 1j * w] + [1.0 / (1.0 + 1j * w * t) for t in tau])
    wt = 1.0 / np.abs(Z)
    Ar = np.vstack([A.real * wt[:, None], A.imag * wt[:, None]])
    br = np.concatenate([Z.real * wt, Z.imag * wt])
    p, *_ = np.linalg.lstsq(Ar, br, rcond=None)
    return A @ p.astype(complex), M
