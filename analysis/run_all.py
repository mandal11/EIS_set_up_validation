#!/usr/bin/env python3
"""Regenerate every number and figure reported in the paper.

    python3 analysis/run_all.py

Writes the derived tables to data/derived/ and the figures to figures/, and
prints the values quoted in the abstract, in Tables I, IV to VI and IX, and in
the conclusions.
"""
import json
import numpy as np
import pandas as pd
from eis import CELLS, BANDS, DERIVED, load, summary, r_squared, lin_kk
from temperature import temperature_table, arrhenius_fit


def band_and_full_metrics():
    rows, agg = [], {}
    for cell in CELLS:
        d = load(cell)
        agg[cell] = {}
        for label, lo, hi in list(BANDS) + [("FULL SWEEP", 0.0, 1e9)]:
            m = (d.f >= lo) & (d.f < hi)
            s = d[m]
            ar = summary(s.e_re, s.ref_re)
            ai = summary(s.e_im, s.ref_im)
            rows.append(dict(
                Cell=cell, Band=label, n=int(m.sum()),
                Zre_RMSE=ar["RMSE"], Zre_MAE=ar["MAE"], Zre_Max=ar["MaxAbs"],
                Zim_RMSE=ai["RMSE"], Zim_MAE=ai["MAE"], Zim_Max=ai["MaxAbs"],
                NormRes_re=float((s.e_re / s.ref_mag).abs().mean() * 100),
                NormRes_im=float((s.e_im / s.ref_mag).abs().mean() * 100)))
            if label == "FULL SWEEP":
                agg[cell] = {"Zreal": ar, "Zimg": ai}
                agg[cell]["Zreal"]["R2"] = r_squared(s.ref_re, s.dut_re)
                agg[cell]["Zimg"]["R2"] = r_squared(s.ref_im, s.dut_im)
                agg[cell]["mag_pct"] = float(
                    (100 * (s.dut_mag - s.ref_mag) / s.ref_mag).abs().mean())
                agg[cell]["phase_deg"] = float((s.dut_ph - s.ref_ph).abs().mean())
    t = pd.DataFrame(rows)
    t.to_csv(DERIVED / "band_summary.csv", index=False)
    return t, agg


def kk_metrics():
    out = {}
    for cell in CELLS:
        d = load(cell)
        out[cell] = {}
        for tag, re_, im_ in (("DUT", d.dut_re, d.dut_im), ("REF", d.ref_re, d.ref_im)):
            Z = re_.values + 1j * im_.values
            fit, M = lin_kk(d.f.values, Z)
            rr = (Z.real - fit.real) / np.abs(Z)
            ri = (Z.imag - fit.imag) / np.abs(Z)
            out[cell][tag] = dict(
                M=M,
                sigma_re_pct=float(np.std(rr, ddof=1) * 100),
                sigma_im_pct=float(np.std(ri, ddof=1) * 100),
                sd_abs_re=float(np.std(Z.real - fit.real, ddof=1)),
                sd_abs_im=float(np.std(Z.imag - fit.imag, ddof=1)))
    json.dump(out, open(DERIVED / "linkk.json", "w"), indent=1)
    return out


def decomposition(agg, kk):
    rows = []
    for cell in CELLS:
        for comp, key, tag in (("Zreal", "Zreal", "re"), ("Zimg", "Zimg", "im")):
            a = agg[cell][key]
            sd, sr = kk[cell]["DUT"][f"sd_abs_{tag}"], kk[cell]["REF"][f"sd_abs_{tag}"]
            comb = float(np.hypot(sd, sr))
            rows.append(dict(Cell=cell, Component=comp, RMSE=a["RMSE"], MAE=a["MAE"],
                             Max=a["MaxAbs"], Bias=a["Bias"], SD_obs=a["SD"],
                             sigma_DUT=sd, sigma_REF=sr, sigma_comb=comb,
                             excess=float(np.sqrt(max(a["SD"] ** 2 - comb ** 2, 0.0))),
                             R2=a["R2"]))
    t = pd.DataFrame(rows)
    t.to_csv(DERIVED / "error_decomposition.csv", index=False)
    return t


def per_frequency_tables():
    for cell in CELLS:
        d = load(cell)
        t = pd.DataFrame({
            "f_Hz": d.f, "Zre_DUT": d.dut_re, "Zre_REF": d.ref_re, "dZre": d.e_re,
            "Zim_DUT": d.dut_im, "Zim_REF": d.ref_im, "dZim": d.e_im,
            "absZ_DUT": d.dut_mag, "absZ_REF": d.ref_mag,
            "ph_DUT": d.dut_ph, "ph_REF": d.ref_ph, "dph_deg": d.dut_ph - d.ref_ph})
        t["res_re_pct"] = 100 * t.dZre / t.absZ_REF
        t["res_im_pct"] = 100 * t.dZim / t.absZ_REF
        t.round(4).to_csv(DERIVED / f"SI_full_error_{cell}.csv", index=False)


def main():
    pd.set_option("display.width", 200)
    bands, agg = band_and_full_metrics()
    kk = kk_metrics()
    dec = decomposition(agg, kk)
    per_frequency_tables()

    print("=" * 78)
    print("Table V, decade-resolved and full-sweep error metrics (mOhm)")
    print(bands.round(3).to_string(index=False))

    print("\n" + "=" * 78)
    print("Abstract, Table I and conclusions")
    for cell in CELLS:
        a = agg[cell]
        print(f"  {cell:9s} Zreal RMSE={a['Zreal']['RMSE']:.2f}  Zimg RMSE={a['Zimg']['RMSE']:.2f}"
              f"   |Z| error={a['mag_pct']:.2f}%   phase error={a['phase_deg']:.2f} deg"
              f"   R2={a['Zreal']['R2']:.4f}/{a['Zimg']['R2']:.4f}")

    print("\n" + "=" * 78)
    print("Kramers-Kronig residual standard deviations (% of |Z|)")
    for cell in CELLS:
        for tag, name in (("DUT", "this work"), ("REF", "EC301   ")):
            v = kk[cell][tag]
            print(f"  {cell:9s} {name}  real={v['sigma_re_pct']:.3f}  imag={v['sigma_im_pct']:.3f}")

    print("\n" + "=" * 78)
    print("Table VI, decomposition of the inter-instrument difference (mOhm)")
    print(dec.round(3).to_string(index=False))
    print("\n" + "=" * 78)
    print("Table IX, ohmic intercept and polarisation resistance from the temperature series (mOhm)")
    t = temperature_table()
    t.round(3).to_csv(DERIVED / "temperature_table.csv", index=False)
    print(t.round(2).to_string(index=False))
    for cell in sorted(set(t.cell)):
        s = t[t.cell == cell]
        for key in ("R_ohm", "R_pol"):
            ea, r2, _ = arrhenius_fit(s.T_C.values, s[key].values)
            ratio = s[key].values[0] / s[key].values[-1]
            print(f"  {cell:9s} {key:5s}  Ea = {1000 * ea:4.0f} meV   R2 = {r2:.3f}   "
                  f"ratio -10/50 degC = {ratio:.2f}")

    print("\nWritten to", DERIVED)


if __name__ == "__main__":
    main()
