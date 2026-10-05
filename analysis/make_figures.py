#!/usr/bin/env python3
"""Generate the figures added to the paper at revision.

    python3 analysis/make_figures.py

Writes PDF and PNG files to figures/:
  capability_map    Fig. 1   frequency range and lowest measured impedance of reported systems
  bandwidth_design  Fig. 5   limits of the measurement band
  err_full_sweep    Fig. 12  deviation from the reference at all 72 frequencies
  band_rmse         Fig. 13  RMSE by decade band
  linkk_residuals   Fig. 14  linear Kramers-Kronig residuals of both instruments
  arrhenius         Fig. 16  ohmic and polarisation resistance against 1/T
  nyquist_*         Fig. 11  Nyquist plots of the validation sweeps
  nyquist_temp_*    Fig. 15  Nyquist plots of the temperature series
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from eis import CELLS, BANDS, FIGURES, RAW, load, lin_kk, summary
from temperature import temperature_table, arrhenius_fit

# Okabe-Ito palette, readable with the common forms of colour-vision deficiency
C = dict(p="#0072B2", r="#D55E00", g="#009E73", v="#CC79A7", o="#E69F00",
         grid="#DDDDDD", ink="#1A1A1A", mute="#666666")
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"], "font.size": 9,
    "mathtext.fontset": "stix", "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 7.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.edgecolor": "#444444",
    "axes.linewidth": 0.8, "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.5,
    "lines.linewidth": 1.4, "lines.markersize": 4.0, "figure.dpi": 150,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "text.color": C["ink"], "axes.labelcolor": C["ink"],
    "xtick.color": C["mute"], "ytick.color": C["mute"]})
NAME = {"Samsung": "Samsung INR18650-35E", "LGM50": "LG INR21700-M50"}


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIGURES / f"{name}.pdf", metadata={"CreationDate": None})
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    plt.close(fig)
    print("  figures/" + name)


def capability_map():
    # name, f_min (Hz), f_max (Hz), lowest measured impedance, validated on commercial cells
    rows = [("Juang 2023", 8e3, 1e4, "1 k$\\Omega$", False),
            ("Jenkins 2019", 60, 2e3, "100 $\\Omega$", False),
            ("Burgos 2022", 1, 1e4, "100 $\\Omega$", False),
            ("Matsubara 2021", 0.1, 1e6, "400 $\\Omega$", False),
            ("Grassini 2015", 1e-3, 1e5, "100 $\\Omega$", False),
            ("Zhang 2024", 0.1, 8e2, "m$\\Omega$ range", True),
            ("Carbone 2023", 0.05, 1e2, "m$\\Omega$ range", True),
            ("Wu 2024", 0.1, 2e3, "m$\\Omega$ range", True),
            ("Yue 2026", 1e-2, 5.62e3, "20 m$\\Omega$", True),
            ("This work", 2e-2, 1e4, "30 m$\\Omega$", True)]
    fig, ax = plt.subplots(figsize=(6.7, 3.6))
    ax.axvspan(2e-2, 1e4, color=C["p"], alpha=0.06, zorder=0)
    for i, (lab, f0, f1, zmin, cell) in enumerate(rows):
        this = lab == "This work"
        col = C["p"] if this else (C["g"] if cell else C["mute"])
        ax.plot([f0, f1], [i, i], "-", color=col, lw=5 if this else 3.2,
                solid_capstyle="butt", zorder=4)
        ax.plot([f0, f1], [i, i], "|", color=col, ms=9, mew=1.5, zorder=5)
        ax.text(1.5e6, i, zmin, fontsize=7.5, va="center", color=col,
                fontweight="bold" if this else "normal")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    for t, r in zip(ax.get_yticklabels(), rows):
        if r[0] == "This work":
            t.set_fontweight("bold")
            t.set_color(C["p"])
    ax.set_xscale("log")
    ax.set_xlim(5e-4, 1.2e6)
    ax.set_xlabel("Frequency range of reported measurements (Hz)")
    ax.grid(axis="y", visible=False)
    ax.text(1.5e6, len(rows) - 0.45, "Lowest\nimpedance", fontsize=7.5, va="bottom", color=C["ink"])
    ax.set_ylim(-0.7, len(rows) + 0.4)
    ax.legend(handles=[Patch(color=C["p"], label="This work"),
                       Patch(color=C["g"], label="Tested on commercial lithium-ion cells"),
                       Patch(color=C["mute"], label="Tested on resistive or electrode loads")],
              frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=7.3)
    save(fig, "capability_map")


def bandwidth_design():
    f = np.logspace(-3, 6, 3000)
    r_loop = 15.0
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.8))
    a = ax[0]
    for cs, ls in ((0.05, ":"), (0.53, "--"), (3.7, "-")):
        x = 2j * np.pi * f * cs * r_loop
        fl = 1 / (2 * np.pi * cs * r_loop)
        a.semilogx(f, 20 * np.log10(np.abs(x / (1 + x))), ls, color=C["p"],
                   label=rf"$C_s$ = {cs:g} F, $f_L$ = {1e3 * fl:.3g} mHz")
    a.axvspan(2e-2, 1e4, color=C["g"], alpha=0.08)
    a.axhline(-3, color=C["mute"], lw=0.8, ls="-.")
    a.set_ylim(-30, 2)
    a.set_xlabel("Frequency (Hz)")
    a.set_ylabel(r"$|I_{\rm cell}|/|I_{\rm cell}|_{f\to\infty}$ (dB)")
    a.set_title(r"(a) Blocking network, $R_{\rm loop}$ = 15 $\Omega$", loc="left")
    a.legend(frameon=False, loc="lower right", fontsize=6.8)
    a = ax[1]
    fc = 1e4
    s = 1j * f / fc
    filt = 1 / (s ** 2 + np.sqrt(2) * s + 1)
    stage = 1 / (1 + 1j * f / 5e5)
    a.semilogx(f, 20 * np.log10(np.abs(filt)), "-", color=C["r"],
               label="second-order sensing filter, $f_c$ = 10 kHz")
    a.semilogx(f, 20 * np.log10(np.abs(stage)), "--", color=C["p"],
               label="TL064 stage, gain 2 (500 kHz)")
    a.axvspan(2e-2, 1e4, color=C["g"], alpha=0.08)
    a.axhline(-3, color=C["mute"], lw=0.8, ls="-.")
    a.set_ylim(-30, 2)
    a.set_xlabel("Frequency (Hz)")
    a.set_ylabel("Normalised gain (dB)")
    a.set_title("(b) Sensing filter and drive stage", loc="left")
    a.legend(frameon=False, loc="lower left", fontsize=6.8)
    save(fig, "bandwidth_design")


def err_full_sweep():
    fig, ax = plt.subplots(2, 2, figsize=(7.1, 4.8))
    for k, cell in enumerate(CELLS):
        d = load(cell)
        a = ax[0, k]
        a.axhspan(-1, 1, color=C["g"], alpha=0.07, zorder=0)
        a.axhline(0, color=C["mute"], lw=0.7)
        a.semilogx(d.f, d.e_re, "o-", color=C["p"], mfc="white", mew=1.0, label=r"$\varepsilon'$ (real part)")
        a.semilogx(d.f, d.e_im, "s-", color=C["r"], mfc="white", mew=1.0, label=r"$\varepsilon''$ (imaginary part)")
        a.set_title(f"({'ab'[k]}) {NAME[cell]}", loc="left")
        a.set_ylabel(r"Deviation from EC301 (m$\Omega$)")
        a.set_xlim(0.015, 1.3e4)
        lo, hi = a.get_ylim()
        a.set_ylim(lo, hi + 0.35 * (hi - lo))
        a.legend(frameon=False, loc="upper left", ncol=2)
        a = ax[1, k]
        a.axhline(0, color=C["mute"], lw=0.7)
        a.semilogx(d.f, 100 * (d.dut_mag - d.ref_mag) / d.ref_mag, "o-", color=C["g"], mfc="white",
                   mew=1.0, label=r"modulus error (%)")
        a.semilogx(d.f, d.dut_ph - d.ref_ph, "^-", color=C["v"], mfc="white", mew=1.0,
                   label=r"phase error ($^\circ$)")
        a.set_title(f"({'cd'[k]}) {NAME[cell]}", loc="left")
        a.set_xlabel("Frequency (Hz)")
        a.set_ylabel("Error (% or $^\\circ$)")
        a.set_xlim(0.015, 1.3e4)
        lo, hi = a.get_ylim()
        a.set_ylim(lo, hi + 0.35 * (hi - lo))
        a.legend(frameon=False, loc="upper left", ncol=2)
    save(fig, "err_full_sweep")


def band_rmse():
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.7), sharey=True)
    labels = ["20 mHz\nto 1 Hz", "1 Hz\nto 100 Hz", "100 Hz\nto 1 kHz", "1 kHz\nto 10 kHz"]
    x = np.arange(len(BANDS))
    w = 0.36
    for k, cell in enumerate(CELLS):
        d = load(cell)
        re, im = [], []
        for _, lo, hi in BANDS:
            m = (d.f >= lo) & (d.f < hi)
            re.append(summary(d.e_re[m])["RMSE"])
            im.append(summary(d.e_im[m])["RMSE"])
        a = ax[k]
        a.bar(x - w / 2, re, w, color=C["p"], label="real part", zorder=3, edgecolor="white")
        a.bar(x + w / 2, im, w, color=C["r"], label="imaginary part", zorder=3, edgecolor="white")
        for xi, v in zip(x - w / 2, re):
            a.text(xi, v + 0.03, f"{v:.2f}", ha="center", fontsize=6.8)
        for xi, v in zip(x + w / 2, im):
            a.text(xi, v + 0.03, f"{v:.2f}", ha="center", fontsize=6.8)
        a.set_xticks(x)
        a.set_xticklabels(labels, fontsize=7)
        a.set_title(f"({'ab'[k]}) {NAME[cell]}", loc="left")
        a.grid(axis="x", visible=False)
        a.set_ylim(0, 1.65)
        if k == 0:
            a.set_ylabel(r"RMSE from EC301 (m$\Omega$)")
            a.legend(frameon=False, ncol=2, loc="upper left")
    save(fig, "band_rmse")


def linkk_residuals():
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.7), sharey=True)
    for k, cell in enumerate(CELLS):
        d = load(cell)
        a = ax[k]
        a.axhline(0, color=C["mute"], lw=0.7)
        for tag, re_, im_, col, mk, lab in (("p", d.dut_re, d.dut_im, C["p"], "o", "This work"),
                                            ("r", d.ref_re, d.ref_im, C["r"], "s", "EC301")):
            z = re_.values + 1j * im_.values
            fit, _ = lin_kk(d.f.values, z)
            a.semilogx(d.f, 100 * (z.real - fit.real) / np.abs(z), mk, color=col, mfc="white",
                       mew=1.0, label=f"{lab}, real")
            a.semilogx(d.f, 100 * (z.imag - fit.imag) / np.abs(z), mk, color=col, alpha=0.5,
                       label=f"{lab}, imaginary")
        a.set_title(f"({'ab'[k]}) {NAME[cell]}", loc="left")
        a.set_xlabel("Frequency (Hz)")
        a.set_xlim(0.015, 1.3e4)
        if k == 0:
            a.set_ylabel(r"Residual (% of $|Z|$)")
            a.legend(frameon=False, ncol=2, loc="upper left", fontsize=7)
    save(fig, "linkk_residuals")


def arrhenius():
    t = temperature_table()
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.8))
    for cell, col, mk in (("Samsung", C["p"], "o"), ("Panasonic", C["r"], "s")):
        s = t[t.cell == cell]
        x = 1000 / (s.T_C + 273.15)
        for a, key in ((ax[0], "R_ohm"), (ax[1], "R_pol")):
            ea, r2, fit = arrhenius_fit(s.T_C.values, s[key].values)
            inside = ~s.at_edge.values
            a.semilogy(x[inside], s[key][inside], mk, color=col, mfc="white", mew=1.2)
            a.semilogy(x[~inside], s[key][~inside], mk, color=col, mfc=col, mew=1.2)
            xx = np.linspace(x.min(), x.max(), 50)
            a.semilogy(xx, fit(xx), "-", color=col, lw=1.1,
                       label=f"{cell}: $E_a$ = {1000 * ea:.0f} meV ($R^2$ = {r2:.2f})")
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, NullLocator, FormatStrFormatter
    edge = Line2D([], [], ls="none", marker="o", color=C["mute"], mfc=C["mute"],
                  label="sweep ended before the sign change")
    for a, title, lab, ticks, lim in ((ax[0], r"(a) Ohmic intercept $R_\Omega$", r"$R_\Omega$ (m$\Omega$)",
                                       [30, 35, 40, 45, 50, 55, 60, 70], (30, 78)),
                                      (ax[1], r"(b) Polarisation resistance $R_{\rm pol}$",
                                       r"$R_{\rm pol}$ (m$\Omega$)", [5, 10, 20, 50, 100, 200, 500], (4, 900))):
        a.set_ylim(*lim)
        a.set_xlabel(r"$1000/T$ (K$^{-1}$)")
        a.set_ylabel(lab)
        a.set_title(title, loc="left", pad=22)
        a.yaxis.set_major_locator(FixedLocator(ticks))
        a.yaxis.set_minor_locator(NullLocator())
        a.yaxis.set_major_formatter(FormatStrFormatter("%g"))
        h, l = a.get_legend_handles_labels()
        a.legend(handles=h + [edge], frameon=False, loc="upper left", fontsize=6.8)
        top = a.secondary_xaxis("top", functions=(lambda v: 1000 / np.maximum(v, 1e-9) - 273.15,
                                                   lambda c: 1000 / (c + 273.15)))
        top.set_xticks([50, 25, 0, -10])
        top.set_xlabel(r"Temperature ($^\circ$C)", fontsize=8, labelpad=2)
    save(fig, "arrhenius")


def nyquist_validation():
    """Fig. 11: Nyquist plots of the validation sweeps, one file per cell."""
    lims = {"Samsung": ((42, 66), (-8, 9)), "LGM50": ((29, 44), (-4, 6))}
    for cell, tag in (("Samsung", "samsung"), ("LGM50", "lg")):
        d = load(cell)
        fig, a = plt.subplots(figsize=(3.4, 2.9))
        a.axhline(0, color=C["mute"], lw=0.7)
        a.plot(d.ref_re, -d.ref_im, "s-", color=C["r"], mfc="white", mew=1.0, ms=3.6, lw=1.0,
               label="EC301 potentiostat")
        a.plot(d.dut_re, -d.dut_im, "o-", color=C["p"], mfc="white", mew=1.0, ms=3.4, lw=1.0,
               label="proposed front-end")
        (x0, x1), (y0, y1) = lims[cell]
        a.set_xlim(x0, x1)
        a.set_ylim(y0, y1)
        a.set_aspect("equal", adjustable="box")
        a.set_xlabel(r"$Z'$ (m$\Omega$)")
        a.set_ylabel(r"$-Z''$ (m$\Omega$)")
        a.set_title(NAME[cell], loc="left")
        a.legend(frameon=False, loc="upper left", fontsize=7)
        hidden = int(((-d.dut_im) < y0).sum())
        print(f"    {cell}: {hidden} front-end points below the plotted range (f >= {d.f[(-d.dut_im) < y0].min():g} Hz)")
        save(fig, f"nyquist_{tag}")


def nyquist_temperature():
    """Fig. 15: Nyquist plots of the temperature series, one file per cell."""
    import matplotlib.cm as cm
    temps = ["-10", "0", "15", "25", "40", "50"]
    cols = dict(zip(temps, ["#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7"]))
    marks = dict(zip(temps, ["o", "s", "^", "D", "v", "P"]))
    for cell, tag in (("Samsung", "samsung"), ("Panasonic", "panasonic")):
        sheets = pd.read_excel(RAW / f"{cell}.xlsx", sheet_name=None)
        fig, a = plt.subplots(figsize=(3.4, 2.7))
        a.axhline(0, color=C["mute"], lw=0.7)
        for t in temps:
            s = sheets[t].sort_values("Frequency")
            a.plot(s.Zreal, s.Zimg, marks[t] + "-", color=cols[t], mfc="white", mew=0.9, ms=1.5,
                   lw=0.9, label=f"${t}$ $^\\circ$C")
        a.set_xlim(25, 200)
        a.set_ylim(-5, 45)
        a.set_aspect("equal", adjustable="box")
        a.set_xlabel(r"$Z'$ (m$\Omega$)")
        a.set_ylabel(r"$-Z''$ (m$\Omega$)")
        a.set_title(cell, loc="left")
        a.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.55), fontsize=6.8, ncol=6,
                 handlelength=1.2, columnspacing=0.8)
        save(fig, f"nyquist_temp_{tag}")


if __name__ == "__main__":
    FIGURES.mkdir(exist_ok=True)
    capability_map()
    bandwidth_design()
    err_full_sweep()
    band_rmse()
    linkk_residuals()
    arrhenius()
    nyquist_validation()
    nyquist_temperature()
