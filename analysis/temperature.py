"""Ohmic and polarisation resistance from the temperature series (Section VI-B8, Table IX).

The temperature files store the real part and the NEGATIVE imaginary part of the
impedance (column `Zimg` = -Im Z, positive in the capacitive region), in milliohms.
"""
import numpy as np
import pandas as pd

from eis import RAW

K_B = 8.617333e-5  # eV/K
T_CELLS = ("Samsung", "Panasonic")


def crossing(f, zre, nzim):
    """Real part and frequency at which -Im Z first changes sign (capacitive to inductive).

    Linear interpolation between the two adjacent points. If the sweep ends before the
    sign change, the last point is returned and `at_edge` is True: R_ohm is then an upper
    bound and the crossing frequency a lower bound.
    """
    for i in range(len(f) - 1):
        if nzim[i] > 0 >= nzim[i + 1]:
            t = nzim[i] / (nzim[i] - nzim[i + 1])
            return zre[i] + t * (zre[i + 1] - zre[i]), f[i] + t * (f[i + 1] - f[i]), False
    return zre[-1], f[-1], True


def temperature_table():
    rows = []
    for cell in T_CELLS:
        sheets = pd.read_excel(RAW / f"{cell}.xlsx", sheet_name=None)
        for name, s in sheets.items():
            s = s.sort_values("Frequency")
            f, zre, nzim = s.Frequency.values, s.Zreal.values, s.Zimg.values
            r_ohm, f_x, edge = crossing(f, zre, nzim)
            z_dc = zre[np.argmin(f)]
            rows.append(dict(cell=cell, T_C=float(name), f_max=float(f.max()), R_ohm=r_ohm,
                             Zre_20mHz=z_dc, R_pol=z_dc - r_ohm, f_cross=f_x, at_edge=edge))
    return pd.DataFrame(rows).sort_values(["cell", "T_C"]).reset_index(drop=True)


def arrhenius_fit(t_c, y):
    """Apparent activation energy (eV) from ln y = Ea/(kB T) + c, with R^2 and the fitted curve
    as a function of 1000/T."""
    x = 1.0 / (K_B * (np.asarray(t_c, float) + 273.15))
    ly = np.log(np.asarray(y, float))
    p = np.polyfit(x, ly, 1)
    r2 = 1 - np.sum((ly - np.polyval(p, x)) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return p[0], r2, lambda inv_t_k: np.exp(np.polyval(p, inv_t_k / 1000.0 / K_B))
