# Milliohm impedance measurement of commercial lithium-ion cells

Data, analysis scripts and circuit files for:

> S. Mandal, B. Sah, S. Pandey, S. K. Mulpuri and P. Kumar,
> "Reference-Validated Milliohm Impedance Measurement of Commercial Lithium-Ion Cells
> Using a Compact Current-Excitation Analogue Front-End",
> *Energy Storage*, manuscript 3584481 (under review).

The error metrics, the Kramers-Kronig analysis, the temperature table and the figures added to
the paper at revision can be reproduced from the files in this repository.

```bash
pip install -r requirements.txt
python3 analysis/run_all.py        # error metrics, Kramers-Kronig test, temperature table
python3 analysis/make_figures.py   # figures added at revision, written to figures/
```

`run_all.py` prints the values quoted in the abstract, in Tables I, IV to VI and IX and in the
conclusions, and writes the derived tables to `data/derived/`.

## Contents

```
data/raw/          measured impedance spectra
  Samsung_pointwise.csv    Samsung INR18650-35E, 72 frequencies, 20 mHz to 10 kHz
  LGM50_pointwise.csv      LG INR21700-M50, 72 frequencies, 20 mHz to 10 kHz
  Samsung.xlsx             temperature series, -10 to 50 degC, one sheet per temperature
  Panasonic.xlsx           temperature series, -10 to 50 degC, one sheet per temperature
data/derived/      tables written by run_all.py
analysis/          eis.py (shared functions), temperature.py, run_all.py, make_figures.py
hardware/          schematic, circuit-board views, bill of materials, LTspice model
figures/           the figures added at revision, as PDF and PNG
```

## Column definitions

`data/raw/*_pointwise.csv`, impedances in milliohms, phases in degrees:

| column | meaning |
|---|---|
| `f` | excitation frequency, Hz |
| `dc_re`, `dc_im` | real and imaginary impedance measured with the analogue front-end |
| `ec_re`, `ec_im` | real and imaginary impedance measured with the EC301 potentiostat |
| `dc_ph`, `ec_ph` | impedance phase |
| `dc_nim`, `ec_nim` | negative imaginary part, as plotted in a Nyquist diagram |

A capacitive response has `im < 0`. The change of sign near 1 kHz marks the transition to
inductive behaviour.

The temperature files (`Samsung.xlsx`, `Panasonic.xlsx`) have the columns `Frequency` (Hz),
`Zreal` and `Zimg` (milliohms). Note that `Zimg` in these files is the negative imaginary part,
positive in the capacitive region.

## Measurement frequencies

The 72 frequencies are spaced by 10 mHz from 20 mHz to 0.1 Hz, by 0.1 Hz up to 1 Hz, by 1 Hz up
to 10 Hz, by 10 Hz up to 100 Hz, by 50 Hz up to 1 kHz and by 500 Hz up to 10 kHz. The same cell
was measured with both instruments at the same frequencies.

## Methods

`analysis/eis.py` contains the two methods used in the paper.

**Error metrics.** At a single frequency the root-mean-square error and the mean absolute error
are both equal to the absolute deviation. The tables therefore give the absolute deviation at
single frequencies, and RMSE and MAE only for sets of frequencies (Eqs. (15) to (19) of the
paper).

**Linear Kramers-Kronig test** (`lin_kk`). A Voigt model with a series resistance, a series
inductance and M parallel RC elements with fixed, logarithmically spaced time constants is fitted
by modulus-weighted linear least squares. The model satisfies the Kramers-Kronig relations for
any parameter values, so the scatter of the residuals estimates the random error of the
measurement without an assumed equivalent circuit. A linear, time-invariant systematic error,
such as a gain error or a series lead inductance, also satisfies the Kramers-Kronig relations and
is not detected by this test. The systematic part is estimated separately by the uncertainty
budget in Table VII of the paper.

Method references: B. A. Boukamp, *J. Electrochem. Soc.* **142**, 1885 (1995);
M. Schonleber, D. Klotz and E. Ivers-Tiffee, *Electrochim. Acta* **131**, 20 (2014);
P. Agarwal, O. D. Crisalle, M. E. Orazem and L. H. Garcia-Rubio,
*J. Electrochem. Soc.* **142**, 4149 (1995).

## Still to be added

- [ ] five-run repeat data for the LGM50 cell (Section VI-B6 of the paper, SI Fig. S3)
- [ ] value of the DC-blocking capacitors; the LTspice model uses 4.7 uF and has no 10 ohm shunt
- [ ] bill of materials: the TL064CN line carries the part number of a DIP socket, and the
      amplifiers and the DC-blocking capacitors are missing
- [ ] Altium source files of the circuit board

## Licence

Code in `analysis/` is MIT (see `LICENSE`). Data, hardware files and figures are CC BY 4.0
(see `LICENSE-DATA`).

## Citation

See `CITATION.cff`. Please cite both the paper and the archived release.
