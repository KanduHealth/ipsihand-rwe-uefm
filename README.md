# Upper Extremity Fugl-Meyer outcomes with home-based BCI therapy in routine care

Analysis code for the manuscript *A Retrospective Real-World Assessment of Post-Stroke Motor Function Improvements Using a
Home-Based Brain Computer Interface*. The study is a retrospective, single-arm cohort of adults prescribed the IpsiHand
brain-computer interface (BCI) system in routine care; the outcome is the Upper Extremity Fugl-Meyer (UEFM) total score.

## Contents

```
rwe_uefm_analysis.ipynb   the analysis, one section per part of the Results, run top to bottom
R/mmrm_functions.R        mixed model for repeated measures (MMRM): model, covariance structure, contrasts
R/mmrm_fit.R              primary MMRM                     Rscript R/mmrm_fit.R <input.csv> <output_prefix>
R/mmrm_mi.R               multiple-imputation sensitivity analyses (carry-forward imputation; tipping point)
requirements.txt          Python packages
```

| Notebook section | Manuscript | Output (in `outputs/`) |
|---|---|---|
| 1 | Methods: patient variables, device use, follow-up weeks, study population | — |
| 2 | Figure 1B, patient flow | `figure1b_patient_flow.csv`, `figures/figure1b_patient_flow.*` |
| 2 | Table 1, patient characteristics | `table1_characteristics.csv` |
| 2 | Supplemental Table 1, follow-up UEFM availability | `supplemental_table1_followup_availability.csv` |
| 2 | Supplemental Table 2, primary payer type | `supplemental_table2_payer.csv` |
| 2 | Supplemental Table 3, active users at each assessment window | `supplemental_table3_active_users.csv` |
| 3 | Change from baseline at each follow-up week (Results 3.2 and 3.3) | `paired_change_by_week.csv` |
| 3 | Figure 2, change after 12 weeks | `figures/figure2_week12.*` |
| 3 | Supplemental Figure 1, baseline characteristics by MCID response | `supplemental_figure1_mcid_baseline.csv`, `figures/supplemental_figure1_mcid_baseline.*` |
| 3 | Supplemental Table 4, baseline characteristics with and without follow-up | `supplemental_table4_followup_comparison.csv` |
| 4 | Supplemental Figure 2, observed change at each follow-up week | `figures/supplemental_figure2_observed_change.*` |
| 4 | Supplemental Figures 3 and 4, change after 16 and 24 weeks | `figures/supplemental_figure3_week16.*`, `figures/supplemental_figure4_week24.*` |
| 4 | Figure 3, MMRM change over one year | `mmrm_lsmeans.csv`, `mmrm_baseline_coef.csv`, `figures/figure3_mmrm_trajectory.*` |
| 4 | Supplemental Figure 5, carry-forward multiple imputation | `mi_carry_forward_lsmeans.csv`, `figures/supplemental_figure5_carry_forward_mi.*` |
| 4 | Supplemental Figure 6, tipping-point analysis | `tipping_point.csv`, `figures/supplemental_figure6_tipping_point.*` |
| 4 | Exploratory analysis, change at weeks 12, 16, and 24 for a baseline UEFM of 20 to 50 (Results 3.3) | `paired_change_baseline_20_50.csv` |

MCID: minimal clinically important difference (5.25 UEFM points). Figures are written as PNG and SVG. Reported values are
rounded half up.

## Requirements

Python 3.9 or later with the packages in `requirements.txt`, and R 4.x with the `mmrm` and `jsonlite` packages (`Rscript` on the
`PATH`).

## Data

The data are patient-level records from the program's clinical, device, and payer systems, extracted on September 6, 2026.
Direct identifiers were removed; dates are retained, so the files form a limited data set. They are not publicly available to
protect patient information; the manuscript's data availability statement describes access. The notebook reads the files below from `data/` (or from the folder named by the environment variable
`RWE_DATA_DIR`). Every file is keyed on `PATIENT_ID`, a coded patient key; only the listed columns are used.

| File | Columns used |
|---|---|
| `rwe_uefm.csv` | `PATIENT_ID`, `SCORE_TIMELINE` (`baseline` marks the baseline assessment), `SCORE_DATE`, `UEFM` |
| `rwe_patient.csv` | `PATIENT_ID`, `AGE`, `GENDER`, `DATE_OF_STROKE`, `DATE_OF_STROKE_PRECISION` (`day`, `month`, or `year`), `SIDE_AFFECTED`, `SDI_SCORE` (Social Deprivation Index), `PAYER_TYPE_PRIMARY`, `FIRST_RX_RECEIVED_DATE`, `FIRST_DEVICE_SHIP_DATE`, `FIRST_SALES_ORDER_SHIP_DATE`, `ANY_SALES_ORDER_DELIVERED`, `FIRST_DEVICE_ASSIGNED_DATE`, `LAST_DEVICE_ASSOCIATION_ENDED_DATE`, `N_OPEN_DEVICE_ASSOCIATIONS`, `LAST_RENTAL_END_DATE`, `N_ACTIVE_RENTAL_SESSIONS` |
| `rwe_device_sessions.csv` | `PATIENT_ID`, `SESSION_START`, `DURATION_SECONDS` |
| `rwe_prescribed.csv` | `PATIENT_ID`, `FIRST_RX_RECEIVED_DATE`, `DEVICE_RECEIVED` |
| `rwe_patient_core.csv` | `PATIENT_ID`, `DOMINANT_HAND` |
| `mde_review_deidentified.csv` | `PATIENT_ID`, `MDE?` (clinical review of significant detrimental medical events: `No MDE present`, blank when not yet reviewed, otherwise the confirmed event) |
| `rwe_prescribed_dropoff.csv` | `PATIENT_ID`, `PAYER_TYPE_PRIMARY_LIVE` (payer on the insurance case), `PAYER_TYPE_OPPORTUNITY` (payer on the case record) |
| `prescribed_dropoff_reasons_2026-09-14.csv` | `PATIENT_ID`, `REASON_BIN`: why a prescribed patient did not receive a device (`A` payer denial or coverage barrier; `B` screening or clinical ineligibility; `C` patient withdrew or was non-responsive; `D` documentation or administrative step not completed; `E` case open without a payer decision; `F` closed without a recorded reason; `R` device received) |
| `recipient_flow_step_eeg_status_2026-09-14.csv` | `PATIENT_ID`, `EEG_SCREEN_STATUS`: pre-therapy EEG screen of device recipients (`F` failed; `X` failed, then passed a retest; `P` passed; `H` held without a recorded result; `N` no screen on record) |

## Run

```
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m ipykernel install --user --name rwe-uefm
.venv/bin/jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=rwe-uefm --ExecutePreprocessor.timeout=3600 rwe_uefm_analysis.ipynb
```

Outputs go to `outputs/` (or the folder named by `RWE_OUTPUT_DIR`). Besides the tables and figures listed above, `outputs/`
holds the two patient-level files the notebook passes to R (`mmrm_input.csv`, `mmrm_impute_cells.csv`), so keep it as private as
the input data; `.gitignore` excludes `data/` and `outputs/`. The tipping-point analysis is the slow step (about 15 minutes). To
refit the primary model alone after a run: `Rscript R/mmrm_fit.R outputs/mmrm_input.csv outputs/mmrm`.

## License

The code is released under the MIT License; see `LICENSE`.
