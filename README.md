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
protect patient information; the manuscript's data availability statement describes access. The notebook reads the files below
from `data/` (or from the folder named by the environment variable `RWE_DATA_DIR`). The files are linked by a coded patient
identifier.

| File | Contents |
|---|---|
| `rwe_uefm.csv` | UEFM assessments: the date and total score of each assessment, and which assessment is the patient's baseline |
| `rwe_patient.csv` | Patients who received a device or had a UEFM assessment: age, sex, date of stroke, affected side, Social Deprivation Index, primary payer, and the dates of prescription, device shipment, device assignment, and rental |
| `rwe_device_sessions.csv` | Therapy sessions logged by the device: the start time and duration of each session |
| `rwe_prescribed.csv` | Patients prescribed the IpsiHand system: the date of the first prescription and whether a device was received |
| `rwe_patient_core.csv` | The dominant hand of each patient |
| `mde_review_deidentified.csv` | The clinical review of significant detrimental medical events, which identifies the patients excluded for such an event |
| `rwe_prescribed_dropoff.csv` | The primary payer recorded for each prescribed patient |
| `prescribed_dropoff_reasons_2026-09-14.csv` | For each prescribed patient, whether a device was received and, if not, the reason |
| `recipient_flow_step_eeg_status_2026-09-14.csv` | The result of each device recipient's pre-therapy electroencephalography (EEG) screen |

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
