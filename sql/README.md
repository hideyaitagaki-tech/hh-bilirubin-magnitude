# SQL

BigQuery Standard SQL, written against `physionet-data.mimiciv_3_1_hosp`,
`physionet-data.mimiciv_3_1_icu` and `physionet-data.eicu_crd`.

Set `@@dataset_project_id` and `@@dataset_id` to your own project and dataset
before running, or replace `{{DATASET}}` with it. Each script writes one table.

| Script | Writes | Used for |
|--------|--------|----------|
| `00_cohort_mimic.sql` | `cohort_mimic` | derivation cohort, one row per admission |
| `01_cohort_eicu.sql` | `cohort_eicu` | validation cohort |
| `02_labs_whole_admission_mimic.sql` | `labs_mimic` | whole-admission peaks, observation density |
| `03_labs_whole_admission_eicu.sql` | `labs_eicu` | the same for validation |
| `04_cholestatic_markers_mimic.sql` | `chol_mimic` | ALP, GGT, direct bilirubin, R ratio inputs, timing anchored to the AST peak |
| `05_cholestatic_markers_eicu.sql` | `chol_eicu` | the same for validation |
| `06_landmark_competing_mimic.sql` | `landmark_mimic` | bilirubin at each landmark, time of death, competing diagnoses |
| `07_apache_competing_eicu.sql` | `apache_eicu` | APACHE IVa and its components, competing diagnoses, hospital id |
| `08_transplant_exclusion_validation.sql` | `tx_check` | post-hoc check of the transplant exclusion |
| `09_landmark_eicu.sql` | `landmark_eicu` | bilirubin at each landmark and time of death, validation cohort |
| `10_timeupdated_sofa_mimic.sql` | `sofa_landmark_mimic` | SOFA at each landmark, total and with the liver subscore removed |
| `11_alt_peak_timing_mimic.sql` | `alt_timing_mimic` | when the qualifying ALT and the AST peak fall relative to time zero |
| `12_alt_peak_timing_eicu.sql` | `alt_timing_eicu` | the same for validation |

Export each as CSV; `analysis/config.py` lists the file names the Python
expects.

## Two notes that cost us time

**eICU uses −1 as a missing-value sentinel** in `apachepatientresult` and
`apacheapsvar`. Convert it to NULL before computing anything. `07` does this in
SQL; `analysis/config.py` does it again defensively on load.

**MIMIC-IV v3.1 lives in `mimiciv_3_1_hosp`, not `mimiciv_hosp`.** The older
dataset name returns `Access Denied` rather than `Not found`, which reads like
a permissions problem and is not one. The same applies to the derived tables:
the hourly SOFA score used by `10` is in `mimiciv_3_1_derived.sofa`, and
`mimiciv_derived.sofa` returns `Access Denied`.

## A divergence to settle before release

The cohort's own time zero (`allshock_t0`, `eicu_t0`) is not the definition
used by three later extractions. In MIMIC the cohort uses itemids 220052 and
220181 with a floor of zero, while `landmark_mimic` adds 225312 and raises the
floor to 20 mmHg, and `sofa_landmark_mimic` and `alt_timing_mimic` add 225312
but keep the floor at zero. In eICU the cohort uses `vitalperiodic` alone,
while `landmark_eicu` and `alt_timing_eicu` also read `vitalaperiodic`.

Comparing the second and third MIMIC definitions, time zero differs in 21 of
1,352 admissions. The divergence from the cohort definition has not been
quantified. Before release, re-run `09` through `12` against the cohort's own
definition and confirm that the affected tables (eTables 9, 10 and 11) are
unchanged, or report the corrected values.

All other queries read the cohort tables and inherit `t0_shock` or `t0_offset`
from them, so they are unaffected.
