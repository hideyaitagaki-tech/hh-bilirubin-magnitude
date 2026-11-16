# Variable derivations

All laboratory variables are **peaks over the whole hospital admission**
unless stated otherwise. An earlier draft of this work computed them over a
restricted window; that was corrected, and the correction changed several
reported descriptive values (see `docs/exclusions.md`, "Corrections").

## Exposure

**Marked hyperbilirubinemia** — peak total bilirubin ≥10 mg/dL, defined a
priori. In the categorical dose–response analysis the same variable is cut at
`[0, 3)`, `[3, 10)`, `[10, ∞)`; note the half-open intervals, so a value of
exactly 10.0 falls in the exposed category and is consistent with the binary
exposure.

## Outcome

**In-hospital mortality.** MIMIC-IV: `admissions.hospital_expire_flag`.
eICU-CRD: `patient.hospitaldischargestatus = 'Expired'`.

## Laboratory variables

### MIMIC-IV — `mimiciv_3_1_hosp.labevents`

| Variable | itemid |
|----------|--------|
| Alanine aminotransferase (ALT) | 50861 |
| Alkaline phosphatase (ALP) | 50863 |
| Aspartate aminotransferase (AST) | 50878 |
| Bilirubin, direct | 50883 |
| Bilirubin, total | 50885 |
| Gamma-glutamyl transferase (GGT) | 50927 |
| Lactate dehydrogenase (LDH) | 50954 |

Rows with NULL or non-positive `valuenum` are discarded.

### eICU-CRD — `eicu_crd.lab`

Matched on `labname`:

| Variable | labname |
|----------|---------|
| ALT | `ALT (SGPT)` |
| ALP | `alkaline phos.` |
| AST | `AST (SGOT)` |
| Bilirubin, direct | `direct bilirubin` |
| Bilirubin, total | `total bilirubin` |
| LDH | `LDH` |

GGT is not recorded in eICU-CRD. It was measured in only 81 of 1,352
derivation admissions (6.0%) and was not analysed.

## Derived biochemical variables

**AST/ALT ratio** — peak AST divided by peak ALT, both over the whole
admission.

**R ratio** — `(ALT / ALT-ULN) ÷ (ALP / ALP-ULN)`, computed from the values
measured **at the time of the bilirubin peak**, not from the peaks. Upper
limits of normal: ALT 33 U/L (women) and 45 U/L (men); ALP 105 and 130 U/L.
Conventional bands: R ≤2 cholestatic, 2–5 mixed, ≥5 hepatocellular.

**Direct/total bilirubin ratio** — direct bilirubin at the bilirubin peak
divided by the peak total bilirubin. Ratios outside `[0, 1.2]` are treated as
recording errors and dropped; the 1.2 bound allows for the two values coming
from different assays on the same draw.

**Timing relative to the AST peak** — hours from the AST peak to the
bilirubin peak (`hours_ast_to_bili`) and to the ALP peak (`hours_ast_to_alp`).
Negative values mean the marker peaked first. This is the anchoring used by
Van den Broecke et al.

## Severity

**SOFA (derivation)** — first-day SOFA. Note that SOFA contains a bilirubin
component; because the exposure is a late peak and the score is taken on day
one, the two are only weakly correlated (Spearman ρ = 0.261). The manuscript
reports the model with SOFA omitted alongside the primary model.

**APACHE IVa (validation)** — `eicu_crd.apachepatientresult`, restricted to
`apacheversion = 'IVa'`, giving `apachescore` and `predictedhospitalmortality`.
Physiological components come from `eicu_crd.apacheapsvar`: `bilirubin`,
`creatinine`, `meanbp`, `ph`, `pao2`, `intubated`, `vent`, `dialysis`. In both
tables **−1 is a missing-value sentinel** and must be converted to NULL before
use; failing to do so silently biases every summary of these variables.

## Other variables

**Transfusion (derivation)** — volumes of packed red blood cells and of all
blood products, from `mimiciv_3_1_icu.inputevents`. In eICU-CRD transfusion is
recorded only as unit categories, so volumes are not comparable and are not
reported.

**Contributing hospitals (validation)** — `patient.hospitalid`.

**Observation density** — the number of ALT measurements per admission and the
interval between the first and last laboratory result, reported to show that
the two databases sample at different rates.

## Recovered from the cohort-building jobs

The following derivations were retrieved verbatim from the BigQuery jobs that
built `allshock_final` and `eicu_extended`, and are recorded here because the
cohort tables carry the resulting columns.

**Enzyme peaks in the cohort table are post-shock.** `peak_ast`, `peak_alt` and
`peak_ldh` in `allshock_final` were taken from `charttime >= t0_shock`
(itemids 50878, 50861, 50954). These are *not* the whole-admission peaks used
in the manuscript; those are `alt_peak_all`, `ast_peak_all` and `ldh_peak_all`
from `labs_mimic` (`sql/02`). Both sets of columns exist in the exported data,
and `analysis/` uses the whole-admission ones throughout.

**Exposure.** `bili_ge10 = CASE WHEN peak_bili >= 10 THEN 1 ELSE 0 END`,
applied in SQL at cohort construction.

**Cohort filter.** Both cohorts keep only admissions with
`bili_peak_hours >= 0`, that is, a bilirubin measurement at or after time zero;
the validation cohort additionally requires a non-null
`hospitaldischargestatus`.

| Variable | Derivation |
|----------|------------|
| `peak_creat` | labevents itemid 50912, `valuenum < 20`, at or after time zero |
| `peak_inr` | labevents itemid 51237, `valuenum < 30`, at or after time zero |
| `prbc_ml` | inputevents itemid 225168, summed from time zero |
| `total_tx_ml` | inputevents itemids 225168, 220970, 225170, 225171, summed from time zero |
| `icu_los_days` | `los` of the first ICU stay of the admission, by `intime` |
| `dialysis` | procedureevents itemids 225802 (CRRT), 225441 (haemodialysis), 225805 (peritoneal dialysis) |
| `peak_dbili` (eICU) | labname `direct bilirubin`, at or after time zero |
| `massive_transfusion` (eICU) | treatmentstring containing `transfusion of > 2 units prbc` |
| `any_transfusion` (eICU) | treatmentstring containing `packed red blood cells` |

The ceilings on creatinine and INR (`< 20`, `< 30`) are value filters applied
at extraction, not winsorisation of the analysis variables.
