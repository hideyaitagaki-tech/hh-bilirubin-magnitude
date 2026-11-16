# Delayed hyperbilirubinemia in hypoxic hepatitis — analysis code

Code and cohort logic for:

> Itagaki H, Endo T. *Delayed, Cholestatic Hyperbilirubinemia in Hypoxic
> Hepatitis: Magnitude Predicts Mortality in Two Critical Care Databases.*

Derivation cohort: MIMIC-IV v3.1. Validation cohort: eICU Collaborative
Research Database v2.0. Both are accessed on Google BigQuery.

## What is here, and what is not

This repository contains the **code**: the SQL that defines the cohorts and
extracts every variable, the Python that produces every number, table and
figure in the manuscript, and prose documentation of the variable
derivations, the exclusion logic and the definition of time zero.

It does **not** contain patient data, and it does not contain the lists of
`hadm_id` or `patientunitstayid` that make up the cohorts. MIMIC-IV and
eICU-CRD are credentialed-access datasets and the PhysioNet Credentialed
Health Data Use Agreement does not permit redistribution. A list of the
admissions in a cohort is derived data from those datasets, so it is not
published here either. Running `sql/` against the source tables regenerates
the cohorts exactly.

## Reproducing the study

1. Obtain credentialed access to MIMIC-IV v3.1 and eICU-CRD v2.0 through
   PhysioNet and complete the required training.
2. Run the scripts in `sql/` in order. Each writes one table into your own
   BigQuery dataset; `sql/README.md` lists them and what each produces.
3. Export those tables as CSV into `data/` (see `analysis/config.py` for the
   expected file names).
4. `pip install -r requirements.txt`
5. Run the scripts in `analysis/` in order. They write the tables and figures
   into `output/`.

The analysis scripts read only the exported CSVs, so step 5 is reproducible
without database access once step 3 has been done by a credentialed user.

## Layout

```
sql/        cohort definition and variable extraction (BigQuery Standard SQL)
analysis/   statistical analysis, tables and figures (Python)
  exploratory/  ad-hoc checks run during development; nothing reported
                depends on them (see its README)
docs/       variable derivations, exclusion logic, time zero, diagnosis codes
```

## Software

Python 3.11 with pandas, numpy, scipy, statsmodels, scikit-learn,
matplotlib and lifelines; versions in `requirements.txt`. The regression
models are fitted with `statsmodels`; the time-varying Cox model uses
`lifelines`.

## Licence

Code is released under the MIT Licence (`LICENSE`). The licence covers this
code only; it confers no rights in MIMIC-IV or eICU-CRD.
