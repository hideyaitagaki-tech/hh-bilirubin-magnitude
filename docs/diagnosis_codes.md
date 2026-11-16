# Diagnosis codes and diagnosis strings

MIMIC-IV records diagnoses as ICD-9 and ICD-10 codes; eICU-CRD records them as
free-text `diagnosisstring` values. The two are therefore matched differently,
and the manuscript states that under-ascertainment is unavoidable in both.

Matching is case-insensitive throughout. In eICU-CRD, `diagnosisstring` is a
hierarchical path such as
`gastrointestinal|hepatic disease|cirrhosis`, so a substring match on the leaf
term also catches the parent categories that contain it.

---

## Liver transplantation (exclusion)

Recovered verbatim from the BigQuery job that created the cohort table.

**MIMIC-IV, `procedures_icd`** — exact codes, not a prefix match

```sql
WHERE (icd_version = 10 AND icd_code IN ('0FY00Z0','0FY00Z1'))
   OR (icd_version = 9  AND icd_code IN ('5059','5051'))
```

| Code | Meaning |
|------|---------|
| `0FY00Z0` | ICD-10-PCS, transplantation of liver, allogeneic |
| `0FY00Z1` | ICD-10-PCS, transplantation of liver, syngeneic |
| `5051` | ICD-9-CM 50.51, auxiliary liver transplant |
| `5059` | ICD-9-CM 50.59, other transplant of liver |

Only **procedure** codes were used. An earlier version of this file also listed
diagnosis codes for a transplanted liver (V42.7, Z94.4) and transplant
complications; those were never applied, and the entry has been corrected
against the job that was run.

**eICU-CRD** — `LOWER(diagnosisstring) LIKE '%liver transplant%'`, nothing else.

---

## Competing causes of hyperbilirubinemia (sensitivity analyses)

These are **not** exclusion criteria. See `docs/exclusions.md`.

### eICU-CRD — `diagnosis.diagnosisstring`

The regular expressions below are the ones used, verbatim:

| Variable | Pattern |
|----------|---------|
| `dx_chronic_liver` | `cirrhosis\|chronic liver\|portal hypertension\|esophageal varice` |
| `dx_biliary_obstruction` | `cholangitis\|biliary obstruction\|common bile duct\|choledocholithiasis\|biliary tract` |
| `dx_hemolysis` | `hemolysis\|hemolytic` |
| `dx_rhabdomyolysis` | `rhabdomyolysis` |

`dx_hemolysis` matched no admission in the validation cohort. This is reported
in the manuscript as evidence that diagnostic entry in eICU-CRD is incomplete,
not as a prevalence of zero.

### MIMIC-IV — `diagnoses_icd.icd_code`

Recovered verbatim from the BigQuery job. Each flag is 1 if **any** diagnosis
of the admission matches. Matching is a prefix regular expression applied
within the stated ICD version.

| Variable | ICD-10 | ICD-9 |
|----------|--------|-------|
| `dx_chronic_liver` | `^(K703\|K704\|K717\|K721\|K74\|K766\|K767\|B18\|I85)` | `^(571\|5723\|4560\|4561\|4562\|07032\|07044\|07054)` |
| `dx_biliary_obstruction` | `^(K80\|K830\|K831\|K82\|C24\|C221)` | `^(5740\|5743\|5745\|5748\|5761\|5762\|1561)` |
| `dx_hemolysis` | `^(D55\|D56\|D57\|D58\|D59)` | `^(282\|283)` |
| `dx_rhabdomyolysis` | `^(M6282\|T795)` | `^(72888\|9584)` |

Chronic liver disease covers alcoholic and other fibrosis and cirrhosis (K70.3,
K70.4, K71.7, K72.1, K74), portal hypertension and hepatorenal syndrome (K76.6,
K76.7), chronic viral hepatitis (B18) and oesophageal varices (I85); the ICD-9
set is the corresponding 571.x, 572.3, 456.0–456.2 and the chronic viral
hepatitis codes. Biliary obstruction covers cholelithiasis and its
complications (K80), cholangitis and biliary obstruction (K83.0, K83.1), other
gallbladder disease (K82) and biliary malignancy (C24, C22.1).

---

## A code set that was removed

An earlier draft also extracted `dx_acute_liver_failure` from ICD codes for
acute liver failure and hepatic necrosis (ICD-10 `^(K72)`, ICD-9 `^(5701|570)`). It
matched 60.4% of the derivation cohort, because those codes **are the index
condition**: hypoxic hepatitis is itself coded as hepatic necrosis. Adjusting
for it produced non-converging models. The variable is not a competing cause
and was dropped. It is recorded here so that anyone repeating the analysis
does not reintroduce it.

---

## Shock

Recovered verbatim from the jobs that built `allshock_cohort` and `eicu_t0`.

### MIMIC-IV — `diagnoses_icd.icd_code`, exact match

| Etiology | Codes |
|----------|-------|
| Septic | `R6521`, `78552`, `R6520`, `T8112XA`, `99802` |
| Cardiogenic | `R570`, `78551`, `T8111XA`, `99801` |
| Hypovolemic | `R571`, `T794XXA`, `9584` |

`R65.21` septic shock, `785.52` septic shock, `R65.20` severe sepsis without
shock, `T81.12XA` and `998.02` postprocedural septic shock; `R57.0` and
`785.51` cardiogenic shock, `T81.11XA` and `998.01` postprocedural cardiogenic
shock; `R57.1` hypovolemic shock, `T79.4XXA` traumatic shock, `958.4` traumatic
shock. Matching is on the exact code, not a prefix.

### eICU-CRD — `diagnosis.diagnosisstring`, LIKE patterns

Rows are first filtered to `'%shock%' OR '%|sepsis%'`, then typed:

| Etiology | Pattern |
|----------|---------|
| Cardiogenic | `%cardiogenic shock%` |
| Septic | `%septic shock%` or `%\|sepsis%` |
| Hypovolemic | `%hypovolemi%` |

### The exclusive hierarchy

An admission with codes for more than one etiology is assigned a single type,
in the order **septic > cardiogenic > hypovolemic**, in both cohorts. This is
applied after typing, by taking the highest-priority type present:

```sql
CASE
  WHEN MAX(CASE WHEN shock_type='sepsis'      THEN 1 ELSE 0 END)=1 THEN 'sepsis'
  WHEN MAX(CASE WHEN shock_type='cardiogenic' THEN 1 ELSE 0 END)=1 THEN 'cardiogenic'
  ELSE 'hypovolemic'
END
```

Note that the eICU row-level `CASE` tests cardiogenic before septic, so a
single diagnosis string naming both is typed cardiogenic at the row level; the
hierarchy above then restores septic priority at the admission level.

---

## Inclusion on ALT

In both cohorts the ALT criterion was applied over the **whole stay**, with no
requirement that the value follow shock onset:

```sql
-- MIMIC
SELECT DISTINCT hadm_id FROM labevents WHERE itemid = 50861 AND valuenum > 500
-- eICU
SELECT DISTINCT patientunitstayid FROM lab
WHERE labname = 'ALT (SGPT)' AND labresult > 500
```

This is the operationalization that the referees identified, and it is why the
manuscript now reports how many admissions qualified only on a pre-shock value
and repeats the analysis with the qualifying value required to fall at or after
time zero.
