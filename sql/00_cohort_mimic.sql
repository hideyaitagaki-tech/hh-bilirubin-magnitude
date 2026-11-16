-- =====================================================================
-- Derivation cohort (MIMIC-IV v3.1)
--
-- RECOVERED VERBATIM from the BigQuery jobs of 2026-08-31 and 2026-09-01.
-- Four steps, run in this order. Code lists and the exclusive hierarchy are
-- explained in docs/diagnosis_codes.md; time zero in docs/time_zero.md.
-- =====================================================================

-- ---------- 1. allshock_cohort: shock type, ALT criterion, baseline ----------
CREATE OR REPLACE TABLE `{{DATASET}}.allshock_cohort` AS
WITH shock_typed AS (
  SELECT hadm_id,
    CASE
      WHEN icd_code IN ('R6521','78552','R6520','T8112XA','99802') THEN 'sepsis'
      WHEN icd_code IN ('R570','78551','T8111XA','99801')          THEN 'cardiogenic'
      WHEN icd_code IN ('R571','T794XXA','9584')                   THEN 'hypovolemic'
    END AS shock_type
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE icd_code IN ('R6521','78552','R6520','T8112XA','99802',
                     'R570','78551','T8111XA','99801',
                     'R571','T794XXA','9584')
),
-- one etiology per admission: sepsis > cardiogenic > hypovolemic
typed_unique AS (
  SELECT hadm_id,
    CASE
      WHEN MAX(CASE WHEN shock_type='sepsis' THEN 1 ELSE 0 END)=1 THEN 'sepsis'
      WHEN MAX(CASE WHEN shock_type='cardiogenic' THEN 1 ELSE 0 END)=1 THEN 'cardiogenic'
      ELSE 'hypovolemic'
    END AS shock_type
  FROM shock_typed WHERE shock_type IS NOT NULL
  GROUP BY hadm_id
),
-- hypoxic hepatitis: ALT > 500 U/L anywhere in the admission
high_alt AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.labevents`
  WHERE itemid = 50861 AND valuenum > 500
),
sofa1 AS (
  SELECT icu.hadm_id, fs.sofa,
         ROW_NUMBER() OVER (PARTITION BY icu.hadm_id ORDER BY icu.intime) AS rn
  FROM `physionet-data.mimiciv_3_1_icu.icustays` icu
  JOIN `physionet-data.mimiciv_3_1_derived.first_day_sofa` fs ON icu.stay_id = fs.stay_id
)
SELECT
  t.hadm_id,
  t.shock_type,
  p.anchor_age AS age,
  p.gender,
  a.hospital_expire_flag AS hospital_death,
  s.sofa
FROM typed_unique t
JOIN high_alt h ON t.hadm_id = h.hadm_id
JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON t.hadm_id = a.hadm_id
JOIN `physionet-data.mimiciv_3_1_hosp.patients` p ON a.subject_id = p.subject_id
LEFT JOIN sofa1 s ON t.hadm_id = s.hadm_id AND s.rn = 1;

-- ---------- 2. allshock_t0: time zero ----------
-- NOTE: two itemids only, and no floor beyond zero. See docs/time_zero.md.
CREATE OR REPLACE TABLE `{{DATASET}}.allshock_t0` AS
SELECT
  c.hadm_id,
  c.shock_type,
  MIN(ce.charttime) AS t0_shock
FROM `{{DATASET}}.allshock_cohort` c
JOIN `physionet-data.mimiciv_3_1_icu.chartevents` ce
  ON ce.hadm_id = c.hadm_id
WHERE ce.itemid IN (220052, 220181)
  AND ce.valuenum < 65
  AND ce.valuenum > 0
GROUP BY c.hadm_id, c.shock_type;

-- ---------- 3. allshock_bili: peak bilirubin and its timing ----------
CREATE OR REPLACE TABLE `{{DATASET}}.allshock_bili` AS
WITH bili_ts AS (
  SELECT
    t.hadm_id,
    t.shock_type,
    le.valuenum AS bili,
    TIMESTAMP_DIFF(le.charttime, t.t0_shock, HOUR) AS hours_from_shock,
    ROW_NUMBER() OVER (PARTITION BY t.hadm_id ORDER BY le.valuenum DESC) AS rn
  FROM `{{DATASET}}.allshock_t0` t
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le
    ON le.hadm_id = t.hadm_id
  WHERE le.itemid = 50885
    AND le.valuenum IS NOT NULL
    AND le.charttime >= t.t0_shock
)
SELECT hadm_id, shock_type, bili AS peak_bili, hours_from_shock AS bili_peak_hours
FROM bili_ts
WHERE rn = 1;

-- =====================================================================
-- RECOVERED, VERBATIM: the final step of the derivation cohort.
--
-- Retrieved from the BigQuery job of 2026-09-01. This is the step that
-- applies the liver-transplant exclusion and adds the dialysis flag. Its
-- input, `allshock_final`, is the upstream step that applies the ALT and
-- shock criteria and the exclusive etiology hierarchy; that job is still
-- to be recovered (see sql/README.md).
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.allshock_pure` AS
WITH liver_transplant AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.procedures_icd`
  WHERE (icd_version = 10 AND icd_code IN ('0FY00Z0','0FY00Z1'))
     OR (icd_version = 9  AND icd_code IN ('5059','5051'))
),
dialysis AS (
  SELECT DISTINCT icu.hadm_id
  FROM `physionet-data.mimiciv_3_1_icu.procedureevents` pe
  JOIN `physionet-data.mimiciv_3_1_icu.icustays` icu ON pe.stay_id = icu.stay_id
  WHERE pe.itemid IN (225802, 225441, 225805)
)
SELECT
  f.*,
  CASE WHEN dz.hadm_id IS NOT NULL THEN 1 ELSE 0 END AS dialysis
FROM `{{DATASET}}.allshock_final` f
LEFT JOIN liver_transplant lt ON f.hadm_id = lt.hadm_id
LEFT JOIN dialysis dz ON f.hadm_id = dz.hadm_id
WHERE lt.hadm_id IS NULL;
