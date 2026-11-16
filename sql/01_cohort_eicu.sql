-- =====================================================================
-- Validation cohort (eICU-CRD v2.0)
--
-- RECOVERED VERBATIM from the BigQuery jobs of 2026-08-31 and 2026-09-01.
-- =====================================================================

-- ---------- 1. eicu_t0: shock type, ALT criterion, time zero, baseline -------
-- NOTE: time zero uses vitalperiodic only. See docs/time_zero.md.
CREATE OR REPLACE TABLE `{{DATASET}}.eicu_t0` AS
WITH shock_typed AS (
  SELECT DISTINCT patientunitstayid,
    CASE
      WHEN LOWER(diagnosisstring) LIKE '%cardiogenic shock%' THEN 'cardiogenic'
      WHEN LOWER(diagnosisstring) LIKE '%septic shock%'
        OR LOWER(diagnosisstring) LIKE '%|sepsis%' THEN 'sepsis'
      WHEN LOWER(diagnosisstring) LIKE '%hypovolemi%' THEN 'hypovolemic'
    END AS shock_type
  FROM `physionet-data.eicu_crd.diagnosis`
  WHERE LOWER(diagnosisstring) LIKE '%shock%'
     OR LOWER(diagnosisstring) LIKE '%|sepsis%'
),
typed_unique AS (
  SELECT patientunitstayid,
    CASE
      WHEN MAX(CASE WHEN shock_type='sepsis' THEN 1 ELSE 0 END)=1 THEN 'sepsis'
      WHEN MAX(CASE WHEN shock_type='cardiogenic' THEN 1 ELSE 0 END)=1 THEN 'cardiogenic'
      ELSE 'hypovolemic'
    END AS shock_type
  FROM shock_typed WHERE shock_type IS NOT NULL
  GROUP BY patientunitstayid
),
high_alt AS (
  SELECT DISTINCT patientunitstayid
  FROM `physionet-data.eicu_crd.lab`
  WHERE labname = 'ALT (SGPT)' AND labresult > 500
),
cohort AS (
  SELECT t.patientunitstayid, t.shock_type
  FROM typed_unique t
  JOIN high_alt h ON t.patientunitstayid = h.patientunitstayid
),
t0 AS (
  SELECT patientunitstayid, MIN(observationoffset) AS t0_offset
  FROM `physionet-data.eicu_crd.vitalperiodic`
  WHERE systemicmean IS NOT NULL AND systemicmean > 0 AND systemicmean < 65
  GROUP BY patientunitstayid
)
SELECT
  c.patientunitstayid, c.shock_type, t0.t0_offset,
  p.age, p.gender, p.unitdischargestatus, p.hospitaldischargestatus
FROM cohort c
JOIN t0 ON c.patientunitstayid = t0.patientunitstayid   -- requires a time zero
JOIN `physionet-data.eicu_crd.patient` p ON c.patientunitstayid = p.patientunitstayid;

-- ---------- 2. eicu_bili: peak bilirubin and its timing ----------
CREATE OR REPLACE TABLE `{{DATASET}}.eicu_bili` AS
WITH bili_ts AS (
  SELECT
    c.patientunitstayid,
    c.shock_type,
    l.labresult AS bili,
    (l.labresultoffset - c.t0_offset) / 60.0 AS hours_from_shock,
    ROW_NUMBER() OVER (PARTITION BY c.patientunitstayid ORDER BY l.labresult DESC) AS rn
  FROM `{{DATASET}}.eicu_t0` c
  JOIN `physionet-data.eicu_crd.lab` l
    ON l.patientunitstayid = c.patientunitstayid
   AND l.labname = 'total bilirubin'
   AND l.labresult IS NOT NULL
   AND l.labresultoffset >= c.t0_offset
)
SELECT patientunitstayid, shock_type, bili AS peak_bili,
       hours_from_shock AS bili_peak_hours,
       CASE WHEN bili >= 10 THEN 1 ELSE 0 END AS bili_ge10
FROM bili_ts
WHERE rn = 1;

-- =====================================================================
-- RECOVERED, VERBATIM: the final step of the validation cohort.
--
-- Retrieved from the BigQuery job of 2026-09-01. Its input,
-- `eicu_extended`, is the upstream step that applies the ALT and shock
-- criteria and the exclusive etiology hierarchy; that job is still to be
-- recovered (see sql/README.md).
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.eicu_pure` AS
WITH liver_transplant AS (
  SELECT DISTINCT patientunitstayid
  FROM `physionet-data.eicu_crd.diagnosis`
  WHERE LOWER(diagnosisstring) LIKE '%liver transplant%'
)
SELECT e.*
FROM `{{DATASET}}.eicu_extended` e
LEFT JOIN liver_transplant lt ON e.patientunitstayid = lt.patientunitstayid
WHERE lt.patientunitstayid IS NULL;
