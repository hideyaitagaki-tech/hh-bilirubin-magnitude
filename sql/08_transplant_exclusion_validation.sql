-- =====================================================================
-- Post-hoc validation of the liver-transplant exclusion (MIMIC-IV)
--
-- Asks whether the aminotransferase elevation in the excluded admissions
-- arose from the transplant itself rather than from shock, by comparing
-- the date of the transplant procedure with the timing of the AST peak.
--
-- Input: the admissions excluded by 00_cohort_mimic.sql, i.e. those in
-- the pre-exclusion cohort but not the final one.
--
-- Writes: {{DATASET}}.tx_check
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.tx_check` AS
WITH ids AS (
  SELECT hadm_id FROM `{{DATASET}}.cohort_mimic_preexclusion`
  EXCEPT DISTINCT
  SELECT hadm_id FROM `{{DATASET}}.cohort_mimic`
),
adm AS (
  SELECT a.hadm_id, a.admittime, a.dischtime
  FROM `physionet-data.mimiciv_3_1_hosp.admissions` a
  JOIN ids i USING (hadm_id)
),
tx_proc AS (
  SELECT p.hadm_id, p.icd_code, p.icd_version, p.chartdate,
         DATE_DIFF(p.chartdate, DATE(a.admittime), DAY) AS tx_day_from_admit
  FROM `physionet-data.mimiciv_3_1_hosp.procedures_icd` p
  JOIN adm a USING (hadm_id)
  WHERE REGEXP_CONTAINS(p.icd_code, r'^(5051|5059|0FY00Z)')
),
tx_dx AS (
  SELECT d.hadm_id, STRING_AGG(DISTINCT d.icd_code ORDER BY d.icd_code) AS tx_dx_codes
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d
  JOIN ids i USING (hadm_id)
  WHERE REGEXP_CONTAINS(d.icd_code, r'^(V427|Z944|9960|99682|T8641|T8649)')
  GROUP BY d.hadm_id
),
ast AS (
  SELECT l.hadm_id,
         MAX(l.valuenum) AS peak_ast,
         DATETIME_DIFF(
           ARRAY_AGG(l.charttime ORDER BY l.valuenum DESC LIMIT 1)[OFFSET(0)],
           a.admittime, HOUR) AS ast_peak_hr_from_admit
  FROM `physionet-data.mimiciv_3_1_hosp.labevents` l
  JOIN adm a USING (hadm_id)
  WHERE l.itemid = 50878 AND l.valuenum IS NOT NULL AND l.valuenum > 0
  GROUP BY l.hadm_id, a.admittime
)
SELECT
  i.hadm_id,
  p.icd_code  AS tx_procedure_code,
  p.chartdate AS tx_date,
  p.tx_day_from_admit,
  s.peak_ast,
  s.ast_peak_hr_from_admit,
  ROUND(s.ast_peak_hr_from_admit / 24.0, 1) AS ast_peak_day_from_admit,
  x.tx_dx_codes,
  DATE_DIFF(DATE(a.dischtime), DATE(a.admittime), DAY) AS los_days
FROM ids i
LEFT JOIN adm     a USING (hadm_id)
LEFT JOIN tx_proc p USING (hadm_id)
LEFT JOIN tx_dx   x USING (hadm_id)
LEFT JOIN ast     s USING (hadm_id)
ORDER BY p.tx_day_from_admit, i.hadm_id;
