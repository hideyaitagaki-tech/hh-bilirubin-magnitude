-- =====================================================================
-- Whole-stay laboratory peaks and observation density (eICU-CRD)
-- Writes: {{DATASET}}.labs_eicu
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.labs_eicu` AS
WITH coh AS (
  SELECT DISTINCT patientunitstayid
  FROM `{{DATASET}}.cohort_eicu`
),
lab AS (
  SELECT l.patientunitstayid, l.labname, l.labresult, l.labresultoffset
  FROM `physionet-data.eicu_crd.lab` l
  JOIN coh c USING (patientunitstayid)
  WHERE l.labname IN ('ALT (SGPT)', 'AST (SGOT)', 'LDH')
    AND l.labresult IS NOT NULL
    AND l.labresult > 0
)
SELECT
  c.patientunitstayid,
  MAX(IF(l.labname = 'ALT (SGPT)', l.labresult, NULL)) AS alt_peak,
  MAX(IF(l.labname = 'AST (SGOT)', l.labresult, NULL)) AS ast_peak,
  MAX(IF(l.labname = 'LDH',        l.labresult, NULL)) AS ldh_peak,
  MIN(l.labresultoffset) / 60.0                        AS first_lab_hr,
  MAX(l.labresultoffset) / 60.0                        AS last_lab_hr,
  COUNTIF(l.labname = 'ALT (SGPT)')                    AS n_alt
FROM coh c
LEFT JOIN lab l USING (patientunitstayid)
GROUP BY c.patientunitstayid
ORDER BY c.patientunitstayid;
