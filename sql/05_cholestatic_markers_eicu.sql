-- =====================================================================
-- Cholestatic markers and timing anchored to the AST peak (eICU-CRD)
--
-- GGT is not recorded in eICU-CRD and is therefore absent here.
-- Writes: {{DATASET}}.chol_eicu
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.chol_eicu` AS
WITH coh AS (
  SELECT patientunitstayid, shock_type, age, bili_ge10, hospital_death
  FROM `{{DATASET}}.cohort_eicu`
),
pt AS (
  SELECT p.patientunitstayid, p.gender
  FROM `physionet-data.eicu_crd.patient` p
  JOIN coh c USING (patientunitstayid)
),
lab AS (
  SELECT l.patientunitstayid, l.labname, l.labresult, l.labresultoffset
  FROM `physionet-data.eicu_crd.lab` l
  JOIN coh c USING (patientunitstayid)
  WHERE l.labname IN ('ALT (SGPT)', 'alkaline phos.', 'AST (SGOT)',
                      'direct bilirubin', 'total bilirubin')
    AND l.labresult IS NOT NULL
    AND l.labresult > 0
),
pk AS (
  SELECT patientunitstayid,
    MAX(IF(labname = 'total bilirubin',  labresult, NULL)) AS v_peak_bili,
    MAX(IF(labname = 'AST (SGOT)',       labresult, NULL)) AS v_peak_ast,
    MAX(IF(labname = 'alkaline phos.',   labresult, NULL)) AS v_peak_alp,
    MAX(IF(labname = 'direct bilirubin', labresult, NULL)) AS peak_dbili,
    COUNTIF(labname = 'alkaline phos.')                    AS n_alp,
    COUNTIF(labname = 'direct bilirubin')                  AS n_dbili
  FROM lab GROUP BY patientunitstayid
),
t_ast AS (
  SELECT patientunitstayid,
         ARRAY_AGG(labresultoffset ORDER BY labresult DESC, labresultoffset LIMIT 1)[OFFSET(0)] AS t
  FROM lab WHERE labname = 'AST (SGOT)' GROUP BY patientunitstayid
),
t_bili AS (
  SELECT patientunitstayid,
         ARRAY_AGG(labresultoffset ORDER BY labresult DESC, labresultoffset LIMIT 1)[OFFSET(0)] AS t
  FROM lab WHERE labname = 'total bilirubin' GROUP BY patientunitstayid
),
t_alp AS (
  SELECT patientunitstayid,
         ARRAY_AGG(labresultoffset ORDER BY labresult DESC, labresultoffset LIMIT 1)[OFFSET(0)] AS t
  FROM lab WHERE labname = 'alkaline phos.' GROUP BY patientunitstayid
),
at_peak AS (
  SELECT b.patientunitstayid,
    MAX(IF(l.labname = 'alkaline phos.',   l.labresult, NULL)) AS alp_at_bili_peak,
    MAX(IF(l.labname = 'ALT (SGPT)',       l.labresult, NULL)) AS alt_at_bili_peak,
    MAX(IF(l.labname = 'direct bilirubin', l.labresult, NULL)) AS dbili_at_bili_peak
  FROM t_bili b
  JOIN lab l ON l.patientunitstayid = b.patientunitstayid
            AND l.labresultoffset = b.t
  GROUP BY b.patientunitstayid
)
SELECT
  c.*, pt.gender,
  p.v_peak_bili, p.v_peak_ast, p.v_peak_alp, p.peak_dbili, p.n_alp, p.n_dbili,
  a.alp_at_bili_peak, a.alt_at_bili_peak, a.dbili_at_bili_peak,
  (tb.t - ta.t) / 60.0 AS hours_ast_to_bili,
  (tp.t - ta.t) / 60.0 AS hours_ast_to_alp
FROM coh c
LEFT JOIN pt         USING (patientunitstayid)
LEFT JOIN pk      p  USING (patientunitstayid)
LEFT JOIN at_peak a  USING (patientunitstayid)
LEFT JOIN t_ast   ta USING (patientunitstayid)
LEFT JOIN t_bili  tb USING (patientunitstayid)
LEFT JOIN t_alp   tp USING (patientunitstayid)
ORDER BY c.patientunitstayid;
