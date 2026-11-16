-- =====================================================================
-- Query D: the same, for eICU-CRD.
-- hr is hours from time zero (first MAP < 65 mmHg); negative = before.
-- =====================================================================
WITH coh AS (
  SELECT patientunitstayid
  FROM `my-mimic-study-507103.cholecystitis_study.eicu_pure`
),
t0 AS (
  -- verbatim from the job that built `eicu_t0`: vitalperiodic only
  SELECT v.patientunitstayid, MIN(v.observationoffset) AS t_zero_min
  FROM `physionet-data.eicu_crd.vitalperiodic` v
  JOIN coh USING (patientunitstayid)
  WHERE v.systemicmean IS NOT NULL AND v.systemicmean > 0 AND v.systemicmean < 65
  GROUP BY v.patientunitstayid
),
lab AS (
  SELECT l.patientunitstayid,
         l.labname,
         l.labresult,
         (l.labresultoffset - z.t_zero_min) / 60.0 AS hr
  FROM `physionet-data.eicu_crd.lab` l
  JOIN t0 z USING (patientunitstayid)
  WHERE l.labname IN ('ALT (SGPT)', 'AST (SGOT)')
    AND l.labresult IS NOT NULL AND l.labresult > 0
),
alt_pk AS (
  SELECT patientunitstayid, labresult AS alt_peak_all, hr AS hours_t0_to_alt_peak
  FROM (SELECT *, ROW_NUMBER() OVER (PARTITION BY patientunitstayid
                                     ORDER BY labresult DESC, hr ASC) AS rn
        FROM lab WHERE labname = 'ALT (SGPT)')
  WHERE rn = 1
),
ast_pk AS (
  SELECT patientunitstayid, labresult AS ast_peak_all, hr AS hours_t0_to_ast_peak
  FROM (SELECT *, ROW_NUMBER() OVER (PARTITION BY patientunitstayid
                                     ORDER BY labresult DESC, hr ASC) AS rn
        FROM lab WHERE labname = 'AST (SGOT)')
  WHERE rn = 1
),
agg AS (
  SELECT patientunitstayid,
    MAX(IF(labname = 'ALT (SGPT)' AND hr <  0, labresult, NULL)) AS alt_peak_pre_t0,
    MAX(IF(labname = 'ALT (SGPT)' AND hr >= 0, labresult, NULL)) AS alt_peak_post_t0,
    MAX(IF(labname = 'AST (SGOT)' AND hr <  0, labresult, NULL)) AS ast_peak_pre_t0,
    MAX(IF(labname = 'AST (SGOT)' AND hr >= 0, labresult, NULL)) AS ast_peak_post_t0,
    MIN(IF(labname = 'ALT (SGPT)' AND labresult > 500, hr, NULL)) AS hours_t0_to_first_alt_gt500,
    COUNTIF(labname = 'ALT (SGPT)' AND hr <  0)                   AS n_alt_pre_t0,
    COUNTIF(labname = 'ALT (SGPT)' AND hr >= 0)                   AS n_alt_post_t0,
    MIN(IF(labname = 'ALT (SGPT)', hr, NULL))                     AS first_alt_hr
  FROM lab
  GROUP BY patientunitstayid
)
SELECT
  c.patientunitstayid,
  IF(z.t_zero_min IS NULL, 0, 1) AS has_t0,
  a.alt_peak_all,  a.hours_t0_to_alt_peak,
  b.ast_peak_all,  b.hours_t0_to_ast_peak,
  g.alt_peak_pre_t0, g.alt_peak_post_t0,
  g.ast_peak_pre_t0, g.ast_peak_post_t0,
  g.hours_t0_to_first_alt_gt500,
  g.n_alt_pre_t0, g.n_alt_post_t0, g.first_alt_hr
FROM coh c
LEFT JOIN t0     z USING (patientunitstayid)
LEFT JOIN alt_pk a USING (patientunitstayid)
LEFT JOIN ast_pk b USING (patientunitstayid)
LEFT JOIN agg    g USING (patientunitstayid)
ORDER BY c.patientunitstayid;
