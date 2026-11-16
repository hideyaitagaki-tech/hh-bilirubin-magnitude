-- =====================================================================
-- Query C: when did the ALT peak occur relative to shock onset? (MIMIC-IV)
--
-- Referee point B-5: the inclusion criterion (ALT > 500 U/L) is applied to
-- the whole admission, so an admission could qualify on a value measured
-- BEFORE time zero. This query dates every aminotransferase peak against
-- time zero so that the cohort can be restricted to admissions whose
-- qualifying ALT was measured at or after shock onset.
--
-- hr is hours from time zero: negative = before shock onset.
-- =====================================================================
WITH coh AS (
  SELECT hadm_id
  FROM `my-mimic-study-507103.cholecystitis_study.allshock_pure`
),
t0 AS (
  -- verbatim from the job that built `allshock_t0`: two itemids, floor at zero
  SELECT ce.hadm_id, MIN(ce.charttime) AS t_zero
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  JOIN coh USING (hadm_id)
  WHERE ce.itemid IN (220052, 220181)
    AND ce.valuenum < 65
    AND ce.valuenum > 0
  GROUP BY ce.hadm_id
),
lab AS (
  SELECT l.hadm_id,
         l.itemid,
         l.valuenum,
         DATETIME_DIFF(l.charttime, z.t_zero, HOUR) AS hr
  FROM `physionet-data.mimiciv_3_1_hosp.labevents` l
  JOIN t0 z USING (hadm_id)
  WHERE l.itemid IN (50861, 50878)      -- 50861 ALT, 50878 AST
    AND l.valuenum IS NOT NULL AND l.valuenum > 0
),
alt_pk AS (
  SELECT hadm_id, valuenum AS alt_peak_all, hr AS hours_t0_to_alt_peak
  FROM (SELECT *, ROW_NUMBER() OVER (PARTITION BY hadm_id
                                     ORDER BY valuenum DESC, hr ASC) AS rn
        FROM lab WHERE itemid = 50861)
  WHERE rn = 1
),
ast_pk AS (
  SELECT hadm_id, valuenum AS ast_peak_all, hr AS hours_t0_to_ast_peak
  FROM (SELECT *, ROW_NUMBER() OVER (PARTITION BY hadm_id
                                     ORDER BY valuenum DESC, hr ASC) AS rn
        FROM lab WHERE itemid = 50878)
  WHERE rn = 1
),
agg AS (
  SELECT hadm_id,
    MAX(IF(itemid = 50861 AND hr <  0, valuenum, NULL)) AS alt_peak_pre_t0,
    MAX(IF(itemid = 50861 AND hr >= 0, valuenum, NULL)) AS alt_peak_post_t0,
    MAX(IF(itemid = 50878 AND hr <  0, valuenum, NULL)) AS ast_peak_pre_t0,
    MAX(IF(itemid = 50878 AND hr >= 0, valuenum, NULL)) AS ast_peak_post_t0,
    MIN(IF(itemid = 50861 AND valuenum > 500, hr, NULL)) AS hours_t0_to_first_alt_gt500,
    COUNTIF(itemid = 50861 AND hr <  0)                  AS n_alt_pre_t0,
    COUNTIF(itemid = 50861 AND hr >= 0)                  AS n_alt_post_t0,
    MIN(IF(itemid = 50861, hr, NULL))                    AS first_alt_hr
  FROM lab
  GROUP BY hadm_id
)
SELECT
  c.hadm_id,
  IF(z.t_zero IS NULL, 0, 1) AS has_t0,
  a.alt_peak_all,  a.hours_t0_to_alt_peak,
  b.ast_peak_all,  b.hours_t0_to_ast_peak,
  g.alt_peak_pre_t0, g.alt_peak_post_t0,
  g.ast_peak_pre_t0, g.ast_peak_post_t0,
  g.hours_t0_to_first_alt_gt500,
  g.n_alt_pre_t0, g.n_alt_post_t0, g.first_alt_hr
FROM coh c
LEFT JOIN t0     z USING (hadm_id)
LEFT JOIN alt_pk a USING (hadm_id)
LEFT JOIN ast_pk b USING (hadm_id)
LEFT JOIN agg    g USING (hadm_id)
ORDER BY c.hadm_id;
