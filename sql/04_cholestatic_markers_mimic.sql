-- =====================================================================
-- Cholestatic markers and timing anchored to the AST peak (MIMIC-IV)
--
-- Produces the evidence for the cholestatic interpretation: peak ALP,
-- GGT, direct bilirubin, the values at the bilirubin peak needed for the
-- R ratio, and the interval from the AST peak to the bilirubin and ALP
-- peaks. The last of these is the anchoring used by Van den Broecke et
-- al. and is independent of blood-pressure sampling.
--
-- Writes: {{DATASET}}.chol_mimic
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.chol_mimic` AS
WITH coh AS (
  SELECT hadm_id, shock_type, gender, age, sofa, bili_ge10, hospital_death
  FROM `{{DATASET}}.cohort_mimic`
),
lab AS (
  SELECT l.hadm_id, l.itemid, l.valuenum, l.charttime
  FROM `physionet-data.mimiciv_3_1_hosp.labevents` l
  JOIN coh c USING (hadm_id)
  WHERE l.itemid IN (50861,  -- ALT
                     50863,  -- ALP
                     50878,  -- AST
                     50883,  -- bilirubin, direct
                     50885,  -- bilirubin, total
                     50927)  -- GGT
    AND l.valuenum IS NOT NULL
    AND l.valuenum > 0
),
pk AS (
  SELECT hadm_id,
    MAX(IF(itemid = 50885, valuenum, NULL)) AS v_peak_bili,
    MAX(IF(itemid = 50878, valuenum, NULL)) AS v_peak_ast,
    MAX(IF(itemid = 50863, valuenum, NULL)) AS v_peak_alp,
    MAX(IF(itemid = 50927, valuenum, NULL)) AS peak_ggt,
    MAX(IF(itemid = 50883, valuenum, NULL)) AS peak_dbili,
    COUNTIF(itemid = 50863) AS n_alp,
    COUNTIF(itemid = 50885) AS n_bili,
    COUNTIF(itemid = 50883) AS n_dbili,
    COUNTIF(itemid = 50927) AS n_ggt
  FROM lab GROUP BY hadm_id
),
-- time of the AST peak, and of the bilirubin and ALP peaks
t_ast AS (
  SELECT hadm_id, ARRAY_AGG(charttime ORDER BY valuenum DESC, charttime LIMIT 1)[OFFSET(0)] AS t
  FROM lab WHERE itemid = 50878 GROUP BY hadm_id
),
t_bili AS (
  SELECT hadm_id, ARRAY_AGG(charttime ORDER BY valuenum DESC, charttime LIMIT 1)[OFFSET(0)] AS t
  FROM lab WHERE itemid = 50885 GROUP BY hadm_id
),
t_alp AS (
  SELECT hadm_id, ARRAY_AGG(charttime ORDER BY valuenum DESC, charttime LIMIT 1)[OFFSET(0)] AS t
  FROM lab WHERE itemid = 50863 GROUP BY hadm_id
),
-- values measured in the same draw as the bilirubin peak, for the R ratio
at_peak AS (
  SELECT b.hadm_id,
    MAX(IF(l.itemid = 50863, l.valuenum, NULL)) AS alp_at_bili_peak,
    MAX(IF(l.itemid = 50861, l.valuenum, NULL)) AS alt_at_bili_peak,
    MAX(IF(l.itemid = 50883, l.valuenum, NULL)) AS dbili_at_bili_peak
  FROM t_bili b
  JOIN lab l ON l.hadm_id = b.hadm_id AND l.charttime = b.t
  GROUP BY b.hadm_id
)
SELECT
  c.*,
  p.v_peak_bili, p.v_peak_ast, p.v_peak_alp, p.peak_ggt, p.peak_dbili,
  p.n_alp, p.n_bili, p.n_dbili, p.n_ggt,
  a.alp_at_bili_peak, a.alt_at_bili_peak, a.dbili_at_bili_peak,
  DATETIME_DIFF(tb.t, ta.t, HOUR) AS hours_ast_to_bili,
  DATETIME_DIFF(tp.t, ta.t, HOUR) AS hours_ast_to_alp
FROM coh c
LEFT JOIN pk      p  USING (hadm_id)
LEFT JOIN at_peak a  USING (hadm_id)
LEFT JOIN t_ast   ta USING (hadm_id)
LEFT JOIN t_bili  tb USING (hadm_id)
LEFT JOIN t_alp   tp USING (hadm_id)
ORDER BY c.hadm_id;
