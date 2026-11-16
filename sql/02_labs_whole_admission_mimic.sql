-- =====================================================================
-- Whole-admission laboratory peaks and observation density (MIMIC-IV)
--
-- Peaks are taken over the entire hospital admission. An earlier version
-- of this analysis used a restricted window, which understated every
-- reported peak; see docs/exclusions.md, "Corrections".
--
-- The windowed columns (7 d, 14 d) are retained only so that the
-- discrepancy can be reproduced and audited. The manuscript uses the
-- *_all columns.
--
-- Writes: {{DATASET}}.labs_mimic
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.labs_mimic` AS
WITH coh AS (
  SELECT DISTINCT hadm_id
  FROM `{{DATASET}}.cohort_mimic`
),
lab AS (
  SELECT l.hadm_id, l.itemid, l.valuenum, l.charttime
  FROM `physionet-data.mimiciv_3_1_hosp.labevents` l
  JOIN coh c USING (hadm_id)
  WHERE l.itemid IN (50861,   -- ALT
                     50878,   -- AST
                     50954)   -- LDH
    AND l.valuenum IS NOT NULL
    AND l.valuenum > 0
),
-- Anchor relative times to the first laboratory result of the admission
-- rather than to admittime, so that this script needs only one table.
anchor AS (
  SELECT hadm_id, MIN(charttime) AS t0
  FROM lab
  GROUP BY hadm_id
),
rel AS (
  SELECT l.hadm_id, l.itemid, l.valuenum,
         DATETIME_DIFF(l.charttime, a.t0, HOUR) AS hr
  FROM lab l
  JOIN anchor a USING (hadm_id)
)
SELECT
  c.hadm_id,
  MAX(IF(r.itemid = 50861, r.valuenum, NULL))                         AS alt_peak_all,
  MAX(IF(r.itemid = 50861 AND r.hr <= 168, r.valuenum, NULL))         AS alt_peak_7d,
  MAX(IF(r.itemid = 50861 AND r.hr <= 336, r.valuenum, NULL))         AS alt_peak_14d,
  MAX(IF(r.itemid = 50878, r.valuenum, NULL))                         AS ast_peak_all,
  MAX(IF(r.itemid = 50878 AND r.hr <= 168, r.valuenum, NULL))         AS ast_peak_7d,
  MAX(IF(r.itemid = 50954, r.valuenum, NULL))                         AS ldh_peak_all,
  MAX(r.hr)                                                            AS obs_span_hr,
  COUNTIF(r.itemid = 50861)                                            AS n_alt
FROM coh c
LEFT JOIN rel r USING (hadm_id)
GROUP BY c.hadm_id
ORDER BY c.hadm_id;
