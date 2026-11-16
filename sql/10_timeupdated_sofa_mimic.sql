-- =====================================================================
-- Time-updated organ dysfunction at each landmark (MIMIC-IV)
--
-- The landmark models adjust for the SOFA score at time zero. A referee
-- asked whether a bilirubin rise merely marks organ failure that has
-- developed in the meantime, which a baseline score cannot capture. This
-- query takes the SOFA score in force at each landmark from the
-- mimiciv_derived hourly table.
--
-- The liver subscore is reported separately, so that the landmark models
-- can adjust for NON-HEPATIC organ dysfunction. Adjusting for the total
-- score would partly adjust for the exposure itself, because bilirubin
-- is the liver subscore.
--
-- Writes: {{DATASET}}.sofa_landmark_mimic
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.sofa_landmark_mimic` AS
WITH coh AS (
  SELECT hadm_id FROM `{{DATASET}}.cohort_mimic`
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
stays AS (
  SELECT i.hadm_id, i.stay_id
  FROM `physionet-data.mimiciv_3_1_icu.icustays` i
  JOIN coh c USING (hadm_id)
),
s AS (
  SELECT st.hadm_id,
         DATETIME_DIFF(sf.endtime, z.t_zero, HOUR)      AS hr,
         sf.sofa_24hours,
         sf.liver_24hours,
         sf.sofa_24hours - sf.liver_24hours             AS sofa_nonhepatic
  FROM `physionet-data.mimiciv_3_1_derived.sofa` sf
  JOIN stays st USING (stay_id)
  JOIN t0    z  ON z.hadm_id = st.hadm_id
  WHERE sf.sofa_24hours IS NOT NULL
)
SELECT
  hadm_id,
  -- the last score in force at or before each landmark
  MAX(IF(hr BETWEEN 0 AND  48, sofa_24hours,    NULL)) AS sofa_max_48h,
  MAX(IF(hr BETWEEN 0 AND  72, sofa_24hours,    NULL)) AS sofa_max_72h,
  MAX(IF(hr BETWEEN 0 AND  96, sofa_24hours,    NULL)) AS sofa_max_96h,
  MAX(IF(hr BETWEEN 0 AND 168, sofa_24hours,    NULL)) AS sofa_max_168h,
  MAX(IF(hr BETWEEN 0 AND  48, sofa_nonhepatic, NULL)) AS sofa_nh_48h,
  MAX(IF(hr BETWEEN 0 AND  72, sofa_nonhepatic, NULL)) AS sofa_nh_72h,
  MAX(IF(hr BETWEEN 0 AND  96, sofa_nonhepatic, NULL)) AS sofa_nh_96h,
  MAX(IF(hr BETWEEN 0 AND 168, sofa_nonhepatic, NULL)) AS sofa_nh_168h,
  MAX(IF(hr BETWEEN 0 AND 168, liver_24hours,   NULL)) AS liver_max_168h,
  COUNTIF(hr >= 0)                                     AS n_sofa_after_t0
FROM s
GROUP BY hadm_id;
