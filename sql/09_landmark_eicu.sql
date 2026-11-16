-- =====================================================================
-- Landmark exposures and time of death (eICU-CRD)
--
-- The eICU analogue of 06_landmark_competing_mimic.sql. bili_max_<L> is
-- the highest total bilirubin measured in the first L hours after time
-- zero, so that no patient can enter the exposed group by surviving
-- longer, and bili_max_48h allows prevalent and incident
-- hyperbilirubinemia to be separated at a later landmark.
--
-- Time zero reproduces the definition used to build the cohort table
-- `eicu_t0`, so that it is identical to the one every other variable is
-- measured from.
--
-- Writes: {{DATASET}}.landmark_eicu
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.landmark_eicu` AS
WITH coh AS (
  SELECT patientunitstayid, shock_type, age, hospital_death, bili_ge10
  FROM `{{DATASET}}.eicu_pure`
),
t0 AS (
  -- verbatim from the job that built `eicu_t0`: vitalperiodic only
  SELECT v.patientunitstayid, MIN(v.observationoffset) AS t_zero_min
  FROM `physionet-data.eicu_crd.vitalperiodic` v
  JOIN coh USING (patientunitstayid)
  WHERE v.systemicmean IS NOT NULL AND v.systemicmean > 0 AND v.systemicmean < 65
  GROUP BY v.patientunitstayid
),
stay AS (
  SELECT p.patientunitstayid,
         p.unitdischargeoffset,
         p.hospitaldischargeoffset,
         p.hospitaldischargestatus
  FROM `physionet-data.eicu_crd.patient` p
  JOIN coh c USING (patientunitstayid)
),
bili AS (
  SELECT l.patientunitstayid,
         l.labresult,
         (l.labresultoffset - z.t_zero_min) / 60.0 AS hr
  FROM `physionet-data.eicu_crd.lab` l
  JOIN t0 z USING (patientunitstayid)
  WHERE l.labname = 'total bilirubin'
    AND l.labresult IS NOT NULL
    AND l.labresult > 0
)
SELECT
  c.*,
  IF(z.t_zero_min IS NULL, 0, 1)                                AS has_t0,
  -- Hours from time zero to hospital death; NULL for survivors, so the
  -- analysis code can select the patients alive at each landmark.
  IF(s.hospitaldischargestatus = 'Expired',
     (s.hospitaldischargeoffset - z.t_zero_min) / 60.0, NULL)   AS hours_t0_to_death,
  (s.hospitaldischargeoffset - z.t_zero_min) / 60.0             AS hours_t0_to_discharge,
  MAX(IF(b.hr BETWEEN 0 AND  48, b.labresult, NULL))            AS bili_max_48h,
  MAX(IF(b.hr BETWEEN 0 AND  72, b.labresult, NULL))            AS bili_max_72h,
  MAX(IF(b.hr BETWEEN 0 AND  96, b.labresult, NULL))            AS bili_max_96h,
  MAX(IF(b.hr BETWEEN 0 AND 168, b.labresult, NULL))            AS bili_max_168h,
  MAX(IF(b.hr >= 0, b.labresult, NULL))                         AS bili_max_all,
  COUNTIF(b.hr >= 0)                                            AS n_bili_after_t0
FROM coh c
LEFT JOIN t0   z USING (patientunitstayid)
LEFT JOIN stay s USING (patientunitstayid)
LEFT JOIN bili b USING (patientunitstayid)
GROUP BY c.patientunitstayid, c.shock_type, c.age, c.hospital_death, c.bili_ge10,
         z.t_zero_min, s.hospitaldischargestatus, s.hospitaldischargeoffset;
