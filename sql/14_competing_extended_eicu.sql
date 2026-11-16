-- =====================================================================
-- Extended competing causes and the time of first crossing (eICU-CRD)
--
-- The eICU analogue of 13_competing_extended_mimic.sql. Diagnostic
-- entry in eICU-CRD is optional and free text, so these flags are read
-- from diagnosis strings and are certain to under-ascertain; they are
-- reported as such.
--
-- Time zero reproduces the definition the cohort was built on:
-- vitalperiodic.systemicmean only. See docs/time_zero.md.
--
-- Writes: {{DATASET}}.competing_extended_eicu
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.competing_extended_eicu` AS
WITH coh AS (
  SELECT patientunitstayid, shock_type, age, hospital_death, bili_ge10
  FROM `{{DATASET}}.eicu_pure`
),
t0 AS (
  SELECT v.patientunitstayid, MIN(v.observationoffset) AS t_zero_min
  FROM `physionet-data.eicu_crd.vitalperiodic` v
  JOIN coh USING (patientunitstayid)
  WHERE v.systemicmean IS NOT NULL AND v.systemicmean > 0 AND v.systemicmean < 65
  GROUP BY v.patientunitstayid
),
stay AS (
  SELECT p.patientunitstayid, p.hospitaldischargeoffset, p.hospitaldischargestatus
  FROM `physionet-data.eicu_crd.patient` p
  JOIN coh c USING (patientunitstayid)
),
bili AS (
  SELECT l.patientunitstayid, l.labresult,
         (l.labresultoffset - z.t_zero_min) / 60.0 AS hr
  FROM `physionet-data.eicu_crd.lab` l
  JOIN t0 z USING (patientunitstayid)
  WHERE l.labname = 'total bilirubin'
    AND l.labresult IS NOT NULL AND l.labresult > 0
),
dx AS (
  SELECT d.patientunitstayid,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.diagnosisstring),
      r'hepatitis a|hepatitis b|hepatitis c|viral hepatitis|acute hepatitis'), 1, 0))
      AS dx_viral_hepatitis_acute,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.diagnosisstring),
      r'acetaminophen|tylenol|paracetamol'), 1, 0)) AS dx_acetaminophen,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.diagnosisstring),
      r'drug-induced|drug induced|toxic hepat|hepatotoxic'), 1, 0)) AS dx_toxic_liver,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.diagnosisstring),
      r'parenteral nutrition|tpn'), 1, 0)) AS dx_parenteral_nutrition
  FROM `physionet-data.eicu_crd.diagnosis` d
  JOIN coh c USING (patientunitstayid)
  GROUP BY d.patientunitstayid
)
SELECT
  c.patientunitstayid, c.shock_type, c.age, c.hospital_death, c.bili_ge10,
  IF(z.t_zero_min IS NULL, 0, 1)                                AS has_t0,
  IF(s.hospitaldischargestatus = 'Expired',
     (s.hospitaldischargeoffset - z.t_zero_min) / 60.0, NULL)   AS hours_t0_to_death,
  (s.hospitaldischargeoffset - z.t_zero_min) / 60.0             AS hours_t0_to_discharge,
  MIN(IF(b.labresult >= 10 AND b.hr >= 0, b.hr, NULL))          AS first_ge10_hr,
  IFNULL(MAX(x.dx_viral_hepatitis_acute), 0)                    AS dx_viral_hepatitis_acute,
  IFNULL(MAX(x.dx_acetaminophen), 0)                            AS dx_acetaminophen,
  IFNULL(MAX(x.dx_toxic_liver), 0)                              AS dx_toxic_liver,
  IFNULL(MAX(x.dx_parenteral_nutrition), 0)                     AS dx_parenteral_nutrition
FROM coh c
LEFT JOIN t0   z USING (patientunitstayid)
LEFT JOIN stay s USING (patientunitstayid)
LEFT JOIN bili b USING (patientunitstayid)
LEFT JOIN dx   x USING (patientunitstayid)
GROUP BY c.patientunitstayid, c.shock_type, c.age, c.hospital_death, c.bili_ge10,
         has_t0, hours_t0_to_death, hours_t0_to_discharge;

SELECT * FROM `{{DATASET}}.competing_extended_eicu` ORDER BY patientunitstayid;
