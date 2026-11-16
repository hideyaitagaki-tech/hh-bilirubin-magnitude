-- =====================================================================
-- APACHE IVa, its physiological components, competing diagnoses and the
-- contributing hospital (eICU-CRD)
--
-- NOTE: apachepatientresult and apacheapsvar use -1 as a missing-value
-- sentinel. It is converted to NULL here; analysis/config.py repeats the
-- conversion defensively on load.
--
-- apacheapsvar has no 'vasopressors' column. Ventilation is 'vent' and
-- 'intubated'; renal replacement is 'dialysis'.
--
-- Writes: {{DATASET}}.apache_eicu
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.apache_eicu` AS
WITH coh AS (
  SELECT patientunitstayid, shock_type, age, bili_ge10, hospital_death
  FROM `{{DATASET}}.cohort_eicu`
),
pt AS (
  SELECT p.patientunitstayid, p.gender, p.hospitalid,
         p.unitdischargeoffset, p.hospitaldischargeoffset,
         p.hospitaldischargestatus, p.unitdischargestatus
  FROM `physionet-data.eicu_crd.patient` p
  JOIN coh c USING (patientunitstayid)
),
apr AS (
  SELECT patientunitstayid,
    NULLIF(MAX(IF(apacheversion = 'IVa', apachescore, NULL)), -1)                AS apache_iva_score,
    NULLIF(MAX(IF(apacheversion = 'IVa', predictedhospitalmortality, NULL)), -1) AS apache_iva_predmort,
    NULLIF(MAX(apachescore), -1)                                                  AS apache_any_score
  FROM `physionet-data.eicu_crd.apachepatientresult`
  GROUP BY patientunitstayid
),
aps AS (
  SELECT patientunitstayid,
         NULLIF(bilirubin,  -1) AS aps_bilirubin,
         NULLIF(creatinine, -1) AS aps_creatinine,
         NULLIF(meanbp,     -1) AS aps_meanbp,
         NULLIF(ph,         -1) AS aps_ph,
         NULLIF(pao2,       -1) AS aps_pao2,
         intubated              AS aps_intubated,
         vent                   AS aps_vent,
         dialysis               AS aps_dialysis
  FROM `physionet-data.eicu_crd.apacheapsvar`
),
dx AS (
  SELECT d.patientunitstayid,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.diagnosisstring),
      r'cirrhosis|chronic liver|portal hypertension|esophageal varice'), 1, 0)) AS dx_chronic_liver,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.diagnosisstring),
      r'cholangitis|biliary obstruction|common bile duct|choledocholithiasis|biliary tract'), 1, 0)) AS dx_biliary_obstruction,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.diagnosisstring),
      r'hemolysis|hemolytic'), 1, 0)) AS dx_hemolysis,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.diagnosisstring),
      r'rhabdomyolysis'), 1, 0)) AS dx_rhabdomyolysis
  FROM `physionet-data.eicu_crd.diagnosis` d
  JOIN coh c USING (patientunitstayid)
  GROUP BY d.patientunitstayid
)
SELECT
  c.patientunitstayid, c.shock_type, p.gender, c.age, c.bili_ge10, c.hospital_death,
  p.hospitalid, p.unitdischargeoffset, p.hospitaldischargeoffset,
  p.hospitaldischargestatus, p.unitdischargestatus,
  r.apache_iva_score, r.apache_iva_predmort, r.apache_any_score,
  v.aps_bilirubin, v.aps_creatinine, v.aps_meanbp, v.aps_ph, v.aps_pao2,
  v.aps_intubated, v.aps_vent, v.aps_dialysis,
  IFNULL(x.dx_chronic_liver, 0)       AS dx_chronic_liver,
  IFNULL(x.dx_biliary_obstruction, 0) AS dx_biliary_obstruction,
  IFNULL(x.dx_hemolysis, 0)           AS dx_hemolysis,
  IFNULL(x.dx_rhabdomyolysis, 0)      AS dx_rhabdomyolysis
FROM coh c
LEFT JOIN pt  p USING (patientunitstayid)
LEFT JOIN apr r USING (patientunitstayid)
LEFT JOIN aps v USING (patientunitstayid)
LEFT JOIN dx  x USING (patientunitstayid)
ORDER BY c.patientunitstayid;
