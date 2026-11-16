-- =====================================================================
-- Landmark exposures, time of death, competing diagnoses (MIMIC-IV)
--
-- bili_max_<L> is the highest total bilirubin measured in the first L
-- hours after time zero. The landmark analysis uses only these columns,
-- so that no patient can enter the exposed group by surviving longer.
--
-- Recovered verbatim from the BigQuery job, with one correction: the job
-- as run derived time zero from three itemids with a 20 mmHg floor, which
-- is not the definition the cohort itself was built on. Time zero here
-- reproduces the cohort definition (two itemids, floor at zero), so that
-- one definition holds throughout. See docs/time_zero.md.
--
-- dx_acute_liver_failure is carried because the original job computed it;
-- it is not used in any analysis. See docs/diagnosis_codes.md, "A code set
-- that was removed".
--
-- Writes: {{DATASET}}.landmark_mimic
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.landmark_mimic` AS
WITH coh AS (
  SELECT hadm_id, shock_type, gender, age, sofa, bili_ge10, hospital_death
  FROM `{{DATASET}}.allshock_pure`
),
t0 AS (
  -- the cohort's own definition: itemids 220052 and 220181, floor at zero
  SELECT ce.hadm_id, MIN(ce.charttime) AS t_zero
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  JOIN coh c USING (hadm_id)
  WHERE ce.itemid IN (220052, 220181)
    AND ce.valuenum < 65
    AND ce.valuenum > 0
  GROUP BY ce.hadm_id
),
adm AS (
  SELECT a.hadm_id, a.dischtime, a.deathtime
  FROM `physionet-data.mimiciv_3_1_hosp.admissions` a
  JOIN coh c USING (hadm_id)
),
bili AS (
  SELECT l.hadm_id, l.valuenum,
         TIMESTAMP_DIFF(l.charttime, z.t_zero, HOUR) AS hr
  FROM `physionet-data.mimiciv_3_1_hosp.labevents` l
  JOIN t0 z USING (hadm_id)
  WHERE l.itemid = 50885 AND l.valuenum IS NOT NULL
),
dx AS (
  SELECT d.hadm_id,
    MAX(IF(d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(K703|K704|K717|K721|K74|K766|K767|B18|I85)')
        OR d.icd_version = 9  AND REGEXP_CONTAINS(d.icd_code, r'^(571|5723|4560|4561|4562|07032|07044|07054)'), 1, 0)) AS dx_chronic_liver,
    MAX(IF(d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(K80|K830|K831|K82|C24|C221)')
        OR d.icd_version = 9  AND REGEXP_CONTAINS(d.icd_code, r'^(5740|5743|5745|5748|5761|5762|1561)'), 1, 0)) AS dx_biliary_obstruction,
    MAX(IF(d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(D55|D56|D57|D58|D59)')
        OR d.icd_version = 9  AND REGEXP_CONTAINS(d.icd_code, r'^(282|283)'), 1, 0)) AS dx_hemolysis,
    MAX(IF(d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(K72)')
        OR d.icd_version = 9  AND REGEXP_CONTAINS(d.icd_code, r'^(5701|570)'), 1, 0)) AS dx_acute_liver_failure,
    MAX(IF(d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(M6282|T795)')
        OR d.icd_version = 9  AND REGEXP_CONTAINS(d.icd_code, r'^(72888|9584)'), 1, 0)) AS dx_rhabdomyolysis
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d
  JOIN coh c USING (hadm_id)
  GROUP BY d.hadm_id
)
SELECT
  c.hadm_id, c.shock_type, c.gender, c.age, c.sofa, c.bili_ge10, c.hospital_death,
  IF(z.t_zero IS NULL, 0, 1)                                   AS has_t0,
  TIMESTAMP_DIFF(a.deathtime, z.t_zero, HOUR)                  AS hours_t0_to_death,
  TIMESTAMP_DIFF(a.dischtime, z.t_zero, HOUR)                  AS hours_t0_to_discharge,
  MAX(IF(b.hr BETWEEN 0 AND  48, b.valuenum, NULL))            AS bili_max_48h,
  MAX(IF(b.hr BETWEEN 0 AND  72, b.valuenum, NULL))            AS bili_max_72h,
  MAX(IF(b.hr BETWEEN 0 AND  96, b.valuenum, NULL))            AS bili_max_96h,
  MAX(IF(b.hr BETWEEN 0 AND 168, b.valuenum, NULL))            AS bili_max_168h,
  MAX(IF(b.hr >= 0, b.valuenum, NULL))                         AS bili_max_all,
  COUNTIF(b.hr >= 0)                                           AS n_bili_after_t0,
  IFNULL(MAX(x.dx_chronic_liver), 0)                           AS dx_chronic_liver,
  IFNULL(MAX(x.dx_biliary_obstruction), 0)                     AS dx_biliary_obstruction,
  IFNULL(MAX(x.dx_hemolysis), 0)                               AS dx_hemolysis,
  IFNULL(MAX(x.dx_acute_liver_failure), 0)                     AS dx_acute_liver_failure,
  IFNULL(MAX(x.dx_rhabdomyolysis), 0)                          AS dx_rhabdomyolysis
FROM coh c
LEFT JOIN t0   z USING (hadm_id)
LEFT JOIN adm  a USING (hadm_id)
LEFT JOIN bili b USING (hadm_id)
LEFT JOIN dx   x USING (hadm_id)
GROUP BY c.hadm_id, c.shock_type, c.gender, c.age, c.sofa, c.bili_ge10,
         c.hospital_death, z.t_zero, a.deathtime, a.dischtime
ORDER BY c.hadm_id;
