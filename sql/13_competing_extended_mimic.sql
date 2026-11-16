-- =====================================================================
-- Extended competing causes and the time of first crossing (MIMIC-IV)
--
-- Two things the manuscript was asked for and did not have.
--
-- (1) Competing causes beyond the four already carried. Acute viral
--     hepatitis, acetaminophen poisoning and toxic liver disease, and
--     parenteral nutrition are the causes of a delayed cholestatic
--     picture in the ICU that the original code sets did not cover.
--     Chronic viral hepatitis (B18) is already inside dx_chronic_liver
--     and is deliberately not repeated here.
--
-- (2) first_ge10_hr, the hour after time zero at which total bilirubin
--     first reached 10 mg/dL. The landmark analyses need only the peak
--     within a window; a time-varying model needs the moment the
--     exposure begins, which no exported table carries.
--
-- Time zero is the cohort's own definition (itemids 220052 and 220181,
-- floor at zero); see docs/time_zero.md.
--
-- Writes: {{DATASET}}.competing_extended_mimic
-- =====================================================================
CREATE OR REPLACE TABLE `{{DATASET}}.competing_extended_mimic` AS
WITH coh AS (
  SELECT hadm_id, shock_type, gender, age, sofa, bili_ge10, hospital_death
  FROM `{{DATASET}}.allshock_pure`
),
t0 AS (
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
    -- acute viral hepatitis; chronic (B18 / 07032, 07044, 07054) is
    -- already counted as chronic liver disease and is excluded here
    MAX(IF(d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(B15|B16|B17|B19)')
        OR d.icd_version = 9  AND REGEXP_CONTAINS(d.icd_code, r'^(0700|0701|0702|0703|0704[0-3]|0705[0-3]|0706|0709)'),
        1, 0)) AS dx_viral_hepatitis_acute,
    -- poisoning by 4-aminophenol derivatives (acetaminophen)
    MAX(IF(d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(T391)')
        OR d.icd_version = 9  AND REGEXP_CONTAINS(d.icd_code, r'^(9654|E9354)'), 1, 0)) AS dx_acetaminophen,
    -- toxic liver disease, which is how drug-induced injury is coded
    MAX(IF(d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(K71)')
        OR d.icd_version = 9  AND REGEXP_CONTAINS(d.icd_code, r'^(5736)'), 1, 0)) AS dx_toxic_liver
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d
  JOIN coh c USING (hadm_id)
  GROUP BY d.hadm_id
),
px AS (
  -- parenteral nutrition: ICD-10-PCS introduction of a nutritional
  -- substance into a central vein, ICD-9 parenteral infusion of
  -- concentrated nutritional substances
  SELECT p.hadm_id,
    MAX(IF(p.icd_version = 10 AND REGEXP_CONTAINS(p.icd_code, r'^(3E04)')
        OR p.icd_version = 9  AND REGEXP_CONTAINS(p.icd_code, r'^(9915)'), 1, 0)) AS px_parenteral_nutrition
  FROM `physionet-data.mimiciv_3_1_hosp.procedures_icd` p
  JOIN coh c USING (hadm_id)
  GROUP BY p.hadm_id
)
SELECT
  c.hadm_id, c.shock_type, c.gender, c.age, c.sofa, c.bili_ge10, c.hospital_death,
  IF(z.t_zero IS NULL, 0, 1)                                   AS has_t0,
  TIMESTAMP_DIFF(a.deathtime, z.t_zero, HOUR)                  AS hours_t0_to_death,
  TIMESTAMP_DIFF(a.dischtime, z.t_zero, HOUR)                  AS hours_t0_to_discharge,
  MIN(IF(b.valuenum >= 10 AND b.hr >= 0, b.hr, NULL))          AS first_ge10_hr,
  IFNULL(MAX(x.dx_viral_hepatitis_acute), 0)                   AS dx_viral_hepatitis_acute,
  IFNULL(MAX(x.dx_acetaminophen), 0)                           AS dx_acetaminophen,
  IFNULL(MAX(x.dx_toxic_liver), 0)                             AS dx_toxic_liver,
  IFNULL(MAX(q.px_parenteral_nutrition), 0)                    AS px_parenteral_nutrition
FROM coh c
LEFT JOIN t0   z USING (hadm_id)
LEFT JOIN adm  a USING (hadm_id)
LEFT JOIN bili b USING (hadm_id)
LEFT JOIN dx   x USING (hadm_id)
LEFT JOIN px   q USING (hadm_id)
GROUP BY c.hadm_id, c.shock_type, c.gender, c.age, c.sofa, c.bili_ge10,
         c.hospital_death, has_t0, hours_t0_to_death, hours_t0_to_discharge;

SELECT * FROM `{{DATASET}}.competing_extended_mimic` ORDER BY hadm_id;
