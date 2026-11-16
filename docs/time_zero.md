# Time zero

Time zero is the first recorded mean arterial pressure (MAP) below 65 mmHg
during the hospital admission. Every interval reported in the manuscript that
is described as "from shock onset" is measured from this point.

## MIMIC-IV

Time zero for the cohort is the table `allshock_t0`, built as:

```sql
SELECT hadm_id, MIN(ce.charttime) AS t0_shock
FROM chartevents ce
WHERE ce.itemid IN (220052, 220181)
  AND ce.valuenum < 65 AND ce.valuenum > 0
GROUP BY hadm_id
```

| itemid | label | used in the cohort t0 |
|--------|-------|----------------------|
| 220052 | Arterial Blood Pressure mean | yes |
| 220181 | Non Invasive Blood Pressure mean | yes |
| 225312 | ART BP mean | **no** |

Every quantity carried in the cohort table — `bili_peak_hours`, the post-shock
enzyme, creatinine, INR and transfusion variables — is measured from this
`t0_shock`.

### A divergence, recorded rather than smoothed over

Two later extractions did not reproduce this definition exactly.

| Extraction | itemids | floor |
|------------|---------|-------|
| `allshock_t0` (cohort, canonical) | 220052, 220181 | `> 0` |
| the landmark job (`landmark_mimic`) | 220052, 220181, **225312** | **`> 20`** |
| `sofa_landmark_mimic`, `alt_timing_mimic` | 220052, 220181, **225312** | `> 0` |

Adding 225312 can only move time zero earlier or leave it unchanged; the 20
mmHg floor can only move it later or leave it unchanged. Comparing the second
and third definitions directly, time zero differs in 21 of 1,352 admissions
(1.6%), by a median of 0 hours and by up to 114 hours. The divergence from the
canonical `allshock_t0` has not yet been quantified; `sql/README.md` records
what remains to be checked.

## eICU-CRD

Time zero for the cohort is built from `vitalperiodic` alone:

```sql
SELECT patientunitstayid, MIN(observationoffset) AS t0_offset
FROM vitalperiodic
WHERE systemicmean IS NOT NULL AND systemicmean > 0 AND systemicmean < 65
GROUP BY patientunitstayid
```

`vitalaperiodic.noninvasivemean` is **not** used. A patient with no qualifying
`systemicmean` has no time zero and is not in the cohort, because `eicu_t0`
joins rather than left-joins the result. The later extractions
`landmark_eicu` and `alt_timing_eicu` additionally read `vitalaperiodic`,
which can only move time zero earlier; this divergence is also recorded in
`sql/README.md`.

## Known limitation

The two databases sample blood pressure differently: MIMIC-IV records ICU
monitor data, eICU-CRD records periodic vital-sign entries. Time zero is
therefore not identified with the same granularity in the two cohorts, and the
manuscript states that this may contribute to the earlier peak times observed
in validation.

For this reason the manuscript also reports timing **anchored to the
aspartate aminotransferase peak**, which does not depend on blood-pressure
sampling at all, and which allows direct comparison with Van den Broecke et al.
(Ann Intensive Care 2018;8:15). That anchoring is computed in
`sql/03_cholestatic_markers_mimic.sql` and `sql/04_cholestatic_markers_eicu.sql`
as `hours_ast_to_bili` and `hours_ast_to_alp`.
