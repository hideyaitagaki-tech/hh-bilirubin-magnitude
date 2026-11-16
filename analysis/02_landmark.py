# -*- coding: utf-8 -*-
"""Landmark analysis (Table 4).

An exposure defined over the whole admission cannot be reached by a patient
who dies early, so the unexposed group collects the early deaths and the
association is diluted. At each landmark the analysis is restricted to
patients alive at that time and the exposure is redefined using only values
measured up to it, which removes that route.

The estimate rises rather than falls, which is the opposite of what the bias
would produce, and is why the manuscript describes the primary estimate as
conservative.
"""
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import statsmodels.formula.api as smf
import config as cfg  # not "C": patsy uses C() for categorical terms

LANDMARKS = [(48, 'bili_max_48h'), (72, 'bili_max_72h'),
             (96, 'bili_max_96h'), (168, 'bili_max_168h')]

der = cfg.load('cohort_mimic')
lm = cfg.load('landmark_mimic')

# The validation cohort carries the same design. Its landmark table was missing
# from the manuscript although the data for it were exported, which is why the
# only validation landmark figures quoted were the incident ones; both cohorts
# are now built here so that Table 4 says the same thing twice.
val = (cfg.load('cohort_eicu')
       .merge(cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score']],
              on='patientunitstayid', how='left'))
lme = (cfg.load('landmark_eicu')
       .merge(cfg.load('cohort_eicu')[['patientunitstayid', 'shock_type', 'age',
                                       'hospital_death']],
              on='patientunitstayid', how='left')
       .merge(cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score']],
              on='patientunitstayid', how='left'))

COHORTS = [
    ('Derivation (MIMIC-IV)', der, lm, cfg.F_DER,
     'hospital_death ~ exp + age + sofa + C(shock_type)'),
    ('Validation (eICU-CRD)', val, lme, cfg.F_VAL,
     'hospital_death ~ exp + age + C(shock_type) + apache_iva_score'),
]

rows = []
for cohort, base, land, f_primary, f_land in COHORTS:
    m = smf.logit(f_primary, data=base).fit(disp=0)
    o, ci, p = cfg.or_ci(m)
    rows.append([cohort, 'Primary analysis (no landmark)', len(base),
                 int(base.bili_ge10.sum()),
                 f'{100 * base.hospital_death.mean():.1f}', o, ci, cfg.pfmt(p)])
    for L, col in LANDMARKS:
        s = land[(land.has_t0 == 1) &
                 ((land.hours_t0_to_death.isna()) |
                  (land.hours_t0_to_death > L))].copy()
        s['exp'] = (s[col] >= 10).astype(int)
        mm = smf.logit(f_land, data=s).fit(disp=0)
        o, ci, p = cfg.or_ci(mm, 'exp')
        rows.append([cohort, f'Alive at {L} h', len(s), int(s.exp.sum()),
                     f'{100 * s.hospital_death.mean():.1f}', o, ci, cfg.pfmt(p)])

t4 = pd.DataFrame(rows, columns=['Cohort', 'Landmark', 'Patients at risk',
                                 'Exposed (n)', 'Mortality (%)', 'aOR',
                                 '95% CI', 'p'])
print(t4.to_string(index=False))
cfg.save(t4, 'table4_landmark')

# Where the dilution comes from: early deaths are unexposed by construction.
early = lm[(lm.has_t0 == 1) & (lm.hours_t0_to_death <= 48)]
unexp = lm[(lm.has_t0 == 1) & (lm.bili_ge10 == 0)]
print()
print(f'  died within 48 h:                 {len(early)}')
print(f'    of these, classified unexposed: {int((early.bili_ge10 == 0).sum())}')
print(f'  unexposed mortality, all:         '
      f'{100 * unexp.hospital_death.mean():.1f}%')
print(f'  unexposed mortality, excluding early deaths: '
      f'{100 * unexp[~(unexp.hours_t0_to_death <= 48)].hospital_death.mean():.1f}%')
