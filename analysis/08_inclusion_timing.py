# -*- coding: utf-8 -*-
"""When did the aminotransferase peak occur relative to shock onset?

The inclusion criterion for hypoxic hepatitis (peak ALT > 500 U/L) is applied
over the whole admission, so an admission could in principle qualify on a
value measured before time zero — that is, on liver injury that preceded the
shock rather than followed it. This script dates the qualifying value and the
ALT peak against time zero, and repeats the primary model under progressively
stricter requirements that the injury post-date the shock.

Writes eTable 11.
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import config as cfg  # not "C": patsy uses C() for categorical terms

der = cfg.load('cohort_mimic').merge(cfg.load('alt_timing_mimic'),
                                     on='hadm_id', how='left')
eic = (cfg.load('cohort_eicu')
       .merge(cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score']],
              on='patientunitstayid', how='left')
       .merge(cfg.load('alt_timing_eicu'), on='patientunitstayid', how='left'))
SPECS = [('Derivation', der, cfg.F_DER), ('Validation', eic, cfg.F_VAL)]

# ------------------------------------------------- how the timing falls out
rows = []
for name, d, _ in SPECS:
    q = d.hours_t0_to_first_alt_gt500
    only_pre = (d.alt_peak_post_t0.fillna(0) <= 500) & (d.alt_peak_pre_t0.fillna(0) > 500)
    rows.append([
        name, len(d),
        f'{int((q >= 0).sum())} ({100 * (q >= 0).mean():.1f})',
        f'{int((d.alt_peak_post_t0 > 500).sum())} '
        f'({100 * (d.alt_peak_post_t0 > 500).mean():.1f})',
        f'{int(only_pre.sum())} ({100 * only_pre.mean():.1f})',
        cfg.iqr(d.hours_t0_to_alt_peak, 0),
        f'{int((d.n_alt_pre_t0 > 0).sum())} ({100 * (d.n_alt_pre_t0 > 0).mean():.1f})'])
cfg.save(pd.DataFrame(rows, columns=[
    'Cohort', 'n', 'First ALT >500 U/L at or after time zero, n (%)',
    'Qualifies on a value at or after time zero, n (%)',
    'Qualifies only on a value before time zero, n (%)',
    'Time from time zero to the ALT peak, median [IQR], h',
    'Any ALT measured before time zero, n (%)']),
    'etable11a_inclusion_timing')

# --------------------------------------- the model under stricter inclusion
RESTRICT = [
    ('None (primary analysis)', lambda x: x.index == x.index),
    ('Qualifying ALT >500 U/L measured at or after time zero',
     lambda x: x.alt_peak_post_t0 > 500),
    ('ALT peak occurring at or after time zero',
     lambda x: x.hours_t0_to_alt_peak >= 0),
    # The strictest restriction, which the referee did not ask for: it keeps
    # only admissions whose first blood draw followed the shock, and so
    # selects a different population (community-onset shock presenting to the
    # emergency department). It is reported because it is the one restriction
    # under which the association loses significance, and the reason is a loss
    # of two thirds of the cohort rather than a change of direction.
    ('No ALT measured at all before time zero', lambda x: x.n_alt_pre_t0 == 0),
]
rows = []
for name, d, f in SPECS:
    for lab, fn in RESTRICT:
        s = d[fn(d)]
        m = smf.logit(f, data=s).fit(disp=0)
        o, ci, p = cfg.or_ci(m)
        rows.append([name, lab, len(s), int(s.bili_ge10.sum()),
                     f'{100 * s.hospital_death.mean():.1f}', o, ci, cfg.pfmt(p)])
cfg.save(pd.DataFrame(rows, columns=['Cohort', 'Inclusion requirement', 'n',
                                     'Exposed, n', 'Mortality (%)',
                                     'aOR', '95% CI', 'p']),
         'etable11b_inclusion_sensitivity')

# ------------- the admissions excluded by the restriction are a distinct group
for name, d, _ in SPECS:
    only_pre = (d.alt_peak_post_t0.fillna(0) <= 500) & (d.alt_peak_pre_t0.fillna(0) > 500)
    a, b = d[only_pre], d[~only_pre]
    print(f'  {name}: {int(only_pre.sum())} admissions qualify only on a pre-shock value; '
          f'mortality {100 * a.hospital_death.mean():.1f}% vs '
          f'{100 * b.hospital_death.mean():.1f}%, '
          f'ALT peak {cfg.iqr(a.hours_t0_to_alt_peak, 0)} vs '
          f'{cfg.iqr(b.hours_t0_to_alt_peak, 0)} h from time zero')
