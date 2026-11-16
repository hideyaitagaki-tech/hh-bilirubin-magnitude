# -*- coding: utf-8 -*-
"""Competing causes not covered by the original code sets, and a
time-varying model of the exposure (eTables 14 and 15).

(1) The four conditions already carried — chronic liver disease, biliary
    obstruction, hemolysis, rhabdomyolysis — leave out the two causes of
    ICU cholestasis that are best recognised, parenteral nutrition and
    drugs, and the two causes of acute hepatocellular injury that most
    resemble hypoxic hepatitis, acetaminophen toxicity and acute viral
    hepatitis. Each is excluded in turn and then all together.

(2) The landmark analyses hold observation constant at four fixed times,
    but a peak after 168 h is not covered by any of them. A Cox model in
    which the exposure switches on at the hour the bilirubin first
    reaches 10 mg/dL covers the whole admission: before that hour the
    patient contributes unexposed person-time, which is what the
    immortal-time objection asks for.

Writes eTables 14 and 15.
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from lifelines import CoxTimeVaryingFitter
import config as cfg  # not "C": patsy uses C() for categorical terms

der = (cfg.load('competing_mimic')
       .merge(cfg.load('cohort_mimic')[['hadm_id', 'sofa']].rename(
           columns={'sofa': 'sofa_c'}), on='hadm_id', how='left'))
val = (cfg.load('competing_eicu')
       .merge(cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score']],
              on='patientunitstayid', how='left'))

# ------------------------------------------------- eTable 14: exclusions
SETS = {
    'Derivation': [
        ('None (primary analysis)', None),
        ('Excluding acute viral hepatitis', 'dx_viral_hepatitis_acute'),
        ('Excluding acetaminophen poisoning', 'dx_acetaminophen'),
        ('Excluding toxic (drug-induced) liver disease', 'dx_toxic_liver'),
        ('Excluding parenteral nutrition', 'px_parenteral_nutrition'),
    ],
    'Validation': [
        ('None (primary analysis)', None),
        ('Excluding acute viral hepatitis', 'dx_viral_hepatitis_acute'),
        ('Excluding acetaminophen poisoning', 'dx_acetaminophen'),
        ('Excluding toxic (drug-induced) liver disease', 'dx_toxic_liver'),
        ('Excluding parenteral nutrition', 'dx_parenteral_nutrition'),
    ],
}
FORM = {'Derivation': 'hospital_death ~ bili_ge10 + age + sofa_c + C(shock_type)',
        'Validation': cfg.F_VAL}

rows, freq = [], []
for name, d in [('Derivation', der), ('Validation', val)]:
    cols = [c for _, c in SETS[name] if c]
    for lab, col in SETS[name]:
        s = d if col is None else d[d[col] == 0]
        m = smf.logit(FORM[name], data=s).fit(disp=0)
        o, ci, p = cfg.or_ci(m)
        rows.append([name, lab, len(s), int(s.bili_ge10.sum()), o, ci, cfg.pfmt(p)])
    s = d[(d[cols] == 0).all(axis=1)]
    m = smf.logit(FORM[name], data=s).fit(disp=0)
    o, ci, p = cfg.or_ci(m)
    rows.append([name, 'Excluding all four', len(s), int(s.bili_ge10.sum()),
                 o, ci, cfg.pfmt(p)])
    for col in cols:
        hi, lo = d[d.bili_ge10 == 1], d[d.bili_ge10 == 0]
        freq.append([name, col.replace('dx_', '').replace('px_', '').replace('_', ' '),
                     f'{100 * hi[col].mean():.1f}', f'{100 * lo[col].mean():.1f}',
                     f'{100 * d[col].mean():.1f}'])

cfg.save(pd.DataFrame(freq, columns=['Cohort', 'Condition', 'Bilirubin ≥10 (%)',
                                     'Bilirubin <10 (%)', 'Overall (%)']),
         'etable14a_competing_extended_freq')
cfg.save(pd.DataFrame(rows, columns=['Cohort', 'Exclusion', 'n', 'Exposed, n',
                                     'aOR', '95% CI', 'p']),
         'etable14b_competing_extended_sensitivity')

# ------------------------------------------------- eTable 15: time-varying Cox
# One row per patient per interval. A patient who never reaches 10 mg/dL
# contributes one unexposed interval; a patient who does contributes an
# unexposed interval up to first_ge10_hr and an exposed interval after it.
# Patients whose exposure begins at or after the end of follow-up, or whose
# follow-up has no length, cannot contribute an exposed interval and are
# carried as unexposed for the time they are observed.
def episodes(d, idc, covs):
    out = []
    for r in d.itertuples(index=False):
        stop = getattr(r, 'hours_t0_to_death')
        event = 1 if not pd.isna(stop) else 0
        if not event:
            stop = getattr(r, 'hours_t0_to_discharge')
        if pd.isna(stop) or stop <= 0:
            continue
        base = {'id': getattr(r, idc), **{c: getattr(r, c) for c in covs}}
        t = getattr(r, 'first_ge10_hr')
        if pd.isna(t) or t >= stop:
            out.append({**base, 'start': 0.0, 'stop': float(stop),
                        'exposed': 0, 'event': event})
        else:
            t = max(float(t), 0.0)
            if t > 0:
                out.append({**base, 'start': 0.0, 'stop': t,
                            'exposed': 0, 'event': 0})
                out.append({**base, 'start': t, 'stop': float(stop),
                            'exposed': 1, 'event': event})
            else:
                out.append({**base, 'start': 0.0, 'stop': float(stop),
                            'exposed': 1, 'event': event})
    return pd.DataFrame(out)


rows = []
for name, d, idc, covs in [
        ('Derivation', der, 'hadm_id', ['age', 'sofa_c', 'shock_type']),
        ('Validation', val, 'patientunitstayid',
         ['age', 'apache_iva_score', 'shock_type'])]:
    s = d.dropna(subset=[c for c in covs if c != 'shock_type']).copy()
    ep = episodes(s, idc, covs)
    ep = pd.get_dummies(ep, columns=['shock_type'], drop_first=True, dtype=float)
    ctv = CoxTimeVaryingFitter()
    ctv.fit(ep, id_col='id', event_col='event',
            start_col='start', stop_col='stop', show_progress=False)
    r = ctv.summary.loc['exposed']
    rows.append([name, 'Exposure as a time-varying covariate',
                 int(ep.id.nunique()), int(ep.event.sum()),
                 f'{np.exp(r["coef"]):.2f}',
                 f'{np.exp(r["coef lower 95%"]):.2f}–{np.exp(r["coef upper 95%"]):.2f}',
                 cfg.pfmt(r['p'])])
    print(f'{name}: {int(ep.id.nunique())} patients, {len(ep)} intervals, '
          f'{int(ep.event.sum())} deaths, HR {np.exp(r["coef"]):.2f}')

t15 = pd.DataFrame(rows, columns=['Cohort', 'Model', 'Patients', 'Deaths',
                                  'HR', '95% CI', 'p'])
print()
print(t15.to_string(index=False))
cfg.save(t15, 'etable15_time_varying_cox')
