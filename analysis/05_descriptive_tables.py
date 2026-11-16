# -*- coding: utf-8 -*-
"""Descriptive tables: cohort by etiology (Table 1), cohort by peak bilirubin
(Table 2), cross-cohort reproducibility (Table 5), observation density, and
the survivor-versus-decedent comparison within the exposed group.

All laboratory values are peaks over the whole hospital admission; see
docs/variables.md.
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, fisher_exact
import config as cfg  # not "C": patsy uses C() for categorical terms

der = (cfg.load('cohort_mimic')
       .merge(cfg.load('labs_mimic'), on='hadm_id', how='left')
       .merge(cfg.load('chol_mimic')[['hadm_id', 'v_peak_alp', 'v_peak_bili',
                                    'dbili_at_bili_peak']],
              on='hadm_id', how='left'))
der['ast_alt'] = der.ast_peak_all / der.alt_peak_all


def mw(a, b):
    return cfg.pfmt(mannwhitneyu(pd.Series(a).dropna(), pd.Series(b).dropna()).pvalue)


# ------------------------------------------------- Table 1: by etiology
ORD = [('sepsis', 'Septic'), ('cardiogenic', 'Cardiogenic'), ('hypovolemic', 'Hypovolemic')]
g = {k: der[der.shock_type == k] for k, _ in ORD}
SPEC = [('Age, median, y', lambda s: cfg.med(s.age)),
        ('SOFA, median', lambda s: cfg.med(s.sofa)),
        ('Peak ALT, median, U/L', lambda s: cfg.med(s.alt_peak_all)),
        ('Peak AST, median, U/L', lambda s: cfg.med(s.ast_peak_all)),
        ('Peak ALP, median, U/L', lambda s: cfg.med(s.v_peak_alp)),
        ('Peak LDH, median, U/L', lambda s: cfg.med(s.ldh_peak_all)),
        ('Peak creatinine, median, mg/dL', lambda s: cfg.med(s.peak_creat, 1)),
        ('Bilirubin ≥10, n (%)',
         lambda s: f'{int(s.bili_ge10.sum())} ({100 * s.bili_ge10.mean():.1f})'),
        ('In-hospital death, %', lambda s: f'{100 * s.hospital_death.mean():.1f}')]
t1 = pd.DataFrame([[lab] + [fn(g[k]) for k, _ in ORD] for lab, fn in SPEC],
                  columns=['Characteristic'] + [f'{n} (n={len(g[k])})' for k, n in ORD])
cfg.save(t1, 'table1_by_etiology')

# -------------------------------------------- Table 2: by peak bilirubin
hi, lo = der[der.bili_ge10 == 1], der[der.bili_ge10 == 0]
rows = [['Age, median [IQR], y', cfg.iqr(hi.age), cfg.iqr(lo.age), mw(hi.age, lo.age)],
        ['Male, n (%)',
         f'{int((hi.gender == "M").sum())} ({100 * (hi.gender == "M").mean():.1f})',
         f'{int((lo.gender == "M").sum())} ({100 * (lo.gender == "M").mean():.1f})',
         cfg.pfmt(fisher_exact(pd.crosstab(der.bili_ge10, der.gender).values)[1])]]
for lab, col, dec in [('SOFA, median [IQR]', 'sofa', 0),
                      ('Peak ALT, U/L', 'alt_peak_all', 0),
                      ('Peak AST, U/L', 'ast_peak_all', 0),
                      ('Peak ALP, U/L', 'v_peak_alp', 0),
                      ('AST/ALT ratio', 'ast_alt', 2),
                      ('Peak creatinine, mg/dL', 'peak_creat', 1),
                      ('Peak INR', 'peak_inr', 1),
                      ('Time to bilirubin peak, h', 'bili_peak_hours', 0)]:
    rows.append([lab, cfg.iqr(hi[col], dec), cfg.iqr(lo[col], dec), mw(hi[col], lo[col])])
rows.append(['In-hospital death, n (%)',
             f'{int(hi.hospital_death.sum())} ({100 * hi.hospital_death.mean():.1f})',
             f'{int(lo.hospital_death.sum())} ({100 * lo.hospital_death.mean():.1f})',
             '<0.001'])
cfg.save(pd.DataFrame(rows, columns=['Characteristic',
                                   f'Bilirubin ≥10 (n={len(hi)})',
                                   f'Bilirubin <10 (n={len(lo)})', 'p']),
       'table2_by_bilirubin')

# ----------------------------------- Table 5: cross-cohort reproducibility
eic = (cfg.load('cohort_eicu')
       .merge(cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score',
                                     'hospitalid']],
              on='patientunitstayid', how='left')
       .merge(cfg.load('chol_eicu')[['patientunitstayid', 'v_peak_alp', 'v_peak_bili',
                                   'dbili_at_bili_peak']],
              on='patientunitstayid', how='left'))


def or_str(m):
    o, ci, _ = cfg.or_ci(m)
    return f'{o} ({ci})'


def dratio(df, bcol='v_peak_bili'):
    r = df.dropna(subset=['dbili_at_bili_peak'])
    r = r[r[bcol] > 0]
    r = r.assign(ra=r.dbili_at_bili_peak / r[bcol])
    return r[(r.ra >= 0) & (r.ra <= 1.2)]


def timing(df, sh):
    s = df[df.shock_type == sh]
    return (f'{cfg.med(s[s.bili_ge10 == 1].bili_peak_hours)} vs '
            f'{cfg.med(s[s.bili_ge10 == 0].bili_peak_hours)} h')


# The two rows below must be the same model in both columns. The first is the
# model without the severity score (the derivation column previously carried
# the SOFA-adjusted estimate under the "no severity score" label, which made
# the same number appear in both rows); the second adds it.
md0 = smf.logit(cfg.F_DER_NO_SOFA, data=der).fit(disp=0)
mv = smf.logit('hospital_death ~ bili_ge10 + age + C(shock_type)', data=eic).fit(disp=0)
md = smf.logit(cfg.F_DER, data=der).fit(disp=0)
mva = smf.logit(cfg.F_VAL, data=eic).fit(disp=0)
rd_hi, rd_lo = dratio(der[der.bili_ge10 == 1]), dratio(der[der.bili_ge10 == 0])
re_hi, re_lo = dratio(eic[eic.bili_ge10 == 1]), dratio(eic[eic.bili_ge10 == 0])

t5 = pd.DataFrame([
    ['Cohort size, n', f'{len(der):,}', f'{len(eic):,}'],
    ['Contributing hospitals, n', '1', f'{eic.hospitalid.nunique()}'],
    ['Death, bilirubin ≥10 vs <10',
     f'{100 * hi.hospital_death.mean():.1f}% vs {100 * lo.hospital_death.mean():.1f}%',
     f'{100 * eic[eic.bili_ge10 == 1].hospital_death.mean():.1f}% vs '
     f'{100 * eic[eic.bili_ge10 == 0].hospital_death.mean():.1f}%'],
    ['Adjusted OR, bilirubin ≥10 (age and etiology)', or_str(md0), or_str(mv)],
    ['Adjusted OR with severity score added',
     f'{or_str(md)}, SOFA', f'{or_str(mva)}, APACHE IVa'],
    ['Bilirubin peak, septic (≥10 vs <10)', timing(der, 'sepsis'), timing(eic, 'sepsis')],
    ['Bilirubin peak, cardiogenic (≥10 vs <10)',
     timing(der, 'cardiogenic'), timing(eic, 'cardiogenic')],
    ['Peak ALP, ≥10 vs <10, U/L',
     f'{cfg.med(hi.v_peak_alp)} vs {cfg.med(lo.v_peak_alp)}',
     f'{cfg.med(eic[eic.bili_ge10 == 1].v_peak_alp)} vs '
     f'{cfg.med(eic[eic.bili_ge10 == 0].v_peak_alp)}'],
    ['Direct/total bilirubin ratio, ≥10 vs <10',
     f'{cfg.med(rd_hi.ra, 2)} vs {cfg.med(rd_lo.ra, 2)} (n={len(rd_hi) + len(rd_lo)})',
     f'{cfg.med(re_hi.ra, 2)} vs {cfg.med(re_lo.ra, 2)} (n={len(re_hi) + len(re_lo)})'],
], columns=['Feature', 'MIMIC-IV (derivation)', 'eICU-CRD (validation)'])
cfg.save(t5, 'table5_reproducibility')

# ------------------------------------------------- observation density
le = cfg.load('labs_eicu')
print(f'  ALT measurements, derivation  {cfg.iqr(der.n_alt)}')
print(f'  ALT measurements, validation  {cfg.iqr(le.n_alt)}')
print(f'  observation span, derivation  {cfg.iqr(der.obs_span_hr)} h')
print(f'  observation span, validation  {cfg.iqr(le.last_lab_hr - le.first_lab_hr)} h')

# ------------------------------- survivors within the high-bilirubin group
s, d = hi[hi.hospital_death == 0], hi[hi.hospital_death == 1]
print(f'\n  high-bilirubin group: {len(s)} survivors, {len(d)} decedents')
for lab, col, dec in [('SOFA', 'sofa', 0), ('peak AST', 'ast_peak_all', 0),
                      ('peak ALT', 'alt_peak_all', 0),
                      ('peak bilirubin', 'peak_bili', 1),
                      ('AST/ALT ratio', 'ast_alt', 2)]:
    print(f'    {lab:16s} {cfg.iqr(s[col], dec)} vs {cfg.iqr(d[col], dec)}  '
          f'p={mw(s[col], d[col])}')
