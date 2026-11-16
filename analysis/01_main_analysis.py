# -*- coding: utf-8 -*-
"""Primary analysis: the association between a marked peak bilirubin and
in-hospital mortality, in both cohorts, with the severity adjustment that
each database supports.

Produces the numbers reported in Results, "Bilirubin >=10 mg/dL and mortality"
and "Reproduction of the principal finding and adjustment for admission
severity", and Table 3.
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu, spearmanr
import config as cfg  # not "C": patsy uses C() for categorical terms

print('=' * 70)
print('DERIVATION COHORT (MIMIC-IV)')
print('=' * 70)

der = cfg.load('cohort_mimic')
m = smf.logit(cfg.F_DER, data=der).fit(disp=0)
o, ci, p = cfg.or_ci(m)
print(f'  primary                     aOR {o} ({ci})  p={cfg.pfmt(p)}  n={int(m.nobs)}')

# etiology-specific models
for sh in ('sepsis', 'cardiogenic'):
    s = der[der.shock_type == sh]
    ms = smf.logit('hospital_death ~ bili_ge10 + age + sofa', data=s).fit(disp=0)
    o, ci, p = cfg.or_ci(ms)
    print(f'  {sh:27s} aOR {o} ({ci})  n={int(ms.nobs)}')

# over-adjustment: SOFA contains a bilirubin component (docs/variables.md)
m0 = smf.logit(cfg.F_DER_NO_SOFA, data=der).fit(disp=0)
o, ci, _ = cfg.or_ci(m0)
rho = spearmanr(der.sofa, der.peak_bili).statistic
print(f'  omitting SOFA               aOR {o} ({ci})')
print(f'  Spearman(first-day SOFA, peak bilirubin) = {rho:.3f}')

# Table 3: the full pooled model
NAMES = {'bili_ge10': 'Bilirubin ≥10 mg/dL', 'age': 'Age (per year)',
         'sofa': 'SOFA (per point)',
         'C(shock_type, Treatment(reference="sepsis"))[T.cardiogenic]': 'Cardiogenic (vs septic)',
         'C(shock_type, Treatment(reference="sepsis"))[T.hypovolemic]': 'Hypovolemic (vs septic)'}
mt = smf.logit('hospital_death ~ bili_ge10 + age + sofa + '
               'C(shock_type, Treatment(reference="sepsis"))', data=der).fit(disp=0)
rows = []
for k, lab in NAMES.items():
    o, ci, p = cfg.or_ci(mt, k)
    rows.append([lab, o, ci, cfg.pfmt(p)])
cfg.save(pd.DataFrame(rows, columns=['Variable', 'aOR', '95% CI', 'p']), 'table3')

# collinearity
X = pd.get_dummies(der[['bili_ge10', 'age', 'sofa', 'shock_type']],
                   columns=['shock_type'], drop_first=True).astype(float)
X = sm.add_constant(X)
vifs = [variance_inflation_factor(X.values, i)
        for i, c in enumerate(X.columns) if c != 'const']
print(f'  variance inflation factors  {min(vifs):.2f} to {max(vifs):.2f}')

print()
print('=' * 70)
print('VALIDATION COHORT (eICU-CRD)')
print('=' * 70)

eic = cfg.load('cohort_eicu').merge(
    cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score',
                           'apache_iva_predmort', 'aps_bilirubin',
                           'aps_vent', 'aps_dialysis', 'hospitalid']],
    on='patientunitstayid', how='left')

for lab, f in [('age + etiology', 'hospital_death ~ bili_ge10 + age + C(shock_type)'),
               ('+ APACHE IVa score', cfg.F_VAL),
               ('+ APACHE IVa predicted mortality',
                'hospital_death ~ bili_ge10 + age + C(shock_type) + apache_iva_predmort'),
               ('+ ventilation and dialysis',
                cfg.F_VAL + ' + aps_vent + aps_dialysis')]:
    mm = smf.logit(f, data=eic).fit(disp=0)
    o, ci, p = cfg.or_ci(mm)
    print(f'  {lab:33s} aOR {o} ({ci})  n={int(mm.nobs)}')

# clustering by hospital
dd = eic.dropna(subset=['apache_iva_score', 'age'])
mc = smf.logit(cfg.F_VAL, data=dd).fit(disp=0, cov_type='cluster',
                                     cov_kwds={'groups': dd.hospitalid})
o, ci, p = cfg.or_ci(mc)
print(f'  {"cluster-robust by hospital":33s} aOR {o} ({ci})')
print(f'  hospitals contributing: {eic.hospitalid.nunique()}')

# is admission severity different at all?
hi, lo = eic[eic.bili_ge10 == 1], eic[eic.bili_ge10 == 0]
for col, lab in [('apache_iva_score', 'APACHE IVa score'),
                 ('apache_iva_predmort', 'predicted mortality')]:
    a, b = hi[col].dropna(), lo[col].dropna()
    print(f'  {lab:33s} {cfg.iqr(a, 0 if "score" in col else 3)} vs '
          f'{cfg.iqr(b, 0 if "score" in col else 3)}  p={cfg.pfmt(mannwhitneyu(a, b).pvalue)}')

# does the delayed peak add anything to first-day bilirubin?
s = eic.dropna(subset=['aps_bilirubin', 'apache_iva_score', 'age'])
mb = smf.logit(cfg.F_VAL + ' + aps_bilirubin', data=s).fit(disp=0)
o, ci, p = cfg.or_ci(mb, 'aps_bilirubin')
print(f'  first-day bilirubin, per mg/dL    OR {o} ({ci})  p={cfg.pfmt(p)}')
s2 = s[~((s.bili_ge10 == 1) & (s.aps_bilirubin >= 10))]
m2 = smf.logit(cfg.F_VAL, data=s2).fit(disp=0)
o, ci, _ = cfg.or_ci(m2)
print(f'  excluding day-1 bilirubin >=10    aOR {o} ({ci})  n={int(m2.nobs)}')

# discrimination
for f, lab in [(cfg.F_VAL.replace('bili_ge10 + ', ''), 'base'),
               (cfg.F_VAL, 'base + bilirubin >=10')]:
    mm = smf.logit(f, data=dd).fit(disp=0)
    print(f'  AUC, {lab:28s} {roc_auc_score(dd.hospital_death, mm.predict(dd)):.3f}')
pm = eic.dropna(subset=['apache_iva_predmort'])
print(f'  AUC, APACHE IVa alone             '
      f'{roc_auc_score(pm.hospital_death, pm.apache_iva_predmort):.3f}')
