# -*- coding: utf-8 -*-
"""Sensitivity analyses: dose-response, competing causes, the definition of
hypoxic hepatitis, and the transplant exclusion (eTables 4 to 7).
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr
import config as cfg  # not "C": patsy uses C() for categorical terms

der = cfg.load('cohort_mimic')
eic = cfg.load('cohort_eicu').merge(
    cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score',
                           'dx_chronic_liver', 'dx_biliary_obstruction',
                           'dx_hemolysis', 'dx_rhabdomyolysis']],
    on='patientunitstayid', how='left')
der = der.merge(cfg.load('landmark_mimic')[[
    'hadm_id', 'dx_chronic_liver', 'dx_biliary_obstruction',
    'dx_hemolysis', 'dx_rhabdomyolysis']], on='hadm_id', how='left')

SPECS = [('Derivation', der, cfg.F_DER, ' + sofa'),
         ('Validation', eic, cfg.F_VAL, '')]

# ------------------------------------------------- eTable 4: dose-response
rows = []
for name, d, _, extra in SPECS:
    d = d.copy()
    # Half-open on the left, so a peak of exactly 10.0 is exposed and the
    # categorical and binary exposures agree (docs/exclusions.md).
    d['cat'] = pd.cut(d.peak_bili, [0, 3, 10, np.inf],
                      labels=['<3', '3–9.9', '≥10'], right=False)
    f = (f'hospital_death ~ C(cat, Treatment(reference="<3")) + age{extra} '
         f'+ C(shock_type)' + (' + apache_iva_score' if name == 'Validation' else ''))
    m = smf.logit(f, data=d).fit(disp=0)
    ref = d[d.cat == '<3']
    rows.append([name, '<3', len(ref), f'{100 * ref.hospital_death.mean():.1f}',
                 '1.00 (reference)', '', ''])
    for k in m.params.index:
        if 'cat' in k:
            lab = k.split('T.')[1].rstrip(']')
            s = d[d.cat == lab]
            o, ci, p = cfg.or_ci(m, k)
            rows.append([name, lab, len(s), f'{100 * s.hospital_death.mean():.1f}',
                         o, ci, cfg.pfmt(p)])
    d['log2b'] = np.log2(d.peak_bili.clip(lower=0.1))
    m2 = smf.logit(f.replace('C(cat, Treatment(reference="<3"))', 'log2b'),
                   data=d).fit(disp=0)
    o, ci, p = cfg.or_ci(m2, 'log2b')
    rows.append([name, 'Per doubling', len(d), '', o, ci, cfg.pfmt(p)])
cfg.save(pd.DataFrame(rows, columns=['Cohort', 'Peak bilirubin (mg/dL)', 'n',
                                   'Mortality (%)', 'aOR', '95% CI', 'p']),
       'etable4_dose_response')

# --------------------------------------------- eTable 5: competing causes
freq, sens = [], []
DX = [('dx_chronic_liver', 'Chronic liver disease'),
      ('dx_biliary_obstruction', 'Biliary obstruction'),
      ('dx_hemolysis', 'Hemolysis'),
      ('dx_rhabdomyolysis', 'Rhabdomyolysis')]
EXCL = [('None (primary analysis)', lambda x: x.index == x.index),
        ('Excluding chronic liver disease', lambda x: x.dx_chronic_liver == 0),
        ('Excluding biliary obstruction', lambda x: x.dx_biliary_obstruction == 0),
        ('Excluding rhabdomyolysis', lambda x: x.dx_rhabdomyolysis == 0),
        ('Excluding all of the above',
         lambda x: (x.dx_chronic_liver == 0) & (x.dx_biliary_obstruction == 0)
                   & (x.dx_rhabdomyolysis == 0) & (x.dx_hemolysis == 0))]
for name, d, f, _ in SPECS:
    for c, lab in DX:
        freq.append([name, lab, f'{100 * d[d.bili_ge10 == 1][c].mean():.1f}',
                     f'{100 * d[d.bili_ge10 == 0][c].mean():.1f}',
                     f'{100 * d[c].mean():.1f}'])
    for lab, fn in EXCL:
        s = d[fn(d)]
        m = smf.logit(f, data=s).fit(disp=0)
        o, ci, p = cfg.or_ci(m)
        sens.append([name, lab, int(m.nobs), o, ci, cfg.pfmt(p)])
cfg.save(pd.DataFrame(freq, columns=['Cohort', 'Coded diagnosis',
                                   'Bilirubin ≥10 (%)', 'Bilirubin <10 (%)',
                                   'Overall (%)']), 'etable5a_competing_freq')
cfg.save(pd.DataFrame(sens, columns=['Cohort', 'Exclusion', 'n', 'aOR', '95% CI', 'p']),
       'etable5b_competing_sensitivity')

# ------------------------------------------- eTable 6: definition of HH
dm = der.merge(cfg.load('labs_mimic'), on='hadm_id', how='left')
dm['alt_uln'] = cfg.uln(dm.gender, 'alt')
ev = eic.merge(cfg.load('labs_eicu'), on='patientunitstayid', how='left') \
        .merge(cfg.load('chol_eicu')[['patientunitstayid', 'gender']],
               on='patientunitstayid', how='left')
ev['alt_uln'] = cfg.uln(ev.gender, 'alt')
DEFS = {
    'Derivation': (dm, cfg.F_DER, [
        ('ALT >500 U/L (primary)', lambda x: x.alt_peak_all > 500),
        ('ALT >20 × ULN', lambda x: x.alt_peak_all > 20 * x.alt_uln),
        ('ALT >1000 U/L', lambda x: x.alt_peak_all > 1000),
        ('ALT >2000 U/L', lambda x: x.alt_peak_all > 2000),
        ('ALT >500 U/L and LDH >400 U/L',
         lambda x: (x.alt_peak_all > 500) & (x.ldh_peak_all > 400))]),
    'Validation': (ev, cfg.F_VAL, [
        ('ALT >500 U/L (primary)', lambda x: x.alt_peak > 500),
        ('ALT >20 × ULN', lambda x: x.alt_peak > 20 * x.alt_uln),
        ('ALT >1000 U/L', lambda x: x.alt_peak > 1000),
        ('ALT >2000 U/L', lambda x: x.alt_peak > 2000),
        ('ALT >500 U/L and LDH >400 U/L',
         lambda x: (x.alt_peak > 500) & (x.ldh_peak > 400))]),
}
rows = []
for name, (d, f, defs) in DEFS.items():
    for lab, fn in defs:
        s = d[fn(d).fillna(False)]
        m = smf.logit(f, data=s).fit(disp=0)
        o, ci, p = cfg.or_ci(m)
        rows.append([name, lab, len(s), f'{100 * s.hospital_death.mean():.1f}',
                     f'{100 * s.bili_ge10.mean():.1f}', o, ci, cfg.pfmt(p)])
cfg.save(pd.DataFrame(rows, columns=['Cohort', 'Definition of hypoxic hepatitis',
                                   'n', 'Mortality (%)', 'Bilirubin ≥10 (%)',
                                   'aOR', '95% CI', 'p']), 'etable6_hh_definition')

# ------------------------------- eTable 7b: the transplant exclusion checked
tx = cfg.load('tx_check')
tx['gap_d'] = tx.ast_peak_day_from_admit - tx.tx_day_from_admit
rho = spearmanr(tx.tx_day_from_admit, tx.ast_peak_day_from_admit).statistic
rows = [['Admissions excluded for transplant-related reasons', len(tx)],
        ['Spearman correlation, transplant date vs AST peak day', f'{rho:.3f}'],
        ['AST peak within ±2 days of the procedure',
         f'{int((tx.gap_d.abs() <= 2).sum())} ({100 * (tx.gap_d.abs() <= 2).mean():.1f}%)'],
        ['AST peak >2 days after the procedure', int((tx.gap_d > 2).sum())],
        ['AST peak >2 days before the procedure', int((tx.gap_d < -2).sum())]]
cfg.save(pd.DataFrame(rows, columns=['', 'Value']), 'etable7b_transplant_check')
print(f'  transplant exclusion: rho = {rho:.3f}, '
      f'{int((tx.gap_d.abs() <= 2).sum())}/{len(tx)} within two days')
