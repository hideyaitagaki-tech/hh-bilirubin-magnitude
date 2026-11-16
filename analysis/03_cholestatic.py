# -*- coding: utf-8 -*-
"""The cholestatic pattern (Results, "Biochemical pattern"; eTable 7).

Three observations, all reported in both cohorts:

  1. peak alkaline phosphatase is higher in the high-bilirubin group;
  2. the ALP peak is delayed, and only in that group;
  3. peak ALT does not differ, while the AST/ALT ratio does.

Together these say the high-bilirubin group is not characterised by greater
hepatocellular injury but by a shift toward a cholestatic profile. The R-ratio
classification is reported too, and is *not* concordant between cohorts; the
manuscript states that limitation rather than omitting the analysis.
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import config as cfg  # not "C": patsy uses C() for categorical terms

COHORTS = [('Derivation', 'chol_mimic', 'labs_mimic', 'hadm_id',
            'alt_peak_all', 'ast_peak_all'),
           ('Validation', 'chol_eicu', 'labs_eicu', 'patientunitstayid',
            'alt_peak', 'ast_peak')]

alp_rows, r_rows, d_rows = [], [], []

for name, chol_key, labs_key, idcol, altcol, astcol in COHORTS:
    d = cfg.load(chol_key).merge(cfg.load(labs_key), on=idcol, how='left')
    d['alp_uln'] = cfg.uln(d.gender, 'alp')
    d['alt_uln'] = cfg.uln(d.gender, 'alt')
    hi, lo = d[d.bili_ge10 == 1], d[d.bili_ge10 == 0]

    print('=' * 70); print(name)

    # 1 and 2: alkaline phosphatase, magnitude and timing
    for g, lab in [(hi, 'bilirubin >=10'), (lo, 'bilirubin <10')]:
        v = g.v_peak_alp.dropna()
        alp_rows.append([name, lab, len(v), cfg.iqr(v), cfg.iqr(v / g.alp_uln.iloc[0], 1),
                         f'{100 * (v > g.alp_uln.iloc[0]).mean():.1f}',
                         f'{100 * (v > 2 * g.alp_uln.iloc[0]).mean():.1f}',
                         f'{100 * (v > 3 * g.alp_uln.iloc[0]).mean():.1f}'])
    p = mannwhitneyu(hi.v_peak_alp.dropna(), lo.v_peak_alp.dropna()).pvalue
    print(f'  peak ALP        {cfg.iqr(hi.v_peak_alp)} vs {cfg.iqr(lo.v_peak_alp)}  '
          f'p={cfg.pfmt(p)}')
    ra = d.dropna(subset=['hours_ast_to_alp'])
    print(f'  ALP peak, hours after the AST peak  '
          f'{cfg.med(ra[ra.bili_ge10 == 1].hours_ast_to_alp)} vs '
          f'{cfg.med(ra[ra.bili_ge10 == 0].hours_ast_to_alp)}')

    # 3: ALT does not differ, the AST/ALT ratio does
    d['ast_alt'] = d[astcol] / d[altcol]
    hi, lo = d[d.bili_ge10 == 1], d[d.bili_ge10 == 0]
    for col, lab in [(altcol, 'peak ALT'), (astcol, 'peak AST'), ('ast_alt', 'AST/ALT ratio')]:
        a, b = hi[col].dropna(), lo[col].dropna()
        dec = 2 if col == 'ast_alt' else 0
        print(f'  {lab:14s}  {cfg.iqr(a, dec)} vs {cfg.iqr(b, dec)}  '
              f'p={cfg.pfmt(mannwhitneyu(a, b).pvalue)}')

    # R ratio
    r = d.dropna(subset=['alp_at_bili_peak', 'alt_at_bili_peak'])
    r = r[(r.alp_at_bili_peak > 0) & (r.alt_at_bili_peak > 0)]
    r = r.assign(R=(r.alt_at_bili_peak / r.alt_uln) / (r.alp_at_bili_peak / r.alp_uln))
    for g, lab in [(1, 'Bilirubin ≥10 mg/dL'), (0, 'Bilirubin <10 mg/dL')]:
        s = r[r.bili_ge10 == g]
        r_rows.append([name, lab, len(s), cfg.rh(s.R.median(), 1),
                       f'{100 * (s.R <= 2).mean():.1f}',
                       f'{100 * ((s.R > 2) & (s.R < 5)).mean():.1f}',
                       f'{100 * (s.R >= 5).mean():.1f}'])

    # direct / total bilirubin at the bilirubin peak
    b = d.dropna(subset=['dbili_at_bili_peak'])
    b = b[b.v_peak_bili > 0]
    b = b.assign(ratio=b.dbili_at_bili_peak / b.v_peak_bili)
    b = b[(b.ratio >= 0) & (b.ratio <= 1.2)]
    for g, lab in [(1, 'Bilirubin ≥10 mg/dL'), (0, 'Bilirubin <10 mg/dL')]:
        s = b[b.bili_ge10 == g]
        d_rows.append([name, lab, len(s), cfg.iqr(s.ratio, 2),
                       cfg.iqr(s.dbili_at_bili_peak, 1)])
    print(f'  direct bilirubin measured in {len(b)} of {len(d)} '
          f'({100 * len(b) / len(d):.1f}%)')
    print()

cfg.save(pd.DataFrame(alp_rows, columns=[
    'Cohort', 'Group', 'n', 'Peak ALP (U/L)', 'x ULN',
    '>ULN (%)', '>2x ULN (%)', '>3x ULN (%)']), 'alp_by_group')
cfg.save(pd.DataFrame(r_rows, columns=[
    'Cohort', 'Group', 'n', 'Median R ratio', 'Cholestatic, R ≤2 (%)',
    'Mixed, R 2–5 (%)', 'Hepatocellular, R ≥5 (%)']), 'etable7a_r_ratio')
cfg.save(pd.DataFrame(d_rows, columns=[
    'Cohort', 'Group', 'n', 'Direct/total ratio',
    'Direct bilirubin (mg/dL)']), 'direct_bilirubin')
