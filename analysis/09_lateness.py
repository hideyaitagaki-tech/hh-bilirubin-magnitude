# -*- coding: utf-8 -*-
"""Lateness of the bilirubin peak within the exposed group (eTable 13).

Lateness is the interval from the aminotransferase peak to the bilirubin peak.
It is only defined where the bilirubin peak actually followed: in a little over
a quarter of exposed patients the two peaks fell in the same blood draw or the
bilirubin peak came first, and a logarithm of such an interval does not exist.
An earlier version of this analysis clipped those intervals to one hour, which
assigned them the smallest possible positive lateness and produced a spurious
inverse association; the continuous model is therefore fitted only where the
interval is positive.

Lateness cannot be separated from how long each patient was observed, because a
late peak can only be recorded in a patient who was still being sampled. The
observation span is added in a secondary model to show the size and direction
of that dependence. The span is itself determined by survival, so that model is
not a corrected estimate of lateness; it is a demonstration that the two cannot
be told apart in this design.

Writes eTable 13.
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
import config as cfg  # not "C": patsy uses C() for categorical terms

der = (cfg.load('chol_mimic')  # carries sofa already; merging it again renames it
       .merge(cfg.load('labs_mimic')[['hadm_id', 'obs_span_hr']], on='hadm_id', how='left'))
val = (cfg.load('chol_eicu')
       .merge(cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score']],
              on='patientunitstayid', how='left')
       .merge(cfg.load('labs_eicu')
              .assign(obs_span_hr=lambda x: x.last_lab_hr - x.first_lab_hr)
              [['patientunitstayid', 'obs_span_hr']],
              on='patientunitstayid', how='left'))

COV = {'Derivation': 'age + sofa + C(shock_type)',
       'Validation': 'age + C(shock_type) + apache_iva_score'}

rows = []
for name, d in [('Derivation', der), ('Validation', val)]:
    cov = COV[name]
    ex = d[d.bili_ge10 == 1].dropna(subset=['hours_ast_to_bili']).copy()
    ex['log2peak'] = np.log2(ex.v_peak_bili)
    ex['late168'] = (ex.hours_ast_to_bili >= 168).astype(int)
    ex['log2span'] = np.log2(ex.obs_span_hr.clip(lower=1))
    # The clipped specification, reported only to show where the earlier
    # inverse association came from.
    ex['log2late_clipped'] = np.log2(ex.hours_ast_to_bili.clip(lower=1))

    pos = ex[ex.hours_ast_to_bili > 0].copy()
    pos['log2late'] = np.log2(pos.hours_ast_to_bili)

    n_pos = len(pos)
    n_zero = int((ex.hours_ast_to_bili == 0).sum())
    n_neg = int((ex.hours_ast_to_bili < 0).sum())
    print(f'{name}: exposed with an interval {len(ex)} — '
          f'positive {n_pos} ({100*n_pos/len(ex):.1f}%), '
          f'zero {n_zero} ({100*n_zero/len(ex):.1f}%), '
          f'negative {n_neg} ({100*n_neg/len(ex):.1f}%)')
    s, dd = pos[pos.hospital_death == 0], pos[pos.hospital_death == 1]
    p_mw = mannwhitneyu(s.hours_ast_to_bili, dd.hours_ast_to_bili).pvalue
    print(f'   interval, positive only: {cfg.iqr(pos.hours_ast_to_bili, 0)} h; '
          f'survivors {cfg.iqr(s.hours_ast_to_bili, 0)}, '
          f'decedents {cfg.iqr(dd.hours_ast_to_bili, 0)}, p={cfg.pfmt(p_mw)}')

    SPECS = [
        ('Lateness per doubling, intervals >0 h (primary)',
         pos, f'hospital_death ~ log2late + log2peak + {cov}', 'log2late'),
        ('  Peak magnitude per doubling, same model',
         pos, f'hospital_death ~ log2late + log2peak + {cov}', 'log2peak'),
        ('Lateness ≥168 h, intervals >0 h',
         pos, f'hospital_death ~ late168 + log2peak + {cov}', 'late168'),
        ('Lateness per doubling, observation span added',
         pos, f'hospital_death ~ log2late + log2peak + log2span + {cov}', 'log2late'),
        ('Lateness ≥168 h, all exposed (sensitivity)',
         ex, f'hospital_death ~ late168 + log2peak + {cov}', 'late168'),
        ('Lateness per doubling, intervals ≤0 h clipped to 1 h (sensitivity)',
         ex, f'hospital_death ~ log2late_clipped + log2peak + {cov}', 'log2late_clipped'),
    ]
    for label, data, f, term in SPECS:
        m = smf.logit(f, data=data).fit(disp=0)
        o, ci, p = cfg.or_ci(m, term)
        rows.append([name, label, int(m.nobs), o, ci, cfg.pfmt(p)])

t = pd.DataFrame(rows, columns=['Cohort', 'Model', 'n', 'aOR', '95% CI', 'p'])
print()
print(t.to_string(index=False))
cfg.save(t, 'etable13_lateness')
