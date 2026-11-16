import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, statsmodels.formula.api as smf, config as cfg

cm = pd.read_csv('/home/claude/alp_mimic.csv')
ce = pd.read_csv('/home/claude/alp_eicu.csv')
coh = cfg.load('cohort_mimic')[['hadm_id','sofa','age','shock_type','hospital_death','bili_ge10']]
ap  = cfg.load('apache_eicu')[['patientunitstayid','apache_iva_score']]

for name, d in [('Derivation', cm), ('Validation', ce)]:
    ex = d[d.bili_ge10==1].dropna(subset=['hours_ast_to_bili']).copy()
    h = ex.hours_ast_to_bili
    print(f'--- {name}: exposed with interval available n={len(ex)}')
    print(f'    negative: {(h<0).sum()} ({100*(h<0).mean():.1f}%)   zero: {(h==0).sum()} ({100*(h==0).mean():.1f}%)   positive: {(h>0).sum()}')
    print(f'    median [IQR]: {h.median():.0f} [{h.quantile(.25):.0f}-{h.quantile(.75):.0f}]  min {h.min():.0f} max {h.max():.0f}')
    print(f'    survivors median {h[ex.hospital_death==0].median():.0f}   decedents {h[ex.hospital_death==1].median():.0f}')
