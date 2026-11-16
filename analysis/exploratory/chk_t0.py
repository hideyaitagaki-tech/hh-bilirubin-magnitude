import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, config as cfg
for name, path in [('Derivation','/home/claude/alp_mimic.csv'), ('Validation','/home/claude/alp_eicu.csv')]:
    d = pd.read_csv(path).dropna(subset=['hours_ast_to_bili'])
    h = d.hours_ast_to_bili
    print(f'--- {name}: n={len(d)}')
    print(f'    whole cohort median [IQR]: {h.median():.0f} [{h.quantile(.25):.0f}–{h.quantile(.75):.0f}]')
    print(f'    exactly 0 h: {(h==0).sum()} ({100*(h==0).mean():.1f}%)   <0: {(h<0).sum()} ({100*(h<0).mean():.1f}%)   >0: {(h>0).sum()} ({100*(h>0).mean():.1f}%)')
    for lab, sub in [('bili <10', d[d.bili_ge10==0]), ('bili >=10', d[d.bili_ge10==1])]:
        hh = sub.hours_ast_to_bili
        print(f'    {lab:10s} n={len(sub):5d}  median [IQR] {hh.median():.0f} [{hh.quantile(.25):.0f}–{hh.quantile(.75):.0f}]  ==0: {100*(hh==0).mean():.1f}%')
