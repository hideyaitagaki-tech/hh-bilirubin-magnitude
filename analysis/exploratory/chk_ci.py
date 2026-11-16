import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, config as cfg
d = cfg.load('cohort_mimic'); e = cfg.load('cohort_eicu')
for name, x in [('Derivation', d), ('Validation', e)]:
    ex = x[x.bili_ge10==1]
    b = pd.cut(ex.peak_bili, [10,15,20,30,1e9], right=False,
               labels=['10-14.9','15-19.9','20-29.9','>=30'])
    print(f'--- {name}: exposed {len(ex)}, survivors {int((ex.hospital_death==0).sum())} ({100*(ex.hospital_death==0).mean():.1f}%)')
    for lab,g in ex.groupby(b, observed=True):
        print(f'    {lab:9s} n={len(g):4d}  mortality {100*g.hospital_death.mean():.1f}%')
# monitoring: still <10 at 72h, later crossed
lm = cfg.load('landmark_mimic')
s = lm[(lm.bili_max_72h < 10)]
cross = s[s.bili_max_all >= 10]
print(f'\nDerivation: below 10 through 72 h n={len(s)}; later crossed n={len(cross)} ({100*len(cross)/len(s):.1f}%)')
print(f'   mortality crossed {100*cross.hospital_death.mean():.1f}%  vs not crossed {100*s[s.bili_max_all<10].hospital_death.mean():.1f}%')
le = cfg.load('landmark_eicu').merge(cfg.load('cohort_eicu')[['patientunitstayid','hospital_death']], on='patientunitstayid', how='left', suffixes=('','_c'))
s2 = le[le.bili_max_72h < 10]; c2 = s2[s2.bili_max_all>=10]
print(f'Validation: below 10 through 72 h n={len(s2)}; later crossed n={len(c2)} ({100*len(c2)/len(s2):.1f}%)')
print(f'   mortality crossed {100*c2.hospital_death.mean():.1f}%  vs not {100*s2[s2.bili_max_all<10].hospital_death.mean():.1f}%')
