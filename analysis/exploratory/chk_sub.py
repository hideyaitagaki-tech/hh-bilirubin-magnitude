import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
from scipy.stats import spearmanr
cm = pd.read_csv('/home/claude/alp_mimic.csv'); ce = pd.read_csv('/home/claude/alp_eicu.csv')
for name,d in [('Derivation',cm),('Validation',ce)]:
    ex = d[d.bili_ge10==1].dropna(subset=['hours_ast_to_bili'])
    tot = d[d.bili_ge10==1]
    late = ex[ex.hours_ast_to_bili>0]
    print(f'{name}: whole exposed n={len(tot)} mort {100*tot.hospital_death.mean():.1f}%')
    print(f'   bilirubin peak AFTER AST peak n={len(late)} mort {100*late.hospital_death.mean():.1f}%')
    sim = ex[ex.hours_ast_to_bili<=0]
    print(f'   peak at/before AST peak n={len(sim)} mort {100*sim.hospital_death.mean():.1f}%')
    l7 = ex[ex.hours_ast_to_bili>=168]
    print(f'   >=168 h after AST peak n={len(l7)} mort {100*l7.hospital_death.mean():.1f}%')
d=pd.read_csv('/home/claude/repo/data/cohort_mimic.csv')
print('rho SOFA/bili %.5f' % spearmanr(d.sofa,d.peak_bili).statistic)
s=d.dropna(subset=['peak_inr','total_tx_ml'])
print('rho INR/tx   %.5f' % spearmanr(s.peak_inr,s.total_tx_ml).statistic)
