import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, statsmodels.formula.api as smf, config as cfg
cm = pd.read_csv('/home/claude/alp_mimic.csv')
ce = pd.read_csv('/home/claude/alp_eicu.csv').merge(
     cfg.load('apache_eicu')[['patientunitstayid','apache_iva_score']], on='patientunitstayid', how='left')
def fit(d,f,t):
    m=smf.logit(f,data=d).fit(disp=0); return cfg.or_ci(m,t)+(int(m.nobs),)
for name,d,cov in [('Derivation',cm,'age + sofa + C(shock_type)'),
                   ('Validation',ce,'age + C(shock_type) + apache_iva_score')]:
    ex=d[d.bili_ge10==1].dropna(subset=['hours_ast_to_bili']).copy()
    ex['log2peak']=np.log2(ex.v_peak_bili)
    ex['log2late']=np.log2(ex.hours_ast_to_bili.clip(lower=1))
    for thr_lab,col in [('>=168',(ex.hours_ast_to_bili>=168)),('>168',(ex.hours_ast_to_bili>168))]:
        ex['b']=col.astype(int)
        o,ci,p,n=fit(ex,'hospital_death ~ b + log2peak','b')
        o2,ci2,p2,n2=fit(ex,f'hospital_death ~ b + log2peak + {cov}','b')
        print(f'{name} binary {thr_lab}: 2-term {o} ({ci}) p={cfg.pfmt(p)} n={n} | adjusted {o2} ({ci2}) p={cfg.pfmt(p2)} n={n2}')
    o,ci,p,n=fit(ex,f'hospital_death ~ log2late + log2peak + {cov}','log2late')
    om,cim,pm,_=fit(ex,f'hospital_death ~ log2late + log2peak + {cov}','log2peak')
    print(f'{name} log2late adjusted: {o} ({ci}) p={cfg.pfmt(p)} n={n}; magnitude {om} ({cim}) p={cfg.pfmt(pm)}')
    pos=ex[ex.hours_ast_to_bili>0].copy(); pos['log2late']=np.log2(pos.hours_ast_to_bili)
    o,ci,p,n=fit(pos,f'hospital_death ~ log2late + log2peak + {cov}','log2late')
    om,cim,pm,_=fit(pos,f'hospital_death ~ log2late + log2peak + {cov}','log2peak')
    print(f'{name} POSITIVE-ONLY adjusted: late {o} ({ci}) p={cfg.pfmt(p)} n={n}; magnitude {om} ({cim}) p={cfg.pfmt(pm)}')
