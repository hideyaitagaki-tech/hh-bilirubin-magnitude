import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, statsmodels.formula.api as smf, config as cfg

cm = pd.read_csv('/home/claude/alp_mimic.csv').merge(
     cfg.load('labs_mimic')[['hadm_id','obs_span_hr']], on='hadm_id', how='left')
ce = pd.read_csv('/home/claude/alp_eicu.csv').merge(
     cfg.load('apache_eicu')[['patientunitstayid','apache_iva_score']], on='patientunitstayid', how='left').merge(
     cfg.load('labs_eicu').assign(obs_span_hr=lambda x: x.last_lab_hr-x.first_lab_hr)[['patientunitstayid','obs_span_hr']], on='patientunitstayid', how='left')

def fit(d,f,t):
    m=smf.logit(f,data=d).fit(disp=0); return cfg.or_ci(m,t)+(int(m.nobs),)

for name,d,cov in [('Derivation',cm,'age + sofa + C(shock_type)'),
                   ('Validation',ce,'age + C(shock_type) + apache_iva_score')]:
    ex=d[d.bili_ge10==1].dropna(subset=['hours_ast_to_bili']).copy()
    ex['log2peak']=np.log2(ex.v_peak_bili)
    ex['span']=np.log2(ex.obs_span_hr.clip(lower=1))
    ex['late168']=(ex.hours_ast_to_bili>=168).astype(int)
    pos=ex[ex.hours_ast_to_bili>0].copy(); pos['log2late']=np.log2(pos.hours_ast_to_bili)
    print(f'=== {name}  exposed n={len(ex)}  positive-interval n={len(pos)}')
    for lab,d2,f,t in [
      ('continuous, adjusted',          pos, f'hospital_death ~ log2late + log2peak + {cov}','log2late'),
      ('  + observation span',          pos, f'hospital_death ~ log2late + log2peak + span + {cov}','log2late'),
      ('  magnitude in that model',     pos, f'hospital_death ~ log2late + log2peak + {cov}','log2peak'),
      ('binary >=168h, adjusted (all)', ex,  f'hospital_death ~ late168 + log2peak + {cov}','late168'),
      ('binary >=168h, positive only',  pos, f'hospital_death ~ late168 + log2peak + {cov}','late168'),
    ]:
        o,ci,p,n=fit(d2,f,t)
        print(f'   {lab:32s} {o} ({ci})  p={cfg.pfmt(p)}  n={n}')
    print(f'   median interval, positive only: {cfg.iqr(pos.hours_ast_to_bili,0)} h')
    print(f'      survivors {cfg.iqr(pos[pos.hospital_death==0].hours_ast_to_bili,0)}  decedents {cfg.iqr(pos[pos.hospital_death==1].hours_ast_to_bili,0)}')
