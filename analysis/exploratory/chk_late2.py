import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, statsmodels.formula.api as smf, config as cfg

cm = pd.read_csv('/home/claude/alp_mimic.csv')
ce = pd.read_csv('/home/claude/alp_eicu.csv')
coh = cfg.load('cohort_mimic')[['hadm_id','sofa']]
ap  = cfg.load('apache_eicu')[['patientunitstayid','apache_iva_score']]
ap.loc[ap.apache_iva_score==-1,'apache_iva_score']=np.nan

def fit(d, f, term):
    m = smf.logit(f, data=d).fit(disp=0)
    return cfg.or_ci(m, term) + (int(m.nobs),)

cm = cm.merge(coh, on='hadm_id', how='left')
ce = ce.merge(ap, on='patientunitstayid', how='left')

for name, d, cov in [('Derivation', cm, 'age + sofa + C(shock_type)'),
                     ('Validation', ce, 'age + C(shock_type) + apache_iva_score')]:
    ex = d[d.bili_ge10==1].dropna(subset=['hours_ast_to_bili']).copy()
    ex['log2late'] = np.log2(ex.hours_ast_to_bili.clip(lower=1))
    ex['log2peak'] = np.log2(ex.v_peak_bili)
    ex['late168']  = (ex.hours_ast_to_bili >= 168).astype(int)
    print(f'=== {name} (n={len(ex)})')
    for lab, f, t in [
        ('log2late + log2peak only',      'hospital_death ~ log2late + log2peak', 'log2late'),
        ('  ...magnitude term',           'hospital_death ~ log2late + log2peak', 'log2peak'),
        ('log2late + log2peak + covars',  f'hospital_death ~ log2late + log2peak + {cov}', 'log2late'),
        ('  ...magnitude term',           f'hospital_death ~ log2late + log2peak + {cov}', 'log2peak'),
        ('late168 + log2peak only',       'hospital_death ~ late168 + log2peak', 'late168'),
        ('late168 + log2peak + covars',   f'hospital_death ~ late168 + log2peak + {cov}', 'late168'),
    ]:
        try:
            o,ci,p,n = fit(ex, f, t)
            print(f'   {lab:32s} {t:9s} aOR {o} ({ci}) p={cfg.pfmt(p)}  n={n}')
        except Exception as e:
            print('   ERR', lab, e)
    # positive-only sensitivity
    pos = ex[ex.hours_ast_to_bili>0].copy()
    pos['log2late']=np.log2(pos.hours_ast_to_bili)
    o,ci,p,n = fit(pos, 'hospital_death ~ log2late + log2peak', 'log2late')
    print(f'   {"positive intervals only":32s} log2late  aOR {o} ({ci}) p={cfg.pfmt(p)}  n={n}')
