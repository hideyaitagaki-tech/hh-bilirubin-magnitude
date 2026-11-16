import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, statsmodels.formula.api as smf, config as cfg
c=cfg.load('cohort_eicu'); ap=cfg.load('apache_eicu')
m=c.merge(ap.drop(columns=[x for x in ap.columns if x in c.columns and x!='patientunitstayid']),
          on='patientunitstayid', how='left')
print('cohort', len(m), 'exposed', int(m.bili_ge10.sum()))
sub = m.dropna(subset=['aps_bilirubin'])
print('with first-day bilirubin:', len(sub), 'exposed', int(sub.bili_ge10.sum()))
for lab,d in [('all with aps_bili', sub)]:
    hi,lo = d[d.bili_ge10==1], d[d.bili_ge10==0]
    print(f'  {lab}: first-day bili hi {cfg.iqr(hi.aps_bilirubin,1)}  lo {cfg.iqr(lo.aps_bilirubin,1)}')
# model subset = complete case on F_VAL + aps_bilirubin
cc = m.dropna(subset=['aps_bilirubin','apache_iva_score','age'])
print('complete-case (age+APACHE+first-day bili):', len(cc), 'exposed', int(cc.bili_ge10.sum()))
hi,lo = cc[cc.bili_ge10==1], cc[cc.bili_ge10==0]
print('  first-day bili hi', cfg.iqr(hi.aps_bilirubin,1), ' lo', cfg.iqr(lo.aps_bilirubin,1))
print('  exposed with first-day >=10:', int((hi.aps_bilirubin>=10).sum()))
# excluding those
ex = cc[~((cc.bili_ge10==1)&(cc.aps_bilirubin>=10))]
mm=smf.logit(cfg.F_VAL, data=ex).fit(disp=0); o,ci,p=cfg.or_ci(mm)
print(f'  excluding first-day>=10 exposed: n={int(mm.nobs)} aOR {o} ({ci}) p={cfg.pfmt(p)}')
# first-day bilirubin itself
mm2=smf.logit(cfg.F_VAL+' + aps_bilirubin', data=cc).fit(disp=0); o,ci,p=cfg.or_ci(mm2,'aps_bilirubin')
print(f'  first-day bilirubin per mg/dL: OR {o} ({ci}) p={cfg.pfmt(p)} n={int(mm2.nobs)}')
