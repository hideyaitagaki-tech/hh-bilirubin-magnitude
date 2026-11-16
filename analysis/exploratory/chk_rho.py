import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, config as cfg
from scipy.stats import spearmanr
d = cfg.load('cohort_mimic')
print('n cohort', len(d))
s = d.dropna(subset=['peak_inr','total_tx_ml'])
print('INR vs total_tx_ml: rho=%.3f n=%d' % (spearmanr(s.peak_inr, s.total_tx_ml).statistic, len(s)))
s2 = d.dropna(subset=['peak_inr','prbc_ml'])
print('INR vs prbc_ml:     rho=%.3f n=%d' % (spearmanr(s2.peak_inr, s2.prbc_ml).statistic, len(s2)))
s3 = d.dropna(subset=['sofa','peak_bili'])
print('SOFA vs peak_bili:  rho=%.3f n=%d' % (spearmanr(s3.sofa, s3.peak_bili).statistic, len(s3)))
