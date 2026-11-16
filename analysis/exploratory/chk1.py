import pandas as pd, numpy as np, config as cfg
ce = cfg.load('cohort_eicu'); le_ = cfg.load('labs_eicu'); che = cfg.load('chol_eicu')
for n,d in [('cohort_eicu',ce),('labs_eicu',le_),('chol_eicu',che)]:
    print(n, d.shape, list(d.columns))
