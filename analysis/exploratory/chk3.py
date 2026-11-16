import pandas as pd, numpy as np, config as cfg
lm = cfg.load('landmark_mimic')
a = lm[(lm.hours_t0_to_death.isna())|(lm.hours_t0_to_death>168)]
ex = a[a.bili_max_168h>=10]
print('exposed at 168:', len(ex))
print('  of which >=10 by 48h:', int((ex.bili_max_48h>=10).sum()))
print('  of which <10 through 48h:', int((ex.bili_max_48h<10).sum()))
print('  of which 48h missing:', int(ex.bili_max_48h.isna().sum()))
# half of exposed at 168 statement in Discussion
