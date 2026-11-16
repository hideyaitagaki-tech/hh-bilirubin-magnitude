import pandas as pd, numpy as np, config as cfg
lm = cfg.load('landmark_mimic')
L=168
a = lm[(lm.hours_t0_to_death.isna())|(lm.hours_t0_to_death>L)]
print('alive at 168h:', len(a))
print('  bili_max_48h missing:', int(a.bili_max_48h.isna().sum()))
print('  >=10 within 48h:', int((a.bili_max_48h>=10).sum()))
print('  <10 through 48h (incident set):', int((a.bili_max_48h<10).sum()))
print('  exposed at 168 landmark:', int((a.bili_max_168h>=10).sum()))
inc = a[a.bili_max_48h<10]
print('  incident exposed:', int((inc.bili_max_168h>=10).sum()))
