# -*- coding: utf-8 -*-
"""Paths, loading and the shared conventions every analysis script uses.

The analysis reads only the CSVs exported from the tables that `sql/` writes.
Put them in `data/` under the names below, or set HH_DATA to another directory.
"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP

ROOT = Path(__file__).resolve().parent.parent
DATA = Path(os.environ.get('HH_DATA', ROOT / 'data'))
OUT = Path(os.environ.get('HH_OUT', ROOT / 'output'))
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    'cohort_mimic':    'cohort_mimic.csv',
    'cohort_eicu':     'cohort_eicu.csv',
    'labs_mimic':      'labs_mimic.csv',
    'labs_eicu':       'labs_eicu.csv',
    'chol_mimic':      'chol_mimic.csv',
    'chol_eicu':       'chol_eicu.csv',
    'landmark_mimic':  'landmark_mimic.csv',
    'apache_eicu':     'apache_eicu.csv',
    'tx_check':        'tx_check.csv',
    'hemorrhagic':     'cohort_hemorrhagic_mimic.csv',
    'landmark_eicu':   'landmark_eicu.csv',
    'sofa_landmark':   'sofa_landmark_mimic.csv',
    'alt_timing_mimic':'alt_timing_mimic.csv',
    'alt_timing_eicu': 'alt_timing_eicu.csv',
    'competing_mimic': 'competing_extended_mimic.csv',
    'competing_eicu':  'competing_extended_eicu.csv',
}

# eICU uses -1 as a missing-value sentinel in the APACHE tables. Converted in
# SQL as well; repeated here because loading a stale export otherwise biases
# every summary of these variables silently.
SENTINEL_NEG1 = ['apache_iva_score', 'apache_iva_predmort', 'apache_any_score',
                 'aps_bilirubin', 'aps_creatinine', 'aps_meanbp', 'aps_ph', 'aps_pao2']

# Upper limits of normal (docs/variables.md)
ULN = {'alt_f': 33, 'alt_m': 45, 'alp_f': 105, 'alp_m': 130}

# Model formulas. The derivation model adjusts for SOFA, the validation model
# for APACHE IVa; see docs/variables.md on why creatinine and INR are not
# entered separately.
F_DER = 'hospital_death ~ bili_ge10 + age + sofa + C(shock_type)'
F_VAL = 'hospital_death ~ bili_ge10 + age + C(shock_type) + apache_iva_score'
F_DER_NO_SOFA = 'hospital_death ~ bili_ge10 + age + C(shock_type)'


def load(key):
    """Read one exported table, applying the sentinel conversion."""
    df = pd.read_csv(DATA / FILES[key])
    for c in SENTINEL_NEG1:
        if c in df.columns:
            df.loc[df[c] == -1, c] = np.nan
    return df


def uln(gender, which):
    """Sex-specific upper limit of normal; accepts MIMIC 'M'/'F' and eICU words."""
    female = pd.Series(gender).isin(['F', 'Female', 'female'])
    return np.where(female, ULN[f'{which}_f'], ULN[f'{which}_m'])


# ---------------------------------------------------------------- formatting
def rh(x, dec=0):
    """Round half up. Python rounds half to even, which turned a SOFA of 8.5
    into 8 and an LDH of x.5 into x in an earlier version of these tables."""
    return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-dec),
                                                rounding=ROUND_HALF_UP))


def med(s, dec=0):
    return rh(np.median(pd.Series(s).dropna().astype(float)), dec)


def iqr(s, dec=0):
    s = pd.Series(s).dropna().astype(float)
    q1, q2, q3 = np.percentile(s, [25, 50, 75])
    return f'{rh(q2, dec)} [{rh(q1, dec)}–{rh(q3, dec)}]'


def pfmt(p):
    return '<0.001' if p < 0.001 else f'{p:.3f}'


def or_ci(model, var='bili_ge10'):
    """Odds ratio, confidence interval and p for one coefficient."""
    b, se = model.params[var], model.bse[var]
    return (f'{np.exp(b):.2f}',
            f'{np.exp(b - 1.96 * se):.2f}–{np.exp(b + 1.96 * se):.2f}',
            model.pvalues[var])


def save(df, name):
    path = OUT / f'{name}.csv'
    df.to_csv(path, index=False)
    print(f'  wrote {path.relative_to(ROOT)}')
    return path
