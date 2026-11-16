# -*- coding: utf-8 -*-
"""Two analyses that address the objection that a peak taken over the whole
admission is confounded by how long each patient was observed.

(1) Observation window. Patients who reach a bilirubin of 10 mg/dL are
    observed for longer, so the exposed group is selected for long
    trajectories. The direction of that selection is tested here three ways:
    restricting to progressively longer minimum observation spans, entering
    the span in the model, and the landmark analysis, in which every patient
    is observed to the same point by construction.

(2) Incident hyperbilirubinemia. Exposure at a late landmark may simply carry
    forward patients who were already hyperbilirubinemic early. The incident
    analysis restricts to patients whose bilirubin stayed below 10 mg/dL
    through the first 48 h, so the exposure is newly developed.

Writes eTables 8 and 9.
"""
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
import config as cfg  # not "C": patsy uses C() for categorical terms

der = cfg.load('cohort_mimic').merge(cfg.load('labs_mimic'), on='hadm_id', how='left')
eic = (cfg.load('cohort_eicu')
       .merge(cfg.load('labs_eicu'), on='patientunitstayid', how='left')
       .merge(cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score']],
              on='patientunitstayid', how='left'))
eic['obs_span_hr'] = eic.last_lab_hr - eic.first_lab_hr
SPECS = [('Derivation', der, cfg.F_DER), ('Validation', eic, cfg.F_VAL)]


def row(name, label, d, f):
    m = smf.logit(f, data=d).fit(disp=0)
    o, ci, p = cfg.or_ci(m)
    return [name, label, len(d), int(d.bili_ge10.sum()), o, ci, cfg.pfmt(p)]


# ------------------------------- eTable 8: the observation window
rows = []
for name, d, f in SPECS:
    for thr in [0, 72, 168, 336]:
        lab = 'None (primary analysis)' if thr == 0 else f'≥{thr} h of laboratory follow-up'
        rows.append(row(name, lab, d[d.obs_span_hr >= thr], f))
    s = d.copy()
    s['log2span'] = np.log2(s.obs_span_hr.clip(lower=1))
    rows.append(row(name, 'Observation span entered in the model', s, f + ' + log2span'))
cfg.save(pd.DataFrame(rows, columns=['Cohort', 'Restriction', 'n', 'Exposed, n',
                                     'aOR', '95% CI', 'p']),
         'etable8_observation_window')

# The span itself, and mortality across its quartiles: the selection runs the
# other way, because an early death truncates observation.
for name, d, f in SPECS:
    s = d.copy()
    s['log2span'] = np.log2(s.obs_span_hr.clip(lower=1))
    m = smf.logit(f + ' + log2span', data=s).fit(disp=0)
    o, ci, p = cfg.or_ci(m, 'log2span')
    print(f'  {name}: observation span, per doubling  aOR {o} ({ci})')
    s['q'] = pd.qcut(s.obs_span_hr, 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
    for q, g in s.groupby('q'):
        print(f'    {q}  span {cfg.iqr(g.obs_span_hr, 0)} h  '
              f'death {100 * g.hospital_death.mean():.1f}%  '
              f'bilirubin ≥10 {100 * g.bili_ge10.mean():.1f}%')

# Rate of rise: peak bilirubin divided by the time taken to reach it. If the
# association were only a matter of having had longer to accumulate, the peak
# should not survive adjustment for how fast it was reached.
rows = []
for name, d, f in SPECS:
    s = d[d.bili_peak_hours > 0].copy()
    s['rate'] = s.peak_bili / (s.bili_peak_hours / 24)
    hi, lo = s[s.bili_ge10 == 1], s[s.bili_ge10 == 0]
    s['log2rate'] = np.log2(s.rate.clip(lower=0.01))
    m = smf.logit(f + ' + log2rate', data=s).fit(disp=0)
    ob, cb, pb = cfg.or_ci(m, 'bili_ge10')
    orr, cr, pr = cfg.or_ci(m, 'log2rate')
    rows.append([name, cfg.iqr(hi.rate, 2), cfg.iqr(lo.rate, 2),
                 cfg.pfmt(mannwhitneyu(hi.rate.dropna(), lo.rate.dropna()).pvalue),
                 f'{ob} ({cb})', cfg.pfmt(pb), f'{orr} ({cr})', cfg.pfmt(pr)])
cfg.save(pd.DataFrame(rows, columns=[
    'Cohort', 'Rate of rise, ≥10 (mg/dL/day)', 'Rate of rise, <10 (mg/dL/day)', 'p',
    'aOR for peak ≥10, rate in the model', 'p ', 'aOR per doubling of rate', 'p  ']),
         'etable8b_rate_of_rise')

# ------------------------------- eTable 9: incident hyperbilirubinemia
# A late landmark may simply carry forward patients who were already
# hyperbilirubinemic. The incident rows restrict to patients whose bilirubin
# stayed below 10 mg/dL through the first 48 h.
lm = cfg.load('landmark_mimic').merge(cfg.load('sofa_landmark'), on='hadm_id', how='left')
le = (cfg.load('landmark_eicu')
      .merge(cfg.load('cohort_eicu')[['patientunitstayid', 'shock_type', 'age',
                                      'hospital_death']],
             on='patientunitstayid', how='left')
      .merge(cfg.load('apache_eicu')[['patientunitstayid', 'apache_iva_score']],
             on='patientunitstayid', how='left'))


def alive_at(d, L):
    return d[(d.hours_t0_to_death.isna()) | (d.hours_t0_to_death > L)]


rows = []
for name, d, f in [('Derivation', lm, cfg.F_DER), ('Validation', le, cfg.F_VAL)]:
    for L in [72, 96, 168]:
        col = f'bili_max_{L}h'
        for lab, prevalent in [
                ('All patients alive at the landmark', None),
                ('Bilirubin <10 mg/dL through the first 48 h (incident)', False),
                ('Bilirubin ≥10 mg/dL within the first 48 h (prevalent)', True)]:
            s = alive_at(d, L).copy()
            if prevalent is True:
                # exposure is having crossed 10 mg/dL by 48 h, whole cohort kept
                s['bili_ge10'] = (s.bili_max_48h >= 10).astype(int)
            else:
                if prevalent is False:
                    s = s[s.bili_max_48h < 10]
                s['bili_ge10'] = (s[col] >= 10).astype(int)
            if s.bili_ge10.sum() < 5:
                continue
            m = smf.logit(f, data=s).fit(disp=0)
            o, ci, p = cfg.or_ci(m)
            e, u = s[s.bili_ge10 == 1], s[s.bili_ge10 == 0]
            rows.append([name, f'{L} h', lab, len(s), len(e),
                         f'{100 * e.hospital_death.mean():.1f}',
                         f'{100 * u.hospital_death.mean():.1f}', o, ci, cfg.pfmt(p)])
cfg.save(pd.DataFrame(rows, columns=['Cohort', 'Landmark', 'Exposure definition', 'n',
                                     'Exposed, n', 'Mortality, exposed (%)',
                                     'Mortality, unexposed (%)', 'aOR', '95% CI', 'p']),
         'etable9_incident_hyperbilirubinemia')

for name, d in [('Derivation', lm), ('Validation', le)]:
    e = alive_at(d, 168)
    e = e[e.bili_max_168h >= 10]
    print(f'  {name}: exposed at the 168 h landmark {len(e)}; '
          f'already ≥10 mg/dL within 48 h {int((e.bili_max_48h >= 10).sum())}; '
          f'newly developed between 48 and 168 h {int((e.bili_max_48h < 10).sum())}')

# ------------------- eTable 10: time-updated organ dysfunction (derivation)
# The landmark models adjust for the SOFA score at time zero, which cannot
# capture organ failure developing in the meantime. The total SOFA score
# contains bilirubin as its liver subscore, so adjusting for it would adjust
# in part for the exposure; the non-hepatic score is the score minus that
# subscore.
rows = []
for L in [48, 72, 96, 168]:
    s = alive_at(lm, L).copy()
    s['bili_ge10'] = (s[f'bili_max_{L}h'] >= 10).astype(int)
    s['sofa_tu'] = s[f'sofa_max_{L}h']
    s['sofa_nh'] = s[f'sofa_nh_{L}h']
    for lab, f in [
            ('SOFA at time zero (primary)', cfg.F_DER),
            ('Time-updated SOFA, total score',
             'hospital_death ~ bili_ge10 + age + sofa_tu + C(shock_type)'),
            ('Time-updated SOFA, non-hepatic score',
             'hospital_death ~ bili_ge10 + age + sofa_nh + C(shock_type)')]:
        m = smf.logit(f, data=s).fit(disp=0)
        o, ci, p = cfg.or_ci(m)
        rows.append([f'{L} h', lab, len(s), int(s.bili_ge10.sum()), o, ci, cfg.pfmt(p)])
    # the most conservative specification: incident exposure, window matched by
    # the landmark, and concurrent non-hepatic organ failure adjusted for
    si = s[s.bili_max_48h < 10]
    if si.bili_ge10.sum() >= 5:
        m = smf.logit('hospital_death ~ bili_ge10 + age + sofa_nh + C(shock_type)',
                      data=si).fit(disp=0)
        o, ci, p = cfg.or_ci(m)
        rows.append([f'{L} h', 'Incident exposure and time-updated non-hepatic SOFA',
                     len(si), int(si.bili_ge10.sum()), o, ci, cfg.pfmt(p)])
cfg.save(pd.DataFrame(rows, columns=['Landmark', 'Adjustment', 'n', 'Exposed, n',
                                     'aOR', '95% CI', 'p']),
         'etable10_time_updated_sofa')
