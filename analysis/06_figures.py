# -*- coding: utf-8 -*-
"""Figures 1 and 2, at 300 dpi in TIFF (for journal submission) and PNG.

Figure 1  A  time to the bilirubin peak from shock onset, by etiology
          B  in-hospital mortality by etiology and bilirubin group
          C  timing of the bilirubin peak relative to the AST peak
Figure 2  A  time to the bilirubin peak, reproduced in validation
          B  mortality, reproduced
          C  alkaline phosphatase, reproduced
          D  direct bilirubin, both cohorts

Only three hypovolemic patients reached the bilirubin threshold. Neither a
median (1A) nor a mortality proportion (1B) is drawn for that subgroup, so
that the figure does not state visually what the text declines to claim.
"""
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import config as cfg

# Liberation Serif is metric-compatible with Times New Roman, which the
# manuscript uses; the figures therefore match the text they sit in.
TITLE, PANEL, LABEL, ANNOT = 11, 10, 9, 8
plt.rcParams.update({'font.family': 'serif',
                     'font.serif': ['Liberation Serif', 'Times New Roman',
                                    'Nimbus Roman', 'DejaVu Serif'],
                     'font.size': LABEL,
                     'axes.titlesize': PANEL, 'axes.labelsize': LABEL,
                     'xtick.labelsize': LABEL, 'ytick.labelsize': LABEL,
                     'legend.fontsize': ANNOT, 'figure.titlesize': TITLE,
                     'axes.spines.top': False, 'axes.spines.right': False})
LO, HI = '#6BAED6', '#E6703B'
MR, MB = '#C0392B', '#2C6FAD'
PUR, PAL = '#6A51A3', '#CBC9E2'
ET = [('sepsis', 'Septic'), ('cardiogenic', 'Cardiogenic'), ('hypovolemic', 'Hypovolemic')]

der = cfg.load('cohort_mimic')
eic = cfg.load('cohort_eicu')
cm = cfg.load('chol_mimic')
ce = cfg.load('chol_eicu')


def q(df, sh, g, col='bili_peak_hours'):
    v = df[(df.shock_type == sh) & (df.bili_ge10 == g)][col].dropna().astype(float)
    return (np.median(v), *np.percentile(v, [25, 75]), len(v)) if len(v) else (np.nan,) * 3 + (0,)


def mort(df, sh, g):
    s = df[(df.shock_type == sh) & (df.bili_ge10 == g)]
    return 100 * s.hospital_death.mean(), len(s)


# ================================ FIGURE 1 ================================
# A and B on the top row, C spanning the bottom: a 1x3 row at this aspect
# is unreadable once scaled into a portrait text column.
fig = plt.figure(figsize=(11.5, 8.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0], hspace=0.50, wspace=0.22)
a1 = fig.add_subplot(gs[0, 0]); a2 = fig.add_subplot(gs[0, 1])
a3 = fig.add_subplot(gs[1, :])
x, w = np.arange(3), 0.36

for g, col, lab, off in [(0, LO, 'Bilirubin <10 mg/dL', -w / 2),
                         (1, HI, 'Bilirubin ≥10 mg/dL', w / 2)]:
    med, q1, q3, ns = [], [], [], []
    for sh, _ in ET:
        m, a, b, n = q(der, sh, g)
        if sh == 'hypovolemic' and g == 1:      # n=3: no median reported
            m = a = b = np.nan
        med.append(m); q1.append(a); q3.append(b); ns.append(n)
    med, q1, q3 = map(lambda v: np.array(v, float), (med, q1, q3))
    a1.bar(x + off, np.nan_to_num(med), w, color=col, label=lab, zorder=2)
    ok = ~np.isnan(med)
    a1.errorbar(x[ok] + off, med[ok], yerr=[med[ok] - q1[ok], q3[ok] - med[ok]],
                fmt='none', ecolor='#333333', capsize=3, lw=1.0, zorder=3)
    for i, (xi, m, n) in enumerate(zip(x + off, med, ns)):
        if np.isnan(m):
            a1.text(xi, 8, f'not calculated\n(n={n})', ha='center', va='bottom',
                    fontsize=ANNOT, color='#555555', style='italic')
        else:
            a1.text(xi, q3[i] + 10, f'{cfg.rh(m)} h\n(n={n})', ha='center',
                    va='bottom', fontsize=ANNOT)
a1.set_xticks(x); a1.set_xticklabels([l for _, l in ET])
a1.set_ylabel('Time to bilirubin peak (hours from shock onset)')
a1.set_title('A. Time to bilirubin peak', loc='left', fontweight='bold', fontsize=PANEL)
a1.set_ylim(0, 470); a1.legend(frameon=False, fontsize=ANNOT, loc='upper right')
a1.grid(axis='y', color='#E8E8E8', zorder=0)

for g, col, lab, off in [(0, LO, 'Bilirubin <10 mg/dL', -w / 2),
                         (1, HI, 'Bilirubin ≥10 mg/dL', w / 2)]:
    vals, ns = zip(*[mort(der, sh, g) for sh, _ in ET])
    # A group of fewer than ten patients carries no interpretable proportion,
    # and a full-height bar would state visually what the text declines to
    # claim. Panel A already withholds the median for the same group.
    shown = [np.nan if n < 10 else v for v, n in zip(vals, ns)]
    a2.bar(x + off, np.nan_to_num(shown), w, color=col, label=lab, zorder=2)
    for xi, v, n in zip(x + off, shown, ns):
        if np.isnan(v):
            a2.text(xi, 2, f'not calculated\n(n={n})', ha='center', va='bottom',
                    fontsize=ANNOT, style='italic', color='#555555')
        else:
            a2.text(xi, v + 1.5, f'{v:.1f}%\n(n={n})', ha='center', va='bottom', fontsize=ANNOT)
a2.set_xticks(x); a2.set_xticklabels([l for _, l in ET])
a2.set_ylabel('In-hospital mortality (%)'); a2.set_ylim(0, 125)
a2.set_yticks(range(0, 101, 20))
a2.set_title('B. In-hospital mortality', loc='left', fontweight='bold', fontsize=PANEL)
a2.legend(frameon=False, fontsize=ANNOT, loc='upper left')
a2.grid(axis='y', color='#E8E8E8', zorder=0)

CATS = [('Before the AST peak', lambda h: h < 0, '#9ECAE1'),
        ('Same blood draw (0 h)', lambda h: h == 0, '#4292C6'),
        ('>0 to <72 h', lambda h: (h > 0) & (h < 72), '#FDBE85'),
        ('72 to <168 h', lambda h: (h >= 72) & (h < 168), '#E6703B'),
        ('≥168 h (≥1 week)', lambda h: h >= 168, '#A63603')]
ck = cm.dropna(subset=['hours_ast_to_bili'])
for yi, g, glab in [(1, 0, 'Bilirubin <10 mg/dL'), (0, 1, 'Bilirubin ≥10 mg/dL')]:
    hh, left = ck[ck.bili_ge10 == g].hours_ast_to_bili, 0.0
    for name, fn, col in CATS:
        p = 100 * fn(hh).mean()
        a3.barh(yi, p, 0.5, left=left, color=col, zorder=2, label=name if yi == 1 else None)
        if p >= 6:
            a3.text(left + p / 2, yi, f'{p:.1f}', ha='center', va='center', fontsize=ANNOT,
                    fontweight='bold',
                    color='white' if col in ('#4292C6', '#A63603', '#E6703B') else '#333333')
        left += p
    a3.text(-2, yi, glab, ha='right', va='center', fontsize=LABEL)
a3.set_yticks([]); a3.set_xlim(0, 100); a3.set_ylim(-0.6, 1.75)
a3.set_xlabel('Patients (%)')
a3.set_title('C. Bilirubin peak relative to the AST peak', loc='left',
             fontweight='bold', fontsize=PANEL)
a3.legend(frameon=False, fontsize=ANNOT, ncol=5, loc='upper center',
          bbox_to_anchor=(0.5, -0.20), columnspacing=1.4, handlelength=1.4)
a3.spines['left'].set_visible(False)
a3.grid(axis='x', color='#E8E8E8', zorder=0)

fig.suptitle('Bilirubin dynamics and mortality in the derivation cohort '
             f'(MIMIC-IV, n={len(der):,})', fontweight='bold', fontsize=TITLE, y=0.97)
fig.savefig(cfg.OUT / 'Figure_1.tiff', dpi=300, format='tiff',
            pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig)

# ================================ FIGURE 2 ================================
fig, ((b1, b2), (b3, b4)) = plt.subplots(2, 2, figsize=(11.5, 8.6))
xs, ww = np.arange(2), 0.2
for lab, col, vals, off in [
        ('MIMIC-IV ≥10', MR, [q(der, 'sepsis', 1)[0], q(der, 'cardiogenic', 1)[0]], -1.5),
        ('eICU-CRD ≥10', HI, [q(eic, 'sepsis', 1)[0], q(eic, 'cardiogenic', 1)[0]], -0.5),
        ('MIMIC-IV <10', MB, [q(der, 'sepsis', 0)[0], q(der, 'cardiogenic', 0)[0]], 0.5),
        ('eICU-CRD <10', LO, [q(eic, 'sepsis', 0)[0], q(eic, 'cardiogenic', 0)[0]], 1.5)]:
    b1.bar(xs + off * ww, vals, ww, color=col, label=lab, zorder=2)
    for xi, v in zip(xs + off * ww, vals):
        b1.text(xi, v + 3, cfg.rh(v), ha='center', va='bottom', fontsize=ANNOT)
b1.set_xticks(xs); b1.set_xticklabels(['Septic', 'Cardiogenic'])
b1.set_ylabel('Time to bilirubin peak (h)'); b1.set_ylim(0, 215)
b1.set_title('A. Time to peak — reproduced', loc='left', fontweight='bold', fontsize=PANEL)
b1.legend(frameon=False, fontsize=ANNOT, ncol=2); b1.grid(axis='y', color='#E8E8E8', zorder=0)

import statsmodels.formula.api as smf
aors = [cfg.or_ci(smf.logit(cfg.F_DER, data=der).fit(disp=0))[0],
        cfg.or_ci(smf.logit('hospital_death ~ bili_ge10 + age + C(shock_type)',
                            data=eic).fit(disp=0))[0]]
for i, (lab, d) in enumerate([('MIMIC-IV', der), ('eICU-CRD', eic)]):
    for j, (g, col) in enumerate(zip((0, 1), (LO, HI))):
        v = 100 * d[d.bili_ge10 == g].hospital_death.mean()
        b2.bar(i + (j - 0.5) * 0.36, v, 0.36, color=col, zorder=2,
               label=('Bilirubin <10 mg/dL', 'Bilirubin ≥10 mg/dL')[j] if i == 0 else None)
        b2.text(i + (j - 0.5) * 0.36, v + 1.2, f'{v:.1f}%', ha='center', va='bottom', fontsize=ANNOT)
    b2.text(i, 4, f'aOR {aors[i]}', ha='center', fontsize=ANNOT, fontweight='bold')
b2.set_xticks([0, 1]); b2.set_xticklabels(['MIMIC-IV', 'eICU-CRD'])
b2.set_ylabel('In-hospital mortality (%)'); b2.set_ylim(0, 92)
b2.set_title('B. Mortality — reproduced', loc='left', fontweight='bold', fontsize=PANEL)
b2.legend(frameon=False, fontsize=ANNOT); b2.grid(axis='y', color='#E8E8E8', zorder=0)

for i, (lab, d) in enumerate([('MIMIC-IV', cm), ('eICU-CRD', ce)]):
    for j, (g, col) in enumerate(zip((0, 1), (LO, HI))):
        v = d[d.bili_ge10 == g].v_peak_alp.dropna().astype(float)
        q1, q2, q3 = np.percentile(v, [25, 50, 75])
        x0 = i + (j - 0.5) * 0.36
        b3.bar(x0, q2, 0.36, color=col, zorder=2,
               label=('Bilirubin <10 mg/dL', 'Bilirubin ≥10 mg/dL')[j] if i == 0 else None)
        b3.errorbar(x0, q2, yerr=[[q2 - q1], [q3 - q2]], fmt='none', ecolor='#333333',
                    capsize=3, lw=1.0, zorder=3)
        b3.text(x0, q3 + 8, f'{q2:.0f}', ha='center', va='bottom', fontsize=ANNOT,
                fontweight='bold')
b3.axhline(cfg.ULN['alp_f'], ls='--', lw=0.9, color='#888888', zorder=1)
b3.text(0.5, cfg.ULN['alp_f'] + 12, 'ULN (women)', fontsize=ANNOT,
        color='#888888', ha='center', va='bottom', zorder=4,
        bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'))
b3.set_xticks([0, 1]); b3.set_xticklabels(['MIMIC-IV', 'eICU-CRD'])
b3.set_ylabel('Peak alkaline phosphatase (U/L)'); b3.set_ylim(0, 600)
b3.set_title('C. Alkaline phosphatase — reproduced', loc='left',
             fontweight='bold', fontsize=PANEL)
b3.legend(frameon=False, fontsize=ANNOT, loc='upper left', ncol=2, bbox_to_anchor=(0.0, 1.0))
b3.grid(axis='y', color='#E8E8E8', zorder=0)


def dbil(d, g):
    s = d[d.bili_ge10 == g].dropna(subset=['dbili_at_bili_peak'])
    s = s[s.v_peak_bili > 0]
    r = s.dbili_at_bili_peak / s.v_peak_bili
    ok = (r >= 0) & (r <= 1.2)
    s, r = s[ok], r[ok]
    return np.median(s.v_peak_bili), np.median(s.dbili_at_bili_peak), np.median(r), len(s)


for i, (lab, d) in enumerate([('MIMIC-IV', cm), ('eICU-CRD', ce)]):
    for j, g in enumerate((0, 1)):
        tot, dir_, rat, n = dbil(d, g)
        x0 = i + (j - 0.5) * 0.36
        b4.bar(x0, dir_, 0.36, color=PUR, zorder=2,
               label='Direct (conjugated)' if (i == 0 and j == 0) else None)
        b4.bar(x0, tot - dir_, 0.36, bottom=dir_, color=PAL, zorder=2,
               label='Indirect' if (i == 0 and j == 0) else None)
        b4.text(x0, tot + 0.6, f'{rat:.2f}\n(n={n})', ha='center', va='bottom',
                fontsize=ANNOT, color=PUR, fontweight='bold')
b4.set_xticks([-0.18, 0.18, 0.82, 1.18])
b4.set_xticklabels(['MIMIC-IV\n<10', 'MIMIC-IV\n≥10',
                    'eICU-CRD\n<10', 'eICU-CRD\n≥10'], fontsize=ANNOT)
b4.set_ylabel('Bilirubin at the bilirubin peak (mg/dL)'); b4.set_ylim(0, 30)
b4.set_title('D. Direct bilirubin — reproduced', loc='left', fontweight='bold', fontsize=PANEL)
b4.legend(frameon=False, fontsize=ANNOT, loc='upper left'); b4.grid(axis='y', color='#E8E8E8', zorder=0)

fig.suptitle('Independent validation: timing (A), mortality (B), '
             'and the cholestatic pattern (C, D)', fontweight='bold', fontsize=TITLE)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(cfg.OUT / 'Figure_2.tiff', dpi=300, format='tiff',
            pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig)

from PIL import Image
for n in (1, 2):
    im = Image.open(cfg.OUT / f'Figure_{n}.tiff')
    im.convert('RGB').save(cfg.OUT / f'Figure_{n}.png', dpi=(300, 300))
    print(f'  wrote Figure_{n}.tiff and .png  {im.size}')
