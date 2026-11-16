# -*- coding: utf-8 -*-
"""Figure 3: the interval from the aminotransferase peak to the bilirubin
peak, and what it does and does not predict.

A, B  Distribution of the interval within the exposed group of each cohort.
      The region at or below zero is shaded, because an interval there has no
      logarithm and is the part of the data that an earlier version of this
      analysis clipped to one hour. The long right tail is the delayed
      phenotype the paper is about.
C     The models fitted within the exposed group, on one axis. Lateness sits
      on the null in both cohorts and under both parameterizations, while the
      magnitude of the peak does not. The clipped specification is drawn last,
      apart, to show where the apparent inverse association came from.

Writes Figure_3.tiff and Figure_3.png at 300 dpi.
"""
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import config as cfg  # not "C": patsy uses C() for categorical terms

TITLE, PANEL, LABEL, ANNOT = 11, 10, 9, 8
plt.rcParams.update({'font.family': 'serif',
                     'font.serif': ['Liberation Serif', 'Times New Roman',
                                    'Nimbus Roman', 'DejaVu Serif'],
                     'font.size': LABEL,
                     'axes.titlesize': PANEL, 'axes.labelsize': LABEL,
                     'xtick.labelsize': LABEL, 'ytick.labelsize': LABEL,
                     'legend.fontsize': ANNOT, 'figure.titlesize': TITLE,
                     'axes.spines.top': False, 'axes.spines.right': False})
HI, GREY, MR, MB = '#E6703B', '#B0B0B0', '#C0392B', '#2C6FAD'

der = pd.read_csv('/home/claude/alp_mimic.csv')
val = pd.read_csv('/home/claude/alp_eicu.csv')

fig = plt.figure(figsize=(7.1, 6.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.42, wspace=0.24,
                      left=0.085, right=0.975, top=0.93, bottom=0.085)

# ---------------------------------------------------------------- A and B
for col, (name, d, letter) in enumerate(
        [('Derivation (MIMIC-IV)', der, 'A'), ('Validation (eICU-CRD)', val, 'B')]):
    ax = fig.add_subplot(gs[0, col])
    h = d[d.bili_ge10 == 1].hours_ast_to_bili.dropna()
    n_pos, n_zero, n_neg = int((h > 0).sum()), int((h == 0).sum()), int((h < 0).sum())
    lo, hi = -240, 480
    bins = np.arange(lo, hi + 24, 24)
    counts, edges = np.histogram(h.clip(lo, hi - 1), bins=bins)
    # The patients whose two peaks fell in the same blood draw sit in the bin
    # that starts at zero, alongside genuinely positive intervals of under a
    # day. That bin is drawn in two tones so the figure says what the counts do.
    z = int(np.where(edges[:-1] == 0)[0][0])
    for i, c in enumerate(counts):
        if i == z:
            ax.bar(edges[i], n_zero, width=24, align='edge', color=GREY,
                   edgecolor='white', linewidth=0.4, zorder=2)
            ax.bar(edges[i], c - n_zero, width=24, align='edge', bottom=n_zero,
                   color=HI, edgecolor='white', linewidth=0.4, zorder=2)
        else:
            ax.bar(edges[i], c, width=24, align='edge', color=HI,
                   edgecolor='white', linewidth=0.4, zorder=2)
    ax.axvspan(lo, 0, color=GREY, alpha=0.28, zorder=1, linewidth=0)
    ax.axvline(0, color='#555555', linewidth=0.9, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_xticks([-240, -120, 0, 120, 240, 360, 480])
    ax.set_xticklabels(['≤−240', '−120', '0', '120', '240', '360', '≥480'])
    ax.set_xlabel('Hours from the aminotransferase peak', labelpad=2)
    ax.set_ylabel('Patients' if col == 0 else '')
    ax.set_title(f'{letter}. {name}', loc='left', pad=6)
    top = ax.get_ylim()[1]
    ax.text(-132, top * 0.95, f'before\n{n_neg} ({100 * n_neg / len(h):.0f}%)',
            ha='center', va='top', fontsize=ANNOT, color='#444444')
    ax.text(-132, top * 0.60, f'same draw\n{n_zero} ({100 * n_zero / len(h):.0f}%)',
            ha='center', va='top', fontsize=ANNOT, color='#444444')
    ax.annotate('', xy=(8, top * 0.44), xytext=(-70, top * 0.50),
                arrowprops=dict(arrowstyle='-', color='#888888', linewidth=0.7))
    ax.text(300, top * 0.93, f'after\n{n_pos} ({100 * n_pos / len(h):.0f}%)\n'
            f'median {h[h > 0].median():.0f} h',
            ha='center', va='top', fontsize=ANNOT, color='#444444')

# ---------------------------------------------------------------- C
ax = fig.add_subplot(gs[1, :])
_b = ax.get_position()
ax.set_position([0.265, _b.y0, 0.975 - 0.265, _b.height])  # room for the row labels
ROWS = [
    ('Lateness, per doubling',        (0.95, 0.70, 1.30), (0.72, 0.49, 1.07)),
    ('Lateness, ≥168 h',              (1.14, 0.49, 2.65), (0.61, 0.18, 2.11)),
    ('Peak magnitude, per doubling',  (3.13, 1.60, 6.10), (3.10, 1.11, 8.61)),
]
ARTIFACT = ('Lateness, non-positive\nintervals clipped to 1 h', 0.87, 0.77, 0.97)

ticks, labels, y = [], [], 0.0
for lab, d_est, v_est in ROWS[::-1]:
    for (o, lo_, hi_), c in [(v_est, MR), (d_est, MB)]:
        ax.plot([lo_, hi_], [y, y], color=c, linewidth=1.5, solid_capstyle='butt', zorder=2)
        ax.plot([o], [y], 'o', color=c, markersize=5.2, zorder=3)
        y += 0.62
    ticks.append(y - 0.93); labels.append(lab)
    y += 0.62
y_art = -1.05
ax.plot([ARTIFACT[2], ARTIFACT[3]], [y_art, y_art], color='#777777',
        linewidth=1.5, solid_capstyle='butt', zorder=2)
ax.plot([ARTIFACT[1]], [y_art], 'o', color='#777777', markersize=5.2, zorder=3)
ticks.append(y_art); labels.append(ARTIFACT[0])
ax.axvline(1, color='#555555', linewidth=0.9, linestyle='--', zorder=1)
ax.set_xscale('log')
ax.set_xlim(0.15, 12)
ax.set_xticks([0.25, 0.5, 1, 2, 4, 8])
ax.set_xticklabels(['0.25', '0.5', '1', '2', '4', '8'])
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
ax.set_xlabel('Adjusted odds ratio for in-hospital death (log scale)', labelpad=3)
ax.set_title('C. Models fitted within the exposed group', loc='left', pad=6)
ax.axhline(-0.45, color='#CCCCCC', linewidth=0.8, zorder=0)
for c, lab in [(MB, 'Derivation'), (MR, 'Validation')]:
    ax.plot([], [], 'o-', color=c, label=lab, markersize=5.2, linewidth=1.5)
ax.legend(loc='lower right', frameon=False, ncol=2, bbox_to_anchor=(1.0, -0.02))
ax.set_ylim(-1.6, y - 0.3)

fig.savefig(cfg.OUT / 'Figure_3.tiff', dpi=300, format='tiff',
            pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig)
im = Image.open(cfg.OUT / 'Figure_3.tiff')
im.convert('RGB').save(cfg.OUT / 'Figure_3.png', dpi=(300, 300))
print(f'  wrote Figure_3.tiff and .png  {im.size}')
