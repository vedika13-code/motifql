#!/usr/bin/env python3
"""Draws the manuscript figures from the text files in results/.

Usage: python3 tools/make_figures.py results/ figures/
"""
import os, re, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'   # validated categorical slots 1-4
INK, MUTED, GRID = '#1f2328', '#5b6470', '#e3e6ea'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.linewidth': 0.8, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'legend.frameon': False, 'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.facecolor': 'white',
})

def fig_pipeline(out):
    fig, ax = plt.subplots(figsize=(7.0, 2.6))
    ax.set_xlim(0, 100); ax.set_ylim(0, 38); ax.axis('off'); ax.grid(False)
    def box(x, y, w, h, title, sub, accent=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3,rounding_size=1.0',
                     fc='#fbede3' if accent else '#eef2f7', ec=ORANGE if accent else '#1f3a5f', lw=1.1))
        ax.text(x + w / 2, y + h * 0.68, title, ha='center', va='center', fontsize=7.6, weight='bold',
                color='#a8461a' if accent else '#1f3a5f')
        ax.text(x + w / 2, y + h * 0.3, sub, ha='center', va='center', fontsize=5.9, color=INK)
    def arrow(x1, y1, x2, y2, style='->'):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle=style, color=MUTED, lw=1.0))
    y1, h = 22, 11
    box(0.5, y1, 15.5, h, 'MotifQL', 'motifs, budgets,\ncontracts')
    box(20, y1, 16, h, 'Front end', 'lexer, parser,\n6 semantic passes')
    box(40.5, y1, 17, h, 'Automaton', 'Thompson NFA, error\nlayers, minimal DFA')
    box(62, y1, 16, h, 'Matcher', 'table-driven scan\nof FASTA')
    box(82.5, y1, 17, h, 'Hits', 'positions, and\nexpected by chance')
    for a, b in [(15.8, 19.7), (36.3, 40.2), (57.8, 61.7), (78.3, 82.2)]:
        arrow(a, y1 + h / 2, b, y1 + h / 2)
    y2 = 3
    box(28, y2, 18, h, 'Specificity', 'chance hits per position\nfrom the DFA (or a scan)', True)
    box(50, y2, 18, h, 'Sensitivity', 'diverged or known sites\nfound by the DFA', True)
    box(72, y2, 21, h, 'Budget synthesis', 'best sensitivity under the\ncontract; monotone pruning', True)
    arrow(49, y1 - 0.4, 37, y2 + h + 0.4); arrow(49, y1 - 0.4, 59, y2 + h + 0.4)
    arrow(46.3, y2 + h / 2, 49.7, y2 + h / 2); arrow(68.3, y2 + h / 2, 71.7, y2 + h / 2)
    ax.annotate('', xy=(55, y1 + 0.2), xytext=(82.5, y2 + h + 0.3),
                arrowprops=dict(arrowstyle='->', color=ORANGE, lw=1.0, connectionstyle='arc3,rad=0.25'))
    ax.text(77, 17.4, 'chosen budgets', fontsize=6.3, color='#a8461a')
    fig.savefig(os.path.join(out, 'fig1_pipeline.png')); plt.close(fig)

def parse_chr21(path):
    txt = open(path).read()
    part1 = txt.split('Part 1')[1].split('Part 2')[0]
    sur = []
    for line in part1.splitlines():
        m = re.match(r'^(\w+)\s+\S+\s+(exact|allow \d+ \w+)\s+(\d+)\s+(\d+)\s', line)
        if m: sur.append((m.group(1), int(m.group(3)), int(m.group(4))))
    part2 = txt.split('Part 2')[1].split('median')[0]
    real = []
    for line in part2.splitlines():
        m = re.match(r'^(\w+)\s+\S+\s+(exact|allow \d+ \w+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
        if m:
            obs = int(m.group(3)); r = [float(m.group(i)) for i in range(4, 9)]
            real.append((m.group(1), m.group(2), obs, r))
    return sur, real

def fig_rates(results, out):
    sur, real = parse_chr21(os.path.join(results, 'eval_chr21.txt'))
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), sharex=False)
    ax = axes[0]
    pts = [(p, o) for _, p, o in sur if p >= 5 and o >= 5]
    ax.loglog([p for p, _ in pts], [o for _, o in pts], 'o', ms=5, color=BLUE, mec='white', mew=0.8)
    lo, hi = 3, 4e5
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.9, ls='--')
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel('predicted chance hits (order-1 model)'); ax.set_ylabel('observed hits')
    ax.set_title('a  Sequence drawn from the model', loc='left', fontsize=8.5, color=INK, weight='bold')
    ax = axes[1]
    for col, color, marker, label in [(1, ORANGE, 's', 'order 1'), (4, AQUA, '^', 'empirical')]:
        xs, ys = [], []
        for name, tol, obs, r in real:
            if obs >= 20 and r[col] > 0:
                xs.append(obs / r[col]); ys.append(obs)
        ax.loglog(xs, ys, marker, ms=5, color=color, mec='white', mew=0.8, label=label)
    ax.plot([lo, hi * 3], [lo, hi * 3], color=MUTED, lw=0.9, ls='--')
    ax.set_xlim(10, 2e6); ax.set_ylim(10, 2e6)
    ax.set_xlabel('predicted chance hits'); ax.set_ylabel('observed hits on real chr21')
    ax.set_title('b  Real chr21, held-out half', loc='left', fontsize=8.5, color=INK, weight='bold')
    ax.legend(loc='upper left', fontsize=7.5, handletextpad=0.3)
    gc = [(obs / r[1], obs) for name, tol, obs, r in real if name == 'GC_box']
    if gc: ax.annotate('GC box', xy=gc[0], xytext=(gc[0][0] * 2.5, gc[0][1] * 3), fontsize=7, color=INK,
                       arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    fig.tight_layout(w_pad=2.5)
    fig.savefig(os.path.join(out, 'fig2_rates.png')); plt.close(fig)

def parse_e2(path):
    txt = open(path).read().split('E2:')[1]
    blocks = re.split(r'chance-hit rates from: ', txt)[1:]
    out = {}
    for b in blocks:
        mode = b.split(' (')[0]
        rows = []
        for line in b.splitlines()[2:]:
            m = re.match(r'^(.{26}) (.{32}) +([\d.]+) +([\d.]+) +([\d.]+) +([\d.]+)', line)
            if m: rows.append((m.group(1).strip(), m.group(2).strip(), float(m.group(3)), float(m.group(5)), float(m.group(6))))
        out[mode] = rows
    return out

def fig_planted(results, out):
    e2 = parse_e2(os.path.join(results, 'eval_chr21.txt'))
    rows = e2.get('empirical', [])
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    exact = next((r for r in rows if r[0].startswith('exact')), None)
    for label, chosen, recall, false, pred in rows:
        synth = label.startswith('synth')
        x = max(false, 0.02)
        same_as_exact = exact and not label.startswith('exact') and abs(recall - exact[2]) < 1e-9 and abs(false - exact[3]) < 1e-9
        if same_as_exact:
            ax.semilogx(x, recall, 'D', ms=7, color=ORANGE, mec='white', mew=0.8)
            continue
        if label.startswith('exact') and any(r[0].startswith('synth') and abs(r[2] - recall) < 1e-9 and abs(r[3] - false) < 1e-9 for r in rows):
            ax.semilogx(x, recall, 'o', ms=10, color=BLUE, mec='white', mew=0.8, zorder=1)
            ax.annotate('exact (0, 0) = compiler, <= 2/100kb', (x, recall), xytext=(x * 1.4, recall + 0.035), fontsize=6.6, color=INK)
            continue
        ax.semilogx(x, recall, 'D' if synth else 'o', ms=7 if synth else 6, color=ORANGE if synth else BLUE,
                    mec='white', mew=0.8)
        text = label.replace('hand-picked ', '').replace('synthesised, ', 'compiler, ')
        if synth and chosen and chosen != '-': text += f'\n({chosen.replace(" edits", "").replace(" edit", "")})'
        dx, ha = (1.4, 'left') if x < 300 else (0.7, 'right')
        ax.annotate(text, (x, recall), xytext=(x * dx, recall + 0.035), fontsize=6.6, color=INK, ha=ha)
    ax.set_xlabel('false hits per 100 kb (log)'); ax.set_ylabel('planted sites found')
    ax.set_ylim(0, 1.12); ax.set_xlim(0.01, 3e4)
    ax.set_title('Planted promoters, real chr21', loc='left', fontsize=8.5, color=INK, weight='bold')
    fig.savefig(os.path.join(out, 'fig3_planted.png')); plt.close(fig)

def parse_pas(path):
    rows = []
    for line in open(path):
        m = re.match(r'^(.{40}) +([\d.]+)±([\d.]+) +([\d.]+) +([\d.]+) +([\d.]+)', line)
        if m: rows.append((m.group(1).strip(), float(m.group(2)), float(m.group(3)), float(m.group(4))))
    return rows

def fig_pas(results, out):
    rows = parse_pas(os.path.join(results, 'eval_pas.txt'))
    fig, ax = plt.subplots(figsize=(4.6, 3.3))
    learned = sorted([(kb, rec, sd) for lab, rec, sd, kb in rows if lab.startswith('learned budget + positions')])
    seen = []
    for kb, rec, sd in learned:
        if (round(kb, 2), round(rec, 3)) not in seen: seen.append((round(kb, 2), round(rec, 3)))
    pwm = sorted([(kb, rec) for lab, rec, sd, kb in rows if lab.startswith('PWM')])
    ax.plot([k for k, _ in pwm], [r for _, r in pwm], '-^', ms=5, color=YELLOW, mec='white', mew=0.8, lw=1.4,
            label='PWM from the same training sites')
    ax.errorbar([k for k, _ in seen], [r for _, r in seen], fmt='-D', ms=5, color=ORANGE, mec='white', mew=0.8,
                lw=1.6, label='learned budget + positions')
    fixed_style = {'exact AATAAA': (BLUE, 'o'), 'AATAAA, 1 mismatch': (BLUE, 'o'), 'AATAAA, 1 edit': (BLUE, 'o'),
                   'expert: AATAAA | ATTAAA': (AQUA, 's'), 'expert: 12 known variants': (AQUA, 's')}
    names = {'exact AATAAA': 'exact', 'AATAAA, 1 mismatch': '1 mismatch', 'AATAAA, 1 edit': '1 edit',
             'expert: AATAAA | ATTAAA': '2 hexamers', 'expert: 12 known variants': '12 known variants'}
    for lab, rec, sd, kb in rows:
        if lab in fixed_style:
            c, mk = fixed_style[lab]
            ax.plot(kb, rec, mk, ms=6, color=c, mec='white', mew=0.8)
            off = {'12 known variants': (-62, 6), '1 mismatch': (6, -11), 'exact': (6, 4), '2 hexamers': (-12, 8),
                   '1 edit': (6, 4)}[names[lab]]
            ax.annotate(names[lab], (kb, rec), xytext=off, textcoords='offset points', fontsize=6.8, color=INK)
    ax.plot([], [], 'o', color=BLUE, label='uniform budget (hand-picked)')
    ax.plot([], [], 's', color=AQUA, label='literature variant lists')
    ax.set_xscale('log'); ax.set_xlim(0.8, 40); ax.set_ylim(0.45, 0.92)
    ax.set_xlabel('hits per kb in held-out control windows (log)')
    ax.set_ylabel("held-out 3' ends with a hit")
    ax.set_title('Polyadenylation signals at annotated 3\' ends', loc='left', fontsize=8.5, color=INK, weight='bold')
    ax.legend(loc='lower right', fontsize=7)
    fig.savefig(os.path.join(out, 'fig4_pas.png')); plt.close(fig)

def main():
    results, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    fig_pipeline(out); fig_rates(results, out); fig_planted(results, out); fig_pas(results, out)
    print('figures written to', out)

if __name__ == '__main__':
    main()
