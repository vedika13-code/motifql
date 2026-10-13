#!/usr/bin/env python3
"""E3: polyadenylation signals at annotated transcript 3' ends.

Usage: python3 tools/eval_pas.py chr21.fa annotation.gtf [n_controls_per_site] [seed]

Real sites, not planted ones. Around 85% of human mRNA 3' ends carry a
polyadenylation signal (PAS) 10-30 nt upstream of the cleavage site:
AATAAA in about 58%, ATTAAA in about 15%, and ten other single-base
variants of AATAAA in about 15% (Beaudoing et al., Genome Res 2000).
So "AATAAA with one error" is a biologically meaningful pattern, and
the question is which error budget finds real signals without flooding
the output with chance hits.

Data
  positives   the 60 nt upstream of every distinct annotated 3' end of
              a protein-coding or lncRNA transcript (strand-oriented)
  controls    windows of the same length at random positions inside the
              same transcripts, on the transcript's strand, at least
              1 kb from any annotated 3' end; split into a TRAINING set
              (only ever used by the compiler, via --empirical, to
              measure chance-hit rates) and a TEST set (used here to
              measure false hits). Positives are never shown to the
              compiler.

Reported per configuration
  recall      share of positive windows with at least one hit
  ctrl/kb     hits per kb in held-out control windows
  ctrl hit    share of control windows with at least one hit
  enrichment  recall / ctrl hit
"""
import os, random, re, subprocess, sys, tempfile, math

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'build', 'motifql'))
COMP = str.maketrans('ACGTN', 'TGCAN')
UPSTREAM = 60

def revcomp(s): return s.translate(COMP)[::-1]

def read_fasta(path):
    with open(path) as f:
        f.readline()
        return f.read().replace('\n', '').upper()

def parse_gtf(path):
    tx = []
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'transcript': continue
        tt = re.search('transcript_type "([^"]+)"', f[8])
        if not tt or tt.group(1) not in ('protein_coding', 'lncRNA'): continue
        tx.append((f[6], int(f[3]), int(f[4])))
    return tx

def window(seq, strand, end, up):
    """`up` nt ending at the 3' end (1-based, inclusive), oriented 5'->3' on the transcript."""
    if strand == '+':
        w = seq[end - up - 1:end]
    else:
        w = revcomp(seq[end - 1:end + up])
    return w if len(w) == up + 1 and 'N' not in w else None

def build_sets(seq, tx, up, per_site, rng):
    ends = sorted({(s, b if s == '+' else a) for s, a, b in tx})
    pos = [(f'pos{i}', w) for i, (s, e) in enumerate(ends) if (w := window(seq, s, e, up))]
    end_set = {}
    for s, e in ends: end_set.setdefault(s, []).append(e)
    controls = []
    tries = 0
    while len(controls) < 2 * per_site * len(pos) and tries < 200 * per_site * len(pos):
        tries += 1
        s, a, b = rng.choice(tx)
        if b - a < 3000: continue
        e = rng.randint(a + up + 1, b - 1)
        if any(abs(e - x) < 1000 for x in end_set[s]): continue
        w = window(seq, s, e, up)
        if w: controls.append(w)
    half = len(controls) // 2
    train = [(f'tr{i}', w) for i, w in enumerate(controls[:half])]
    test = [(f'ct{i}', w) for i, w in enumerate(controls[half:])]
    return pos, train, test

def write_records(path, recs):
    with open(path, 'w') as f:
        for name, w in recs: f.write(f'>{name}\n{w}\n')

def run(motif_text, tmp, extra=()):
    with open(os.path.join(tmp, 'e.motif'), 'w') as f: f.write(motif_text)
    return subprocess.run([BIN, '--symbols', '--empirical', 'train.fa'] + list(extra) + ['e.motif'],
                          cwd=tmp, capture_output=True, text=True)

def hits_by_record(out, motif_name, fasta):
    block = re.search(rf'^{motif_name} in "{re.escape(fasta)}".*?-> \d+ hit\(s\)', out, re.S | re.M)
    if not block: return None
    return [m for m in re.findall(r'^\s+(\w+):\d+-\d+', block.group(0), re.M)]

def evaluate(label, decls, final, tmp, n_pos, n_test, test_bp, extra=()):
    text = 'alphabet DNA;\n' + decls + f'\nsearch {final} in "pos.fa";\nsearch {final} in "test.fa";\n'
    res = run(text, tmp, extra)
    if 'REJECTED' in res.stdout:
        reason = [l for l in res.stderr.splitlines() if 'error' in l]
        return dict(label=label, rejected=reason[0] if reason else 'rejected')
    pos_hits = hits_by_record(res.stdout, final, 'pos.fa')
    ctl_hits = hits_by_record(res.stdout, final, 'test.fa')
    chosen = re.search(r'note: chose (.*)', res.stdout)
    recall = len(set(pos_hits)) / n_pos
    ctl_win = len(set(ctl_hits)) / n_test
    return dict(label=label, chosen=chosen.group(1) if chosen else '', recall=recall,
                ci=wilson(len(set(pos_hits)), n_pos), ctl_kb=1000 * len(ctl_hits) / test_bp,
                ctl_win=ctl_win, enrich=recall / ctl_win if ctl_win else float('inf'))

def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0, c - h), min(1, c + h))

VARIANTS = ['AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'CATAAA', 'GATAAA',
            'AATATA', 'AATACA', 'AATAGA', 'ACTAAA', 'AAGAAA', 'AATGAA']   # Beaudoing et al. 2000

IDX = {c: i for i, c in enumerate('ACGT')}

def pwm_from_sites(windows, consensus='AATAAA'):
    """Position weight matrix built the usual way: align training sites on
    their best match to the consensus (at most one mismatch), count, add
    pseudocounts. Windows with no such match contribute nothing."""
    import numpy as np
    w = len(consensus); counts = np.full((w, 4), 0.5)
    for _, seq in windows:
        best = None
        for i in range(len(seq) - w + 1):
            d = sum(1 for a, b in zip(seq[i:i+w], consensus) if a != b)
            if d <= 1 and (best is None or d < best[0]): best = (d, i)
        if best:
            for j, ch in enumerate(seq[best[1]:best[1]+w]): counts[j, IDX[ch]] += 1
    return counts / counts.sum(axis=1, keepdims=True)

def pwm_scores(pwm, windows, bg):
    """Log-odds score of every start position of every window."""
    import numpy as np
    w = pwm.shape[0]; lo = np.log2(pwm / bg)
    per_window = []
    for _, seq in windows:
        idx = np.array([IDX[c] for c in seq])
        sc = np.array([lo[np.arange(w), idx[i:i+w]].sum() for i in range(len(seq) - w + 1)])
        per_window.append(sc)
    return per_window

def pwm_eval(train_sites, train_ctrl, test_pos, test_ctrl, targets_per_kb):
    """PWM baseline: threshold chosen on the TRAINING controls so the hit
    rate is at most the target (like the compiler's contract), then scored
    on the test sites and held-out controls."""
    import numpy as np
    comp = np.zeros(4)
    for _, s in train_ctrl:
        for c in s: comp[IDX[c]] += 1
    bg = comp / comp.sum()
    pwm = pwm_from_sites(train_sites)
    tr = np.concatenate(pwm_scores(pwm, train_ctrl, bg)); tr_bp = sum(len(s) for _, s in train_ctrl)
    te_pos = pwm_scores(pwm, test_pos, bg)
    te_ctl = pwm_scores(pwm, test_ctrl, bg); te_bp = sum(len(s) for _, s in test_ctrl)
    te_ctl_all = np.concatenate(te_ctl)
    out = {}
    ordered = np.sort(tr)[::-1]
    for target in targets_per_kb:
        allowed = int(target * tr_bp / 1000)           # hits allowed on the training controls
        t = ordered[allowed] + 1e-9 if allowed < len(ordered) else ordered[-1]
        recall = sum(1 for sc in te_pos if (sc >= t).any()) / len(te_pos)
        out[target] = (recall, 1000 * (te_ctl_all >= t).sum() / te_bp,
                       sum(1 for sc in te_ctl if (sc >= t).any()) / len(te_ctl))
    return out

def split_eval(pos, train, test, splits, rng, tmp):
    """Fixed configurations, plus budgets the compiler learns from half of
    the known 3' ends (--sites) and the training controls (--empirical),
    scored on the other half and on held-out controls. Repeated over
    random splits of the positives."""
    fixed = [
        ('exact AATAAA',             'motif pas = "AATAAA";'),
        ('expert: AATAAA | ATTAAA',  'motif pas { "AATAAA" | "ATTAAA" }'),
        ('expert: 12 known variants', 'motif pas { ' + ' | '.join(f'"{v}"' for v in VARIANTS) + ' }'),
        ('AATAAA, 1 mismatch',       'motif pas = "AATAAA" allow 1 mismatch;'),
        ('AATAAA, 1 edit',           'motif pas = "AATAAA" allow 1 edits;'),
    ]
    learned = []
    for target in (2, 3, 4, 6, 8, 12, 16):
        learned.append((f'learned budget, <= {target}/kb',
                        f'motif pas = "AATAAA" allow ? mismatch require rate <= {target} per kb;'))
        learned.append((f'learned budget + positions, <= {target}/kb',
                        f'motif pas = "AATAAA" allow ? mismatch at ? require rate <= {target} per kb;'))
    results = {label: [] for label, _ in fixed + learned}
    pwm_targets = (2, 3, 4, 6, 8, 12, 16)
    pwm_results = {t: [] for t in pwm_targets}
    chosen = {label: {} for label, _ in learned}
    test_bp = sum(len(w) for _, w in test)
    write_records(os.path.join(tmp, 'train.fa'), train)
    write_records(os.path.join(tmp, 'test.fa'), test)
    for sp in range(splits):
        idx = list(range(len(pos))); rng.shuffle(idx)
        half = len(idx) // 2
        tr = [pos[i] for i in idx[:half]]; te = [pos[i] for i in idx[half:]]
        write_records(os.path.join(tmp, 'sites.fa'), tr)
        write_records(os.path.join(tmp, 'pos.fa'), te)
        for t, v in pwm_eval(tr, train, te, test, pwm_targets).items():
            pwm_results[t].append(v)
        for label, decl in fixed + learned:
            extra = ['--sites', 'sites.fa'] if label.startswith('learned') else []
            r = evaluate(label, decl, 'pas', tmp, len(te), len(test), test_bp, extra)
            if 'rejected' in r: continue
            results[label].append(r)
            if label in chosen:
                key = r['chosen'].replace('pas = ', '')
                chosen[label][key] = chosen[label].get(key, 0) + 1
    print(f'Held-out evaluation: {splits} random splits of the {len(pos)} 3\' ends into '
          f'{len(pos)//2} training sites (shown to the compiler) and {len(pos) - len(pos)//2} test sites\n')
    print(f'{"configuration":40s} {"test recall":>12s} {"ctrl/kb":>9s} {"ctrl hit":>9s} {"enrich":>7s}   most often chosen')
    def mean_sd(xs):
        m = sum(xs) / len(xs); v = sum((x - m) ** 2 for x in xs) / max(1, len(xs) - 1)
        return m, v ** 0.5
    rows = []
    for label, _ in fixed + learned:
        rs = results[label]
        if not rs:
            print(f'{label:40s} rejected in every split'); continue
        rec, rsd = mean_sd([r['recall'] for r in rs])
        kb, _ = mean_sd([r['ctl_kb'] for r in rs])
        win, _ = mean_sd([r['ctl_win'] for r in rs])
        top = ''
        if label in chosen and chosen[label]:
            k, v = max(chosen[label].items(), key=lambda kv: kv[1])
            top = f'{k} ({v}/{len(rs)})'
        rows.append((label, rec, rsd, kb, win, rec / win if win else float('inf'), top, len(rs)))
        print(f'{label:40s} {rec:6.3f}±{rsd:4.3f} {kb:9.2f} {win:9.3f} {rec / win if win else float("inf"):7.2f}   {top}')
    for t in pwm_targets:
        rs = pwm_results[t]
        rec, rsd = mean_sd([r[0] for r in rs]); kb, _ = mean_sd([r[1] for r in rs]); win, _ = mean_sd([r[2] for r in rs])
        label = f'PWM baseline, <= {t}/kb'
        rows.append((label, rec, rsd, kb, win, rec / win if win else float('inf'), '', len(rs)))
        print(f'{label:40s} {rec:6.3f}±{rsd:4.3f} {kb:9.2f} {win:9.3f} {rec / win if win else float("inf"):7.2f}')
    return rows, chosen

def main():
    seq = read_fasta(sys.argv[1]); tx = parse_gtf(sys.argv[2])
    per_site = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 7)
    splits = int(sys.argv[5]) if len(sys.argv) > 5 else 20
    pos, train, test = build_sets(seq, tx, UPSTREAM, per_site, rng)
    test_bp = sum(len(w) for _, w in test)
    print(f'E3: polyadenylation signals at annotated 3\' ends ({os.path.basename(sys.argv[2])})')
    print(f'{len(pos)} positive windows (60 nt upstream of each distinct 3\' end), '
          f'{len(train)} training and {len(test)} held-out control windows ({test_bp:,} bp)\n')
    with tempfile.TemporaryDirectory() as tmp:
        split_eval(pos, train, test, splits, rng, tmp)

if __name__ == '__main__':
    main()
