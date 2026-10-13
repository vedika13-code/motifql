#!/usr/bin/env python3
"""Evaluation on real human sequence (hg38 chromosome 21).

Usage: python3 tools/eval_chr21.py path/to/chr21.fa [outdir]

E1  Specificity on real DNA. For a panel of well-known regulatory and
    restriction-site consensus motifs, compares the chance-hit count the
    compiler predicts (order-0 and order-1 backgrounds fitted to chr21)
    with (a) the count observed on a synthetic sequence drawn from the
    order-1 model itself, which tests the method, and (b) the count
    observed on real chr21, which also contains biology: repeats,
    CpG depletion, real sites.

E2  Budget synthesis on planted sites. Plants mutated copies of a
    composite promoter motif (TATA box, 20-40 bp spacer, Kozak core) into
    a 5 Mb stretch of real chr21, then compares exact matching,
    hand-picked budgets, and budgets synthesised by the compiler from a
    `require rate <= ...` clause: how many planted sites each finds, and
    how many hits fall outside planted sites.
"""
import os, random, re, subprocess, sys, tempfile, time

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'build', 'motifql'))
IUPAC = {'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'CG','W':'AT',
         'K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}

PANEL = [  # name, consensus, tolerance
    ('TATA_box',  'TATAWAWR',   ''),
    ('TATA_box',  'TATAWAWR',   'allow 1 edits'),
    ('CAAT_box',  'CCAAT',      ''),
    ('E_box',     'CACGTG',     ''),
    ('E_box',     'CACGTG',     'allow 1 edits'),
    ('AP1',       'TGASTCA',    ''),
    ('CRE',       'TGACGTCA',   ''),
    ('CRE',       'TGACGTCA',   'allow 1 mismatch'),
    ('NFkB',      'GGGRNTTYCC', ''),
    ('GC_box',    'GGGGCGGGG',  ''),
    ('polyA',     'AATAAA',     ''),
    ('polyA',     'AATAAA',     'allow 1 edits'),
    ('Kozak',     'GCCRCCATGG', ''),
    ('Kozak',     'GCCRCCATGG', 'allow 1 edits'),
    ('p53_half',  'RRRCWWGYYY', ''),
    ('CTCF_core', 'CCGCGNGGNGGCAG', ''),
    ('CTCF_core', 'CCGCGNGGNGGCAG', 'allow 1 edits'),
    ('EcoRI',     'GAATTC',     ''),
]

def read_fasta(path):
    with open(path) as f:
        f.readline()
        return f.read().replace('\n', '').upper()

def write_fasta(path, name, seq):
    with open(path, 'w') as f:
        f.write(f'>{name}\n')
        for i in range(0, len(seq), 80): f.write(seq[i:i+80] + '\n')

def order1_model(seq):
    counts = {a: {b: 1 for b in 'ACGT'} for a in 'ACGT'}
    prev = None
    for c in seq:
        if c not in 'ACGT': prev = None; continue
        if prev: counts[prev][c] += 1
        prev = c
    return {a: [counts[a][b] / sum(counts[a].values()) for b in 'ACGT'] for a in 'ACGT'}

def markov_surrogate(model, n, rng):
    out = [rng.choice('ACGT')]
    cum = {a: [sum(model[a][:i+1]) for i in range(4)] for a in 'ACGT'}
    for _ in range(n - 1):
        u = rng.random(); c = cum[out[-1]]
        out.append('A' if u < c[0] else 'C' if u < c[1] else 'G' if u < c[2] else 'T')
    return ''.join(out)

def run(args, cwd):
    return subprocess.run([BIN] + args, cwd=cwd, capture_output=True, text=True).stdout

def panel_motif_file(path, fasta):
    with open(path, 'w') as f:
        f.write('alphabet DNA;\n')
        for i, (name, pat, tol) in enumerate(PANEL):
            f.write(f'motif m{i} = "{pat}" {tol};\n')
        for i in range(len(PANEL)):
            f.write(f'search m{i} in "{fasta}";\n')

def observed_counts(out):
    return {int(m.group(1)): int(m.group(2))
            for m in re.finditer(r'^m(\d+) in "[^"]+".*?-> (\d+) hit\(s\)', out, re.S | re.M)}

def spec_probs(out):
    probs = {}
    for line in out.split('== specificity')[1].splitlines()[2:]:
        cols = line.split()
        if not cols or not cols[0].startswith('m') or cols[0] == 'motifs': continue
        probs[int(cols[0][1:])] = float(cols[3] if cols[1] == 'exact' else cols[4])
    return probs

def e1(chr21, tmp, rng):
    print('E1: chance-hit prediction on chr21 (hg38)\n')
    # Out-of-sample: backgrounds are fitted on chr21 up to 28 Mb and
    # predictions are checked on chr21 from 28 Mb to the end, so a
    # high-order background cannot simply memorise the counts it is
    # judged against.
    train, test = chr21[:SPLIT], chr21[SPLIT:]
    write_fasta(os.path.join(tmp, 'train.fa'), 'chr21_train', train)
    write_fasta(os.path.join(tmp, 'chr21.fa'), 'chr21_test', test)
    valid = sum(test.count(c) for c in 'ACGT')
    print(f'background fitted on chr21:1-{SPLIT:,} ({sum(train.count(c) for c in "ACGT"):,} non-N positions); '
          f'predictions checked on chr21:{SPLIT+1:,}-end ({valid:,} non-N positions)\n')
    model = order1_model(chr21[:SPLIT])
    surrogate = markov_surrogate(model, 10_000_000, rng)
    write_fasta(os.path.join(tmp, 'sur.fa'), 'order1_surrogate', surrogate)

    panel_motif_file(os.path.join(tmp, 'panel.motif'), 'chr21.fa')
    orders = [0, 1, 3, 5]
    q = {k: spec_probs(run(['--specificity', '--no-search', '--background', 'train.fa', '--markov', str(k), 'panel.motif'], tmp))
         for k in orders}
    q['emp'] = spec_probs(run(['--specificity', '--no-search', '--empirical', 'train.fa', 'panel.motif'], tmp))
    orders = orders + ['emp']
    real = observed_counts(run(['--symbols', '--count-only', 'panel.motif'], tmp))
    panel_motif_file(os.path.join(tmp, 'panel_sur.motif'), 'sur.fa')
    sur = observed_counts(run(['--symbols', '--count-only', '--markov', '1', 'panel_sur.motif'], tmp))

    rows = []
    print('Part 1: order-1 surrogate of chr21 (tests the method)')
    print(f'{"motif":10s} {"pattern":15s} {"tolerance":17s} {"predicted":>9s} {"observed":>9s} {"obs/pred":>8s}')
    for i, (name, pat, tol) in enumerate(PANEL):
        sp = q[1][i] * len(surrogate)
        print(f'{name:10s} {pat:15s} {tol or "exact":17s} {sp:9.0f} {sur[i]:9d} {sur[i]/sp if sp else float("nan"):8.2f}')
    print('\nPart 2: real chr21 (method plus biology), observed / predicted for each background order')
    print(f'{"motif":10s} {"pattern":15s} {"tolerance":17s} {"observed":>9s} ' + ' '.join(f'{("o"+str(k)) if k != "emp" else "emp":>6s}' for k in orders))
    import math
    logerr = {k: [] for k in orders}
    for i, (name, pat, tol) in enumerate(PANEL):
        ratios = [real[i] / (q[k][i] * valid) if q[k][i] > 0 else float('nan') for k in orders]
        if real[i] >= 20:
            for k, r in zip(orders, ratios):
                if r == r and r > 0: logerr[k].append(abs(math.log2(r)))
        rows.append((name, pat, tol or 'exact', q[1][i] * len(surrogate), sur[i], real[i], ratios))
        print(f'{name:10s} {pat:15s} {tol or "exact":17s} {real[i]:9d} ' + ' '.join(f'{r:6.2f}' for r in ratios))
    print('\nmedian |log2(observed/predicted)| over motifs with >= 20 hits: ' +
          ', '.join(f'{("order " + str(k)) if k != "emp" else "empirical"}: {sorted(v)[len(v)//2]:.2f}' for k, v in logerr.items()))
    print(f'(surrogate: {len(surrogate):,} positions drawn from the order-1 model of the training half)\n')
    return rows

SPLIT = 28_000_000
PROMOTER_TATA, PROMOTER_KOZAK = 'TATAWAWR', 'RCCATGG'

def mutate(s, sub, indel, rng):
    out = []
    for c in s:
        u = rng.random()
        if u < indel / 2: continue
        if u < indel: out.append(rng.choice('ACGT'))
        if rng.random() < sub: c = rng.choice([x for x in 'ACGT' if x != c])
        out.append(c)
    return ''.join(out)

def e2(chr21, tmp, rng):
    print('E2: budget synthesis on planted promoters in real chr21 sequence\n')
    region = list(chr21[32_000_000:37_000_000])
    assert 'N' not in region
    write_fasta(os.path.join(tmp, 'bg.fa'), 'chr21_train', chr21[:SPLIT])   # background: the other half of chr21
    sites = []
    for p in range(10_000, len(region) - 1_000, 20_000):
        tata = ''.join(rng.choice(IUPAC[c]) for c in PROMOTER_TATA)
        kozak = ''.join(rng.choice(IUPAC[c]) for c in PROMOTER_KOZAK)
        spacer = ''.join(rng.choice('ACGT') for _ in range(rng.randint(20, 40)))
        site = mutate(tata + spacer + kozak, 0.10, 0.02, rng)
        region[p:p+len(site)] = list(site)        # replace bases; coordinates stay fixed
        sites.append((p, p + len(site)))
    seq = ''.join(region)
    write_fasta(os.path.join(tmp, 'planted.fa'), 'planted', seq)

    configs = [
        ('exact (0, 0)',           'allow 0 edits', 'allow 0 edits', ''),
        ('hand-picked (1, 1)',     'allow 1 edits', 'allow 1 edits', ''),
        ('hand-picked (2, 2)',     'allow 2 edits', 'allow 2 edits', ''),
        ('synthesised, <= 2/100kb',  'allow ? edits', 'allow ? edits', 'require rate <= 2 per 100kb, divergence 10%'),
        ('synthesised, <= 20/100kb', 'allow ? edits', 'allow ? edits', 'require rate <= 20 per 100kb, divergence 10%'),
    ]
    results = {}
    for mode in ('order 1', 'order 5', 'empirical'):
        results[mode] = e2_run(configs, mode, sites, seq, tmp)
    print(f'{len(sites)} planted sites (10% substitution, 2% indel) in chr21:32,000,001-37,000,000; '
          f'background fitted on chr21:1-{SPLIT:,}\n')
    for order in results:
        print(f'chance-hit rates from: {order} (fitted on / scanned over chr21:1-{SPLIT:,})')
        print(f'{"configuration":26s} {"budgets chosen":32s} {"recall":>6s} {"pred":>5s} {"false/100kb":>11s} {"pred":>6s}')
        for label, ch, rec, ps, fr, pr, dt in results[order]:
            print(f'{label:26s} {ch:32s} {rec:6.3f} {ps:5.2f} {fr:11.2f} {pr:6.2f}   ({dt:.1f} s)')
        print()
    return results

def e2_run(configs, order, sites, seq, tmp):
    rows = []
    for label, ta, ko, req in configs:
        ta_s = '' if ta == 'allow 0 edits' else ta
        ko_s = '' if ko == 'allow 0 edits' else ko
        with open(os.path.join(tmp, 'p.motif'), 'w') as f:
            f.write(f'alphabet DNA;\nmotif tata = "{PROMOTER_TATA}" {ta_s};\nmotif kozak = "{PROMOTER_KOZAK}" {ko_s};\n'
                    f'motif promoter {{ @tata "N"{{20,40}} @kozak }} {req};\nsearch promoter in "planted.fa";\n')
        t0 = time.time()
        flags = ['--empirical', 'bg.fa'] if order == 'empirical' else ['--background', 'bg.fa', '--markov', order.split()[1]]
        out = run(['--symbols', '--sensitivity'] + flags + ['p.motif'], tmp)
        dt = time.time() - t0
        chosen = re.search(r'note: chose (.*)', out)
        hits = [int(h) - 1 for h in re.findall(r'^\s+planted:(\d+)-\d+', out, re.M)]
        ends = {}
        found = set(); false = 0
        j = 0
        hits.sort()
        for h in hits:
            while j < len(sites) and sites[j][1] <= h - 120: j += 1
            hit_site = None
            for k in range(max(0, j - 1), min(len(sites), j + 2)):
                s0, s1 = sites[k]
                if s0 - 60 <= h < s1: hit_site = k   # a hit starting in or just before the site
            if hit_site is None: false += 1
            else: found.add(hit_site)
        sens = re.search(r'^promoter\s+\S+(?: \S+)?\s+(\S+)\s+(\S+)', out.split('== sensitivity')[1], re.M)
        pred_rate = float(sens.group(1)) / 10 if sens else float('nan')   # per 100 kb
        pred_sens = float(sens.group(2)) if sens else float('nan')
        rows.append((label, chosen.group(1) if chosen else '-', len(found) / len(sites), pred_sens,
                     false / (len(seq) / 1e5), pred_rate, dt))
    return rows

def main():
    path = sys.argv[1]
    rng = random.Random(2026)
    chr21 = read_fasta(path)
    with tempfile.TemporaryDirectory() as tmp:
        e1(chr21, tmp, rng)
        e2(chr21, tmp, rng)

if __name__ == '__main__':
    main()
