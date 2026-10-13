#!/usr/bin/env python3
"""Checks the compile-time chance-hit estimate against reality.

Generates random sequences, runs real searches with the compiled
matcher, and compares the number of hits actually found with the
number the specificity pass predicts from the automaton alone.

  Experiment A: uniform random DNA, a range of motifs and budgets.
  Experiment B: CpG-depleted DNA (an order-1 Markov source), to show
                why the background model's order matters.

Usage: python3 tools/rate_check.py [megabases] [seed]
"""
import os, random, re, subprocess, sys, tempfile

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'build', 'motifql'))

MOTIFS_A = [
    ('ecori',      '"GAATTC"', ''),
    ('ecori_mm1',  '"GAATTC"', 'allow 1 mismatch'),
    ('ecori_ed1',  '"GAATTC"', 'allow 1 edits'),
    ('polya_ed1',  '"AATAAA"', 'allow 1 edits'),
    ('tata',       '"TATAWAWR"', ''),
    ('tata_mm1',   '"TATAWAWR"', 'allow 1 mismatch'),
    ('tata_ed1',   '"TATAWAWR"', 'allow 1 edits'),
    ('tata_ed2',   '"TATAWAWR"', 'allow 2 edits'),
    ('ebox_sp',    '{ "CANNTG" "N"{0,5} "CANNTG" }', ''),
]
MOTIFS_B = [
    ('hpaii',      '"CCGG"', ''),
    ('gcgc_mm1',   '"GCGCGC"', 'allow 1 mismatch'),
    ('tata_ed1',   '"TATAWAWR"', 'allow 1 edits'),
]

def write_motif(path, motifs, fasta):
    with open(path, 'w') as f:
        f.write('alphabet DNA;\n')
        for name, body, tol in motifs:
            if body.startswith('{'):
                f.write(f'motif {name} {body} {tol};\n')
            else:
                f.write(f'motif {name} = {body} {tol};\n')
        for name, _, _ in motifs:
            f.write(f'search {name} in "{fasta}";\n')

def uniform_seq(n, rng):
    return ''.join(rng.choice('ACGT') for _ in range(n))

def cpg_depleted_seq(n, rng):
    # Order-1 source: after C, G is five times rarer than the other bases.
    out = [rng.choice('ACGT')]
    for _ in range(n - 1):
        if out[-1] == 'C':
            out.append(rng.choices('ACGT', weights=[5, 5, 1, 5])[0])
        else:
            out.append(rng.choice('ACGT'))
    return ''.join(out)

def search_counts(motif_file, cwd):
    res = subprocess.run([BIN, '--no-search', motif_file], cwd=cwd, capture_output=True, text=True)
    res = subprocess.run([BIN, '--symbols', motif_file], cwd=cwd, capture_output=True, text=True)
    blocks = re.findall(r'^(\w+) in "[^"]+".*?-> (\d+) hit\(s\)\n  expected by chance: ([\d.]+)', res.stdout, re.S | re.M)
    return {name: (int(obs), float(exp)) for name, obs, exp in blocks}

def spec_table(motif_file, cwd, bg=None, order=None):
    cmd = [BIN, '--specificity', '--no-search', motif_file]
    if bg: cmd += ['--background', bg, '--markov', str(order)]
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    table = res.stdout.split('== specificity')[1]
    out = {}
    for line in table.splitlines()[2:]:
        cols = line.split()
        if len(cols) < 5 or cols[0] == 'motifs': continue
        # tolerance is one word ("exact") or two ("1 edit")
        q = cols[3] if cols[1] == 'exact' else cols[4]
        out[cols[0]] = float(q)
    return out

def main():
    mb = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
    n = int(mb * 1_000_000)
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, 'u.fa'), 'w') as f: f.write('>uniform\n' + uniform_seq(n, rng) + '\n')
        write_motif(os.path.join(tmp, 'a.motif'), MOTIFS_A, 'u.fa')
        print(f'Experiment A: {mb:g} Mb uniform random DNA')
        print(f'{"motif":12s} {"predicted":>10s} {"observed":>9s} {"obs/pred":>9s}')
        for name, (obs, exp) in search_counts('a.motif', tmp).items():
            print(f'{name:12s} {exp:10.1f} {obs:9d} {obs/exp if exp else float("nan"):9.3f}')

        with open(os.path.join(tmp, 'c.fa'), 'w') as f: f.write('>cpg_depleted\n' + cpg_depleted_seq(n, rng) + '\n')
        write_motif(os.path.join(tmp, 'b.motif'), MOTIFS_B, 'c.fa')
        obs = search_counts('b.motif', tmp)
        q0 = spec_table('b.motif', tmp, 'c.fa', 0)
        q1 = spec_table('b.motif', tmp, 'c.fa', 1)
        print(f'\nExperiment B: {mb:g} Mb CpG-depleted DNA (order-1 source)')
        print(f'{"motif":12s} {"order-0 pred":>12s} {"order-1 pred":>12s} {"observed":>9s}')
        for name in obs:
            print(f'{name:12s} {q0[name]*n:12.1f} {q1[name]*n:12.1f} {obs[name][0]:9d}')

if __name__ == '__main__':
    main()
