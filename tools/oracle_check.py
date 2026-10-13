#!/usr/bin/env python3
"""Brute-force oracle for MotifQL's approximate matching.

For random IUPAC motifs and random DNA, compares the set of hit START
positions reported by the compiled matcher against a direct dynamic-
programming check of the definition:

  mismatch model: position p is a hit if the window of motif length
                  starting at p differs in at most k positions.
  edit model:     position p is a hit if some non-empty substring
                  starting at p is within Levenshtein distance k of an
                  instance of the motif, with no residue inserted
                  before the first motif position (free-start search).

Usage: python3 tools/oracle_check.py [trials] [seed]
"""
import os, random, re, subprocess, sys, tempfile

IUPAC = {'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'CG','W':'AT',
         'K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}
BIN = os.path.join(os.path.dirname(__file__), '..', 'build', 'motifql')

def hamming_hits(pat, text, k):
    m = len(pat); out = set()
    for p in range(len(text) - m + 1):
        mism = sum(1 for i in range(m) if text[p+i] not in IUPAC[pat[i]])
        if mism <= k: out.add(p)
    return out

def edit_hits(pat, text, k):
    m = len(pat); INF = 10**9; out = set()
    for p in range(len(text)):
        L = min(len(text) - p, m + k)
        prev = [0] + [INF] * L                     # row 0: nothing inserted before pos 0
        for i in range(1, m + 1):
            cur = [i] + [0] * L
            allowed = IUPAC[pat[i-1]]
            for j in range(1, L + 1):
                cost = 0 if text[p+j-1] in allowed else 1
                cur[j] = min(prev[j-1] + cost, prev[j] + 1, cur[j-1] + 1)
            prev = cur
        if min(prev[1:L+1], default=INF) <= k: out.add(p)
    return out

def run_motifql(pat, k, model, text, tmp):
    with open(os.path.join(tmp, 'r.fa'), 'w') as f: f.write('>r\n' + text + '\n')
    word = 'edits' if model == 'edit' else 'mismatch'
    with open(os.path.join(tmp, 'r.motif'), 'w') as f:
        f.write(f'alphabet DNA;\nmotif m = "{pat}" allow {k} {word};\nsearch m in "r.fa";\n')
    res = subprocess.run([os.path.abspath(BIN), '--no-search', 'r.motif'], cwd=tmp, capture_output=True, text=True)
    if 'ACCEPTED' not in res.stdout: return None
    res = subprocess.run([os.path.abspath(BIN), '--symbols', 'r.motif'], cwd=tmp, capture_output=True, text=True)
    return {int(s) - 1 for s in re.findall(r'^\s+r:(\d+)-\d+', res.stdout, re.M)}

def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    codes = 'ACGTACGTACGTRYSWKMN'
    ok = bad = skipped = 0
    hits_compared = {'mismatch': 0, 'edit': 0}; cases = {'mismatch': 0, 'edit': 0}
    with tempfile.TemporaryDirectory() as tmp:
        for t in range(trials):
            m = rng.randint(5, 9); k = rng.randint(1, 2); model = rng.choice(['mismatch', 'edit'])
            pat = ''.join(rng.choice(codes) for _ in range(m))
            text = ''.join(rng.choice('ACGT') for _ in range(rng.randint(200, 600)))
            got = run_motifql(pat, k, model, text, tmp)
            if got is None: skipped += 1; continue   # rejected by semantic checks (budget too loose)
            want = edit_hits(pat, text, k) if model == 'edit' else hamming_hits(pat, text, k)
            cases[model] += 1; hits_compared[model] += len(want)
            if got == want: ok += 1
            else:
                bad += 1
                print(f'MISMATCH {model} k={k} {pat}: only-motifql={sorted(got-want)[:5]} only-oracle={sorted(want-got)[:5]}')
    for mdl in cases: print(f'  {mdl:8s}: {cases[mdl]} cases, {hits_compared[mdl]} oracle hit positions')
    print(f'oracle agreement: {ok} agree, {bad} disagree, {skipped} rejected at compile time')
    sys.exit(1 if bad else 0)

if __name__ == '__main__':
    main()
