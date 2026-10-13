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

def positional_hits(pat, k, model, mask, text):
    """Hit starts when errors are only allowed at the 0-based positions in `mask`.
    Edit model: a substitution or deletion of position i, or an insertion just
    before position i, costs 1 only if i is allowed (never before position 0)."""
    m = len(pat); out = set(); INF = 10**9
    for p in range(len(text)):
        if model == 'mismatch':
            if p + m > len(text): break
            ok, cost = True, 0
            for i in range(m):
                if text[p+i] not in IUPAC[pat[i]]:
                    if i in mask: cost += 1
                    else: ok = False; break
            if ok and cost <= k: out.add(p)
            continue
        L = min(len(text) - p, m + k)
        prev = [0] + [INF] * L
        for i in range(1, m + 1):
            allowed = IUPAC[pat[i-1]]; err = (i - 1) in mask
            cur = [prev[0] + 1 if err else INF] + [INF] * L
            for j in range(1, L + 1):
                match = prev[j-1] + (0 if text[p+j-1] in allowed else (1 if err else INF))
                dele = prev[j] + 1 if err else INF
                ins = cur[j-1] + 1 if (i < m and i in mask) else INF   # insertion before position i (0-based), i >= 1
                cur[j] = min(match, dele, ins)
            prev = cur
        if min(prev[1:L+1], default=INF) <= k: out.add(p)
    return out

def run_positional(pat, k, model, mask, text, tmp):
    with open(os.path.join(tmp, 'r.fa'), 'w') as f: f.write('>r\n' + text + '\n')
    word = 'edits' if model == 'edit' else 'mismatch'
    at = ','.join(str(i + 1) for i in sorted(mask))
    with open(os.path.join(tmp, 'r.motif'), 'w') as f:
        f.write(f'alphabet DNA;\nmotif m = "{pat}" allow {k} {word} at {at};\nsearch m in "r.fa";\n')
    res = subprocess.run([os.path.abspath(BIN), '--symbols', 'r.motif'], cwd=tmp, capture_output=True, text=True)
    if 'ACCEPTED' not in res.stdout: return None
    return {int(s) - 1 for s in re.findall(r'^\s+r:(\d+)-\d+', res.stdout, re.M)}

def part_ends(pat, k, model, text, q):
    """End offsets e (exclusive) such that text[q:e] matches pat within budget k.
    Edits count only inside the part's own span: no residue inserted before
    its first position or after its last (the same convention the compiler
    uses for a referenced motif with its own budget)."""
    m = len(pat); n = len(text); ends = set()
    if model == 'mismatch' or k == 0:
        if q + m <= n and sum(1 for i in range(m) if text[q+i] not in IUPAC[pat[i]]) <= k:
            ends.add(q + m)
        return ends
    INF = 10**9; L = min(n - q, m + k)
    if L <= 0: return ends
    prev = [0] + [INF] * L
    last = None
    for i in range(1, m + 1):
        cur = [i] + [0] * L
        nohoriz = [i] + [0] * L        # best cost where the last step is not an insertion
        allowed = IUPAC[pat[i-1]]
        for j in range(1, L + 1):
            cost = 0 if text[q+j-1] in allowed else 1
            nohoriz[j] = min(prev[j-1] + cost, prev[j] + 1)
            cur[j] = min(nohoriz[j], cur[j-1] + 1)
        prev = cur; last = nohoriz
    for j in range(1, L + 1):
        if last[j] <= k: ends.add(q + j)
    return ends

def composite_hits(a, ka, b, kb, model, lo, hi, text):
    out = set()
    for p in range(len(text)):
        found = False
        for e1 in part_ends(a, ka, model, text, p):
            for gap in range(lo, hi + 1):
                if part_ends(b, kb, model, text, e1 + gap):
                    found = True; break
            if found: break
        if found: out.add(p)
    return out

def run_composite(a, ka, b, kb, model, lo, hi, text, tmp):
    with open(os.path.join(tmp, 'r.fa'), 'w') as f: f.write('>r\n' + text + '\n')
    word = 'edits' if model == 'edit' else 'mismatch'
    tol = lambda k: f' allow {k} {word}' if k else ''
    with open(os.path.join(tmp, 'r.motif'), 'w') as f:
        f.write(f'alphabet DNA;\nmotif a = "{a}"{tol(ka)};\nmotif b = "{b}"{tol(kb)};\n'
                f'motif m {{ @a "N"{{{lo},{hi}}} @b }}\nsearch m in "r.fa";\n')
    res = subprocess.run([os.path.abspath(BIN), '--symbols', 'r.motif'], cwd=tmp, capture_output=True, text=True)
    if 'ACCEPTED' not in res.stdout: return None
    return {int(s) - 1 for s in re.findall(r'^\s+r:(\d+)-\d+', res.stdout, re.M)}

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
    # Composite motifs with per-part budgets: @a N{lo,hi} @b
    cok = cbad = 0; chits = 0
    with tempfile.TemporaryDirectory() as tmp:
        for t in range(max(1, trials // 4)):
            model = rng.choice(['mismatch', 'edit'])
            a = ''.join(rng.choice(codes) for _ in range(rng.randint(4, 6)))
            b = ''.join(rng.choice(codes) for _ in range(rng.randint(4, 6)))
            ka, kb = rng.randint(0, 1), rng.randint(0, 1)
            lo = rng.randint(0, 3); hi = lo + rng.randint(0, 3)
            text = ''.join(rng.choice('ACGT') for _ in range(rng.randint(150, 300)))
            got = run_composite(a, ka, b, kb, model, lo, hi, text, tmp)
            if got is None: continue
            want = composite_hits(a, ka, b, kb, model, lo, hi, text); chits += len(want)
            if got == want: cok += 1
            else:
                cbad += 1
                print(f'COMPOSITE MISMATCH {model} @{a}/{ka} N{{{lo},{hi}}} @{b}/{kb}: only-motifql={sorted(got-want)[:5]} only-oracle={sorted(want-got)[:5]}'); print('   text', text)
    print(f'composite per-part budgets: {cok} agree, {cbad} disagree ({chits} oracle hit positions)')
    bad += cbad
    # Position-restricted errors: allow k ... at P
    pok = pbad = 0; phits = 0
    with tempfile.TemporaryDirectory() as tmp:
        for t in range(max(1, trials // 4)):
            model = rng.choice(['mismatch', 'edit'])
            pat = ''.join(rng.choice(codes) for _ in range(rng.randint(5, 8)))
            k = rng.randint(1, 2)
            mask = set(rng.sample(range(len(pat)), rng.randint(1, len(pat) - 1)))
            text = ''.join(rng.choice('ACGT') for _ in range(rng.randint(150, 300)))
            got = run_positional(pat, k, model, mask, text, tmp)
            if got is None: continue
            want = positional_hits(pat, k, model, mask, text); phits += len(want)
            if got == want: pok += 1
            else:
                pbad += 1
                print(f'POSITIONAL MISMATCH {model} k={k} {pat} at {sorted(mask)}: only-motifql={sorted(got-want)[:5]} only-oracle={sorted(want-got)[:5]}')
    print(f'position-restricted errors: {pok} agree, {pbad} disagree ({phits} oracle hit positions)')
    bad += pbad
    sys.exit(1 if bad else 0)

if __name__ == '__main__':
    main()
