# Phase 3: Sensitivity, Per-Part Budgets and Budget Synthesis

Phase 3 completes the research direction: the compiler now reports
both sides of an approximate motif, how many chance hits it costs and
how many diverged real sites it finds, and can choose error budgets
itself from a stated target.

## 1. Per-part budgets (`nfa.c`, `sema.c`)

A referenced motif with its own `allow k ...` is built and layered on
its own, then embedded in the composite with every accepting state
turned into an exit (`frag_embed`). Length inference widens a
reference by its edit budget so the matcher's length cut-off does not
stop inside an insertion.

The brute-force oracle (`tools/oracle_check.py`) found one real
question of definition here: may a part absorb an inserted residue
just after its last position? Inside a composite that changes which
sites match. The convention chosen is that edits count only inside a
part's span; the spacer absorbs residues between parts. With the
oracle using the same definition, `oracle_check.py 2000 2026` gives
exact agreement on all 2,000 single-motif cases (75,287 hit positions)
and all 470 composite cases with per-part budgets (600 hit positions).

## 2. Sensitivity (`sensitivity.c`)

Monte Carlo with a fixed seed: draw an exact instance of the motif
from its syntax tree, mutate it (substitution rate s, indel rate g),
and count it found if the compiled matcher reports any hit in it.
Reported with a 95% Wilson interval. With indels off, the estimate for
`TATAWAWR` matches the closed form 0.9^5 x (0.9 + 0.1/3)^3 = 0.4801
(test suite). Every candidate budget is scored on the same random
sites, so comparisons are paired and monotone.

## 3. Budget synthesis (`synthesis.c`)

For a `require` clause, the free budgets in the motif and the motifs it
references are enumerated (each 0..cap, cap <= 3 and leaving at least
half the motif unedited), smallest total first. Rate and sensitivity
both only grow with any budget, so:

- an assignment at or above a known infeasible one is never compiled
  (monotone pruning: 10 of 16 assignments skipped for the promoter);
- the best feasible assignment is maximal, so sensitivity is only
  estimated for maximal feasible assignments and for the minimal
  infeasible ones, which are shown as "what one more error costs".

## 4. Chance-hit rates on real DNA

Backgrounds are now Markov chains of any order up to 6
(`--markov K`), and `--empirical FILE` measures rates by scanning a
reference sequence instead. See `tools/eval_chr21.py`.

E1, out of sample (background from chr21:1-28 Mb, counts on
chr21:28 Mb-end, hg38):

- On sequence drawn from the order-1 model itself, predictions match
  observed counts to within Poisson noise for every motif with more
  than 100 hits. The method is exact for the model it assumes.
- On real chr21, median |log2(observed/predicted)|: order 0 0.43,
  order 1 0.84, order 3 0.46, order 5 0.59, empirical 0.57. The two
  halves of the chromosome differ in composition, which sets a floor
  of about 1.35x on any global estimate. CpG-rich motifs show the
  largest errors (GC box: 49x under order 0, 3.1x empirical).

E2, planted promoters (250 sites, 10% substitution, 2% indel) in
chr21:32-37 Mb:

| Configuration | Budgets | Recall | False hits / 100 kb |
| --- | --- | --- | --- |
| exact | 0, 0 | 0.20 | 0.04 |
| hand-picked | 1, 1 | 0.73 | 80.1 |
| hand-picked | 2, 2 | 0.99 | 4,666 |
| synthesised, <= 20 / 100 kb, empirical rates | tata 1, kozak 0 | 0.42 | 2.3 |
| synthesised, <= 2 / 100 kb, order-1 rates | tata 0, kozak 1 | 0.37 | 4.0 (target missed) |

Synthesis doubles recall over exact matching while producing 35 times
fewer false hits than the obvious hand-picked budgets. With an order-1
background the 2-per-100-kb target was missed by 2x on real sequence;
with order-5 or empirical rates every target was met. The guarantee is
only as good as the rate estimate, and on real genomes the empirical
or high-order estimate should be preferred.

## 5. Paper extensions

- **Position-restricted errors** (`allow k ... at 1,2,5`): error edges
  are added only at the listed residue positions. Checked by the oracle
  (500 random cases, exact agreement).
- **Learning from known sites** (`allow ? mismatch at ?` with
  `--sites FILE`): synthesis searches budgets and position sets, scoring
  sensitivity as the share of known sites found. The objective is
  lexicographic: highest sensitivity first, then fewest chance hits.
- **E3, polyadenylation signals** (`tools/eval_pas.py`): 203 annotated
  3' ends on chr21 (GENCODE subset), 20 random train/test splits, chance
  rates measured on training control windows and checked on held-out
  ones. Learned position sets trace a recall/false-hit curve close to
  the literature-curated PAS variant lists (Beaudoing et al., 2000), and
  the first position learned is position 2 (AATAAA to ATTAAA), the most
  common human variant. Full table in `results/eval_pas.txt`.
