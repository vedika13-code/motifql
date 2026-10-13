# Phase 2: Automaton Construction, Optimization, and Matching

Phase 2 takes a validated motif (Phase 1's output: a type-checked AST
with resolved references and inferred lengths) and turns it into a
minimized DFA, then a bytecode matcher that scans FASTA sequences.
Every stage below is a distinct, inspectable pass: `--nfa`, `--dfa`,
and `--bytecode` dump their respective stage's output.

## 1. Pipeline

```
AST (Phase 1) --Thompson's construction--> NFA
              --(if mismatch budget > 0)--> layered NFA
              --subset construction------> DFA
              --dead-state elimination---> DFA (pruned)
              --partition-refinement------> minimized DFA
              --prefix extraction---------> literal prefix
              --packaging------------------> MatcherProgram (bytecode)
              --scan--------------------> hits over a FASTA record
```

Source files: `nfa.c`, `dfa.c`, `codegen.c`, `matcher.c`, `fasta.c`.

## 2. Thompson's construction (`nfa.c`)

Standard three-state-kind construction: `SPLIT` (epsilon branch),
`CONSUME` (one residue, guarded by a mask), `MATCH` (accept, no
outgoing edges). Every AST node type reduces to wiring these
together:

| AST node | NFA fragment |
|---|---|
| `NODE_LITERAL` | a chain of `CONSUME` states, one per residue |
| `NODE_CLASS` | a single `CONSUME` with the resolved class mask |
| `NODE_CONCAT` | fragments wired end to end |
| `NODE_ALTERNATE` | a `SPLIT` branching into each alternative |
| `NODE_REPEAT` | see below |
| `NODE_REF` | the referenced motif's fragment, built fresh |

`@reference` expansion is inlining, not a call/return mechanism: the
referenced motif's fragment is built again from its AST every time it
is referenced. This is safe only because Phase 1 already proved the
reference graph acyclic -- an unchecked cyclic reference would make
this recursion infinite, which is precisely the failure mode Phase
1's cycle detection exists to rule out before Phase 2 ever runs.

### Bounded repetition

`{min,max}` is unrolled rather than looped: `min` mandatory fresh
copies of the child fragment, concatenated, followed either by one
more fresh copy wrapped in a Kleene star (if `max` is unbounded) or by
`(max - min)` fresh copies each wrapped in a `?` (zero-or-one),
chained in series. Each copy gets its own states -- Thompson's
construction never shares states across repetitions, which is what
keeps the construction correct without needing a separate counter
mechanism.

## 3. Mismatch tolerance (`nfa_apply_mismatch_tolerance`)

`allow k mismatch` is compiled as a **layered automaton**: the plain
NFA is replicated into `k+1` layers, 0 through k, where layer *i*
means "*i* substitutions used so far". A `CONSUME` state with mask
`M` and target `t`, at layer *i*, becomes two edges:

- on any residue **in** `M` &rarr; layer *i*, state `t` (no cost)
- on any residue **not in** `M` &rarr; layer *i+1*, state `t`
  (substitution, only if *i* < k)

`SPLIT` states are simply replicated per layer with same-layer
targets, since they consume nothing. `MATCH` is accepting at every
layer: reaching the end having used at most k substitutions is
success regardless of which layer got there.

This is Hamming-distance approximate matching (substitutions only),
not general edit distance -- there is no insertion or deletion, which
matches what "mismatch" means for a fixed-frame biological motif and
is considerably simpler to construct than a full k-differences
automaton (no epsilon layers for indels are needed).

One pleasant confirmation this construction is right: the complement
mask on a mismatch edge, printed through `alphabet_format_mask`,
comes out as a real IUPAC code. The mismatch edge for a `CONSUME [A]`
state prints as `CONSUME [B]` -- and `B` is the actual IUPAC code for
"not A". That was not hand-coded; it fell out of reusing the same
mask-formatting logic Phase 1 already had for printing ambiguity
codes.

## 4. Subset construction (`dfa_from_nfa`)

Standard construction: a DFA state is an epsilon-closed set of NFA
states, discovered by a worklist starting from the closure of the
NFA's start state. For each discovered set and each concrete alphabet
residue, `move` collects every `CONSUME` state's target whose mask
contains that residue, and the epsilon-closure of that becomes the
next DFA state (constructing it if new). A residue with no live
targets produces `DFA_DEAD` directly -- no trap state is ever
physically created for the ordinary case of "no transition here".

## 5. Dead-state elimination -- and an honest negative result

`dfa_eliminate_dead_states` computes, by backward fixed-point
propagation from the accepting states, which states can reach
acceptance at all, and prunes anything that cannot.

Across every test motif in this project -- fixed literals, classes,
alternation, bounded and unbounded repetition, composed `@reference`
patterns, and the mismatch-tolerance layering -- **this pass removes
zero states, every time.** That is not a bug in the pass; it is a
genuine property of Thompson's construction. Every dangling pointer
in the construction is patched toward something that eventually
leads to `MATCH`; the construction never wires in a path that goes
nowhere. So a single motif's automaton is dead-state-free by
construction, and this pass is verified to be correctly a no-op
rather than silently broken (a broken implementation could just as
easily report 0 removals by never finding anything, which is why the
pass is still exercised and its count checked, not simply asserted to
be correct).

The pass remains in the pipeline for two reasons: it is still
necessary for general correctness (an automaton built by some other
route -- composing several motifs into one combined search, for
instance, which is a natural Phase 3 extension -- can absolutely
produce reachable dead states), and dropping it would mean the
minimization step in section 6 has to defend against dead states
itself instead of relying on an earlier pass to have ruled them out.

## 6. Minimization (`dfa_minimize`)

Iterative partition refinement, in the tradition of Moore's
algorithm: states start split only into accepting/non-accepting
blocks, and are repeatedly re-split whenever two states in the same
block transition to *different* blocks on some residue, until a full
pass produces no further split. This converges to the same unique
minimal DFA that Hopcroft's algorithm reaches; the difference is
purely asymptotic (this is the O(n<sup>2</sup>&middot;|&Sigma;|)
formulation rather than Hopcroft's O(n log n) worklist-driven one),
which is irrelevant at the state counts this project deals with
(single digits to a few hundred).

Concrete effect: `nglyc` (N-X-[ST], X&ne;P) compiles to exactly 4
states after minimization -- which is optimal, since the pattern has
exactly 4 distinguishable positions (before N, after N, after X, and
accept). `pyrimidine_run` ([^AG]{3,}) compiles to exactly 4 states:
three to count up to the minimum run length, and one accepting
self-loop absorbing any further pyrimidines. Both are hand-verifiable
minimal automata for their languages.

## 7. Literal-prefix extraction (`dfa_extract_literal_prefix`)

Walks forward from the start state as long as exactly one residue has
a live transition at the current state -- i.e. as long as the DFA's
own structure leaves no ambiguity about what has to come next. A
fully-literal motif like the EcoRI site `GAATTC` extracts its entire
six residues as the prefix (the DFA is a straight chain with exactly
one live transition at every step). A motif starting with an
ambiguity code, or an alternation, extracts an empty prefix
immediately, correctly.

This is used by the matcher as a cheap pre-filter (a `strncmp` against
the prefix before running the automaton at all), and is exactly the
kind of literal-prefix analysis Phase 3's native code generation will
build on for a proper skip-search over the hot scan loop.

## 8. Matching policy (`matcher.c`)

For every start position `p` in a sequence (skipping positions the
literal prefix already rules out), the DFA runs forward from `p`
until it dies or the sequence ends, bounded by the motif's Phase-1
inferred maximum length where that is finite. If the DFA passes
through one or more accepting states, the **longest** accepting
length is reported as one hit at `p`; the scan then continues at
`p+1` regardless, so overlapping hits are all reported. This matches
standard motif-search tools, which look for every occurrence rather
than only non-overlapping ones -- a repeat element or a densely
packed set of binding sites is exactly the case where overlapping
hits are the correct answer, not an artifact to suppress.

## 9. A real finding from testing, not just a demonstration

`tests/promoter.motif` searches for `polyA` (`"AATAAA"` allow 1
mismatch) over a sequence containing both an exact `AATAAA` and a run
of plain `A`s. The exact site is found, correctly. But `AAAAAA` is
also reported as a hit -- because it genuinely is Hamming-distance 1
from `AATAAA` (one substitution at the T position). This is not a
bug: it is the mismatch budget doing exactly what it was asked to do,
and it demonstrates a real, well-known problem in approximate
biological sequence matching -- permissive mismatch tolerance over
low-complexity or repetitive sequence produces a burst of technically
correct but biologically uninteresting hits. Real tools handle this
with a low-complexity pre-filter (e.g. dustmasker-style masking)
before searching; that is a natural, explicitly out-of-scope
extension noted here rather than silently avoided in the test data.

## 10. Edit distance (`nfa_apply_error_tolerance`)

`allow k edits` reuses the k+1 layers of Section 3. A CONSUME state
at layer i < k gains two more exits into layer i+1: a deletion
(epsilon, skipping the motif position) and an insertion (consume any
residue of the alphabet, staying at the same motif position). The
insertion mask is the alphabet's full residue set, never "any byte",
so alphabet typing carries through to approximate matching.

Insertion edges are left off the CONSUME states reachable from the
start before anything is consumed. An insertion there would only
reproduce every hit one residue earlier; search already has a free
start, so the edge adds duplicates and no information.

Correctness is checked against a brute-force dynamic program
(`tools/oracle_check.py`): for random IUPAC motifs, random DNA and
k = 1 or 2 under both models, the set of hit start positions must be
identical. 1,000 random cases (35,399 hit positions) agreed exactly.

Edits still produce clusters: the same site minus its first residue
(one deletion) is a separate, overlapping hit. `--merge-contained`
drops a hit wholly contained in the previous one. It is off by
default so the raw output stays comparable with the oracle.

## 11. Compile-time specificity (`specificity.c`)

Given a background model, the pass computes the probability that the
matcher reports a hit at one position of a random sequence, straight
from the minimised DFA: probability mass is pushed through the
transition table one residue at a time, mass reaching an accepting
state is counted and removed, and mass on a dead transition is
dropped. Multiplying by the number of positions gives the expected
number of chance hits. Backgrounds can be uniform, or order-0 or
order-1 Markov models estimated from a FASTA file.

`--specificity` prints the estimate for every motif and warns when a
motif expects more than `--max-rate` chance hits per Mb (default
1,000). Every search also prints the expected count for the searched
file next to the observed count.

Checks: exact agreement with closed forms (1/8192 for `TATAWAWR`,
19/4096 for a 6-mer with one mismatch, 1/8 for `[^AG]{3,}`), and
agreement within a few per cent with observed counts on 5 Mb of
random DNA (`tools/rate_check.py`). On CpG-depleted sequence the
order-0 background overestimates by 1.2x to 4x while order 1 is
within 1.5%, which is why order 1 is the default for searches.
