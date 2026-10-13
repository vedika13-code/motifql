# Compiler Phase Coverage

How each classical compiler phase is realised in MotifQL, and which
project phase implements it.

| Compiler phase | Realisation in MotifQL | Project phase | Source |
|---|---|---|---|
| Lexical analysis | Hand-written scanner: keywords, quoted sequence literals, quantifiers, class brackets, `@` references, two comment styles | 1 | `lexer.c` |
| Syntax analysis | Recursive-descent parser with precedence-free pattern grammar (alternation > concatenation > repetition > atom) | 1 | `parser.c` |
| Symbol table | Flat motif namespace with forward references, duplicate detection, use/def tracking, reference-to-index binding | 1 | `symtab.c` |
| Semantic analysis | Five passes: symbol resolution, cycle detection, alphabet checking, quantifier bound checking, match-length inference | 1 | `sema.c` |
| Error detection & recovery | Panic-mode recovery, shared diagnostic sink across lexical/syntax/semantic phases, caret-positioned messages | 1 | `diagnostics.c`, `parser.c` |
| Intermediate representation | AST to NFA via **Thompson's construction**; epsilon transitions for alternation and repetition | 2 | `nfa.c` |
| Code optimisation | **Subset construction** (NFA to DFA), **partition-refinement minimisation** (Moore's algorithm), dead-state elimination, literal-prefix extraction for scan skipping | 2 | `dfa.c` |
| Target code generation | DFA to matcher bytecode; transition tables laid out per state over the alphabet | 2 | `codegen.c` |
| Runtime / execution | Bytecode matcher scanning FASTA records, reporting hit coordinates | 2 | `matcher.c` |
| Runtime optimisation | Native x86-64 emission for the hot scan loop, benchmarked against the bytecode matcher | 3 | `jit.c` |

## Why the automata work belongs here

The regular-expression-to-NFA-to-DFA pipeline is the same construction
a lexer generator such as Lex or Flex performs internally to turn token
patterns into a scanner. Implementing it directly means the project
covers that theory as a first-class component rather than as something
delegated to a tool.

The biological framing is not decoration. IUPAC ambiguity codes make
alphabet-sensitive symbol resolution a real semantic-analysis problem
rather than a trivial character check, and representing symbols as
residue bitmasks means an ambiguity code costs exactly the same at
match time as a concrete base.

## Phase 1 status

Implemented and tested. The front end accepts a `.motif` file and
produces, for each stage, output that can be inspected independently:

- `--tokens` prints the token stream (lexical analysis)
- `--ast` prints the annotated syntax tree with resolved symbol masks
  and inferred length ranges (syntax + semantic analysis)
- `--symbols` prints the symbol table with resolved references,
  mismatch budgets, and use/def status

Eight test cases in `tests/` cover the accepting path across all three
alphabets and each distinct rejection path. Run with `make test`.

## Phase 2 status

Implemented and tested. Given a Phase-1-validated motif, the compiler
builds an NFA (Thompson's construction, with a layered extension for
`allow k mismatch`), determinises it (subset construction), removes
unreachable-to-accept states, minimises via partition refinement,
extracts a literal prefix from the resulting DFA's own structure, and
packages the whole thing as a `MatcherProgram` that a bytecode matcher
runs against FASTA records.

- `--nfa` prints the constructed NFA (states, kinds, edges)
- `--dfa` prints subset-construction, dead-state-elimination, and
  minimization statistics, then the minimized transition table
- `--bytecode` prints the final `MatcherProgram`: transition table,
  accepting states, inferred length bounds, and literal prefix
- with no dump flags, `search` statements run against real FASTA
  files and report hit coordinates

See `docs/PHASE2.md` for the algorithms, including an honest negative
result (dead-state elimination removes zero states on every test
motif, which is a property of Thompson's construction, not a broken
pass) and a genuine finding from testing (mismatch tolerance
over-matching on repetitive sequence). Six additional test cases
verify exact, hand-predicted hit positions -- not just pass/fail --
for a plain literal, an IUPAC class embedded in a literal, a
composite motif built from `@references` and a bounded spacer, and
mismatch-tolerant matching. Run with `make test`.

Phase 2 was then extended to follow the updated research proposal:

- `allow k edits`: alphabet-typed edit distance (substitutions,
  insertions, deletions), with new semantic checks for edit budgets
- a compile-time specificity pass (`--specificity`) that predicts
  chance hits per Mb from the minimised DFA and a background model
- validation tools: a brute-force oracle and a predicted-vs-observed
  rate check on random DNA

The suite now has 26 tests.

## Phase 3 plan

- sensitivity analysis: probability that a diverged true site is still
  found, under a stated substitution and indel model
- per-sub-motif error budgets and budget synthesis from a
  `require rate <= ...` clause
- evaluation on real sequence (human chromosome 21), final report
