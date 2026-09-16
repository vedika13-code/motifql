# Compiler Phase Coverage

How each classical compiler phase is realised in MotifQL, and which
project phase implements it.

| Compiler phase | Realisation in MotifQL | Project phase | Source |
|---|---|---|---|
| Lexical analysis | Hand-written scanner: keywords, quoted sequence literals, quantifiers, class brackets, `@` references, two comment styles | 1 | `lexer.c` |
| Syntax analysis | Recursive-descent parser with precedence-free pattern grammar (alternation > concatenation > repetition > atom) | 1 | `parser.c` |
| Symbol table | Flat motif namespace with forward references, duplicate detection, use/def tracking | 1 | `symtab.c` |
| Semantic analysis | Five passes: symbol resolution, cycle detection, alphabet checking, quantifier bound checking, match-length inference | 1 | `sema.c` |
| Error detection & recovery | Panic-mode recovery, shared diagnostic sink across lexical/syntax/semantic phases, caret-positioned messages | 1 | `diagnostics.c`, `parser.c` |
| Intermediate representation | AST to NFA via Thompson's construction | 2 (not in this package) | `nfa.c` |
| Code optimisation | Subset construction, partition-refinement minimisation, dead-state elimination, literal-prefix extraction | 2 (not in this package) | `dfa.c` |
| Target code generation | DFA to matcher bytecode | 2 (not in this package) | `codegen.c` |
| Execution | Bytecode matcher scanning FASTA records | 2 (not in this package) | `matcher.c` |
| Runtime optimisation | Native x86-64 emission for the hot scan loop | 3 (planned) | `jit.c` |

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

This package contains Phase 1 only. Phase 2 (automaton construction
and the FASTA matcher) and Phase 3 (native JIT) build on this same
front end but are not included here.
