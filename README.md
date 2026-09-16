# MotifQL — Phase 1

A domain-specific language and compiler front end for biological
sequence motif search — written in C, for a Compiler Design course
project. **This package is Phase 1 only**: the front end (lexer,
parser, symbol table, semantic analyser). Phase 2 (automaton
construction and a FASTA matcher) and Phase 3 (native JIT) build on
this same front end but are separate deliverables.

MotifQL lets you describe a DNA/RNA/protein motif declaratively, with
IUPAC ambiguity codes as native, alphabet-checked symbols.

```
alphabet DNA;

motif tata     = "TATAWAWR";
motif kozak    { ("A" | "G") "CC" "ATG" "G" }
motif polyA    { "AATAAA" } allow 1 mismatch;
motif promoter { @tata "N"{20,40} @kozak }

search promoter in "chr1.fa";
```

## Why

Generic regular expressions don't know biology: a DNA pattern
containing `U`, or a protein pattern with a non-amino-acid letter,
simply compiles and silently never matches. MotifQL makes the
alphabet part of the language, so every symbol is checked against the
declared alphabet (DNA, RNA, or protein) at compile time, and
ambiguity codes like `W` (A or T) or `N` (any base) are native
symbols rather than something you expand by hand.

## What's in this package

Lexical analysis, recursive-descent syntax analysis, a symbol table,
and a five-pass semantic analyser (symbol resolution, cycle
detection, alphabet checking, quantifier-bound checking, and
match-length inference). See `docs/LANGUAGE.md` for the grammar and
design decisions, and `docs/PHASES.md` for how each classical
compiler phase maps onto this project.

## Build and run

Requires a C11 compiler and GNU Make. No external dependencies.

```
make            # builds build/motifql
make test       # runs all 8 front-end tests
make demo       # runs tokens/AST/symbol-table dumps against a sample motif
```

Try it directly:

```
./build/motifql --symbols tests/promoter.motif
./build/motifql --ast tests/promoter.motif
./build/motifql --tokens tests/promoter.motif
```

## Tests

8 tests in `tests/run_tests.sh`, covering the accepting path across
all three alphabets (DNA, RNA, protein) and five distinct rejection
paths (invalid symbol, bad quantifier bounds, dangling reference,
reference cycle, syntax/lexical errors with recovery).
