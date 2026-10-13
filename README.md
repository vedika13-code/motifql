# MotifQL

A domain-specific language and compiler, written in C, for searching
biological sequences for motifs. Motifs are written with IUPAC
ambiguity codes, composed from named parts, checked at compile time
against a declared alphabet (DNA, RNA or protein), compiled through
an NFA and a minimised DFA, and run over FASTA files.

```
alphabet DNA;
motif tata   = "TATAWAWR" allow 1 edits;
motif kozak  = ("A" | "G") "CCATGG";
motif promoter { @tata "N"{20,40} @kozak }
search promoter in "chr1.fa";
```

## Build and test

```
make            # builds build/motifql (C11, no dependencies)
make test       # 26 tests, including a brute-force oracle (needs python3)
make demo
```

## Useful flags

| Flag | Shows |
| --- | --- |
| `--tokens`, `--ast`, `--symbols` | front-end output: tokens, annotated syntax tree, symbol table |
| `--nfa`, `--dfa`, `--bytecode` | the automaton at each stage |
| `--specificity` | expected chance hits per motif, computed from the automaton |
| `--background FILE --markov 0\|1` | background model estimated from a FASTA file |
| `--max-rate R` | warn above R expected chance hits per Mb (default 1000) |
| `--merge-contained` | drop hits lying wholly inside the previous hit |
| `--no-search` | compile only |

## Validation tools

- `tools/oracle_check.py [trials] [seed]` compares hits with a brute-force dynamic program.
- `tools/rate_check.py [Mb] [seed]` compares predicted and observed chance hits on random DNA.

See `docs/LANGUAGE.md` for the language and `docs/PHASE2.md` for the compiler back end.
