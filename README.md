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

Budgets can be chosen by the compiler:

```
motif tata  = "TATAWAWR" allow ? edits;
motif kozak = "RCCATGG"  allow ? edits;
motif promoter { @tata "N"{20,40} @kozak }
    require rate <= 2 per 100kb, divergence 10%;
```

## Build and test

```
make            # builds build/motifql (C11, no dependencies)
make test       # 41 tests, including a brute-force oracle (needs python3)
make demo
```

## Live demo

```
bash demo/run_demo.sh     # 10 steps, press Enter between them
```

## Useful flags

| Flag | Shows |
| --- | --- |
| `--tokens`, `--ast`, `--symbols` | front-end output: tokens, annotated syntax tree, symbol table |
| `--nfa`, `--dfa`, `--bytecode` | the automaton at each stage |
| `--specificity` | expected chance hits per motif, computed from the automaton |
| `--background FILE --markov K` | order-K background (0..6) estimated from a FASTA file |
| `--empirical FILE` | measure chance-hit rates by scanning a reference sequence |
| `--sensitivity` | chance hits and sensitivity (share of diverged sites found) per motif |
| `--divergence D --indel G` | mutation model for sensitivity (default 0.10, 0.02) |
| `--samples N --seed S` | Monte Carlo settings (default 5000, 1) |
| `--sites FILE` | known sites (one per record): sensitivity becomes the share found |
| `--count-only` | print hit counts only |
| `--max-rate R` | warn above R expected chance hits per Mb (default 1000) |
| `--merge-contained` | drop hits lying wholly inside the previous hit |
| `--no-search` | compile only |

## Validation tools

- `tools/oracle_check.py [trials] [seed]` compares hits with a brute-force dynamic program.
- `tools/rate_check.py [Mb] [seed]` compares predicted and observed chance hits on random DNA.
- `tools/eval_chr21.py data/chr21.fa` runs the real-genome evaluation (about 1 minute). chr21 is not
  bundled (47 MB); download hg38 chr21, for example from
  `https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/sequence/genome.fasta`,
  and save it as `data/chr21.fa`.
- `tools/eval_pas.py data/chr21.fa data/chr21_gencode.gtf` evaluates polyadenylation-signal
  detection at annotated 3' ends (GENCODE annotation subset from the same repository,
  `.../chr21/sequence/chr21_gencode.gtf`). Any GENCODE GTF works; the full annotation gives
  more sites.
- `tools/make_figures.py results figures` redraws the manuscript figures.

See `docs/LANGUAGE.md` for the language, `docs/PHASE2.md` for the back end and `docs/PHASE3.md` for sensitivity, budget synthesis and the evaluation.


## Author

Vedika Goyal (24BCB0049), B.Tech. CSE (Bioinformatics), VIT Vellore. Compiler Design Laboratory project, 2026.

## License

MIT. See `LICENSE`.
