#ifndef MOTIFQL_SPECIFICITY_H
#define MOTIFQL_SPECIFICITY_H

#include "common.h"
#include "alphabet.h"
#include "codegen.h"
#include "fasta.h"

/*
 * Compile-time specificity analysis.
 *
 * Question answered: if the searched sequence were random -- drawn
 * from a background model with the same base composition as the real
 * input -- how often would this motif report a hit purely by chance?
 *
 * The answer is computed from the compiled matcher itself, not by
 * scanning anything. The matcher reports a hit at start position p
 * exactly when its DFA, run from the start state on the residues from
 * p onward, reaches an accepting state. So the chance of a hit at one
 * position is the probability that a random walk through the DFA --
 * each step choosing residue r with the background probability of r --
 * reaches acceptance before it dies. That probability is found by
 * pushing probability mass through the transition table one residue
 * at a time: mass that reaches an accepting state is counted and
 * removed, mass that falls into a dead transition is dropped, and the
 * rest carries on. This is the pattern-Markov-chain idea of Nicodeme,
 * Salvy and Flajolet (2002) and Nuel (2008), applied to the automaton
 * MotifQL already builds.
 *
 * Two background models:
 *   order 0 -- residues independent, with fixed frequencies
 *   order 1 -- each residue depends on the one before it (captures
 *              effects such as CpG depletion in vertebrate genomes)
 */

#define BG_MAX 32

typedef struct {
    AlphabetKind alphabet;
    int size;
    int order;                    /* 0 or 1 */
    double p0[BG_MAX];            /* residue frequencies */
    double p1[BG_MAX][BG_MAX];    /* p1[a][b] = P(next = b | previous = a); order 1 only */
    char source[128];             /* human-readable description for reports */
} Background;

void background_uniform(Background *bg, AlphabetKind alphabet);

/* Estimates an order-0 or order-1 background from the residues of a
 * FASTA file. Residues outside the alphabet (N runs, soft gaps) are
 * skipped and also break order-1 context. Returns false if the file
 * contains no usable residues. */
bool background_from_fasta(Background *bg, AlphabetKind alphabet, const FastaFile *f, int order,
                           const char *label);

/* Probability that the matcher reports a hit at one given start
 * position of a random sequence drawn from `bg`. Multiply by the
 * number of start positions for the expected number of chance hits. */
double specificity_hit_probability(const MatcherProgram *mp, const Background *bg);

/* Counts the start positions in a FASTA file that the matcher could
 * use (a residue of the alphabet at that position), so an expected
 * count can be compared like-for-like with the observed one. */
long background_valid_positions(const FastaFile *f, AlphabetKind alphabet);

#endif
