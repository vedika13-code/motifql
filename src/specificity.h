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
 * Background models are Markov chains of order k:
 *   order 0 -- residues independent, with fixed frequencies
 *   order 1 -- each residue depends on the one before it (captures
 *              effects such as CpG depletion in vertebrate genomes)
 *   order k -- each residue depends on the k before it (up to
 *              BG_MAX_ORDER for nucleotides, 1 for protein)
 * For order k the probability mass is tracked per (DFA state, last k
 * residues), so the cost grows as 4^k; k <= 6 keeps it small.
 */

#define BG_MAX 32
#define BG_MAX_ORDER 6

typedef struct {
    AlphabetKind alphabet;
    int size;
    int order;                    /* k */
    int contexts;                 /* size^k */
    double p0[BG_MAX];            /* residue frequencies */
    double *init;                 /* [contexts]: frequency of each k-residue context */
    double *trans;                /* [contexts * size]: P(next residue | context) */
    char source[160];             /* human-readable description for reports */
} Background;

void background_free(Background *bg);

void background_uniform(Background *bg, AlphabetKind alphabet);

/* Estimates an order-k background from the residues of a FASTA file.
 * Residues outside the alphabet (N runs) are skipped and break the
 * context. k is clamped to BG_MAX_ORDER (1 for protein). Returns false
 * if the file contains no usable residues. */
bool background_from_fasta(Background *bg, AlphabetKind alphabet, const FastaFile *f, int order,
                           const char *label);

/* Probability that the matcher reports a hit at one given start
 * position of a random sequence drawn from `bg`. Multiply by the
 * number of start positions for the expected number of chance hits. */
double specificity_hit_probability(const MatcherProgram *mp, const Background *bg);

/* Empirical alternative: the fraction of positions in a reference
 * sequence where the matcher actually reports a hit. Makes no Markov
 * assumption, so it captures repeats, CpG islands and other structure
 * a background model misses -- at the cost of scanning the reference
 * once per evaluation. The reference should NOT be the sequence being
 * searched, or the estimate simply restates the answer. */
double specificity_empirical_rate(const MatcherProgram *mp, const FastaFile *reference);

/* Where a chance-hit rate comes from: a background model, or (when
 * `reference` is set) an empirical scan of a reference sequence. */
typedef struct {
    const Background *background;
    const FastaFile *reference;
    const char *reference_label;
} RateSource;

double rate_of(const MatcherProgram *mp, const RateSource *src);
const char *rate_source_text(const RateSource *src, char *buf, size_t n);

/* Counts the start positions in a FASTA file that the matcher could
 * use (a residue of the alphabet at that position), so an expected
 * count can be compared like-for-like with the observed one. */
long background_valid_positions(const FastaFile *f, AlphabetKind alphabet);

#endif
