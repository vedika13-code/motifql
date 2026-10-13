#ifndef MOTIFQL_SENSITIVITY_H
#define MOTIFQL_SENSITIVITY_H

#include "common.h"
#include "symtab.h"
#include "codegen.h"
#include "fasta.h"

/*
 * Compile-time sensitivity analysis.
 *
 * Question answered: if a real site of this motif has diverged from
 * the consensus -- some residues substituted, a few inserted or lost --
 * how likely is the compiled matcher to still find it?
 *
 * Method (Monte Carlo, fixed seed so results are reproducible):
 *   1. draw a random exact instance of the motif from its syntax tree
 *      (a residue from each IUPAC mask, a branch of each alternation,
 *      a count for each repeat, referenced motifs expanded in place);
 *   2. pass it through a mutation channel: each residue is deleted
 *      with probability indel/2, has a random residue inserted before
 *      it with probability indel/2, and is substituted by a different
 *      residue of the alphabet with probability `substitution`;
 *   3. run the compiled matcher on the result; the site counts as
 *      found if the matcher reports any hit in it.
 * Sensitivity is the fraction found, with a 95% Wilson interval.
 *
 * Because every budget assignment is scored on the SAME random sites
 * (same seed), comparisons between assignments are paired, and a
 * larger budget can never score lower than a smaller one: it accepts
 * a superset of sequences.
 */

typedef struct {
    double substitution;   /* per-residue substitution probability */
    double indel;          /* per-residue insertion + deletion probability */
} MutationModel;

typedef struct {
    double sensitivity;
    double ci_low, ci_high;
    int samples;
    int found;
} SensitivityResult;

/* Draws one random exact instance of a motif into buf (NUL-terminated).
 * Returns its length, or -1 if it would not fit. Exposed for tests. */
int sensitivity_sample_instance(const Program *program, int motif_index, char *buf, int cap, uint64_t *rng);

SensitivityResult sensitivity_estimate(const Program *program, int motif_index, const MatcherProgram *mp,
                                       const MutationModel *model, int samples, uint64_t seed);

/* Empirical alternative: the share of records in a FASTA file of known
 * sites (one site, with a little flanking sequence, per record) in
 * which the matcher reports at least one hit. */
SensitivityResult sensitivity_on_sites(const MatcherProgram *mp, const FastaFile *sites);

#endif
