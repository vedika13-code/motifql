#ifndef MOTIFQL_SYNTHESIS_H
#define MOTIFQL_SYNTHESIS_H

#include "common.h"
#include "symtab.h"
#include "diagnostics.h"
#include "specificity.h"
#include "sensitivity.h"

/*
 * Budget synthesis: settles every `allow ?` budget from the `require`
 * clause that governs it.
 *
 * For a motif M with `require rate <= R per unit`, let F be the free
 * budgets in M and in the motifs it references. Each budget ranges
 * over 0..cap, where cap keeps at least one residue of the motif
 * unedited and never exceeds MAX_SYNTH_BUDGET. The compiler looks for
 * the assignment with the highest sensitivity whose expected chance
 * hit rate is at most R.
 *
 * Both measures only grow when any budget grows (a larger budget
 * accepts a superset of sequences), which gives two shortcuts:
 *   - an assignment at or above a known infeasible one, in every
 *     coordinate, is infeasible too and is never compiled;
 *   - the best feasible assignment is a MAXIMAL feasible one (raising
 *     any single budget breaks the target), so sensitivity is only
 *     estimated for those, plus the minimal infeasible assignments
 *     shown as "what one more error would cost".
 *
 * With `allow ? ... at ?` the positions are free too: each option is a
 * (budget, set of positions) pair, and option A is above option B when
 * A's budget is at least B's and A's positions include B's. Pruning and
 * maximality use that order, so they stay exact.
 *
 * Sensitivity comes from the Monte Carlo mutation model, or, with
 * --sites FILE, from the share of known sites the motif finds -- which
 * lets the compiler learn WHERE errors are tolerated from examples.
 *
 * A `require` clause with no free budgets is simply checked: a motif
 * that exceeds its required rate is a compile-time error.
 */

#define MAX_SYNTH_BUDGET 3
#define MAX_FREE_BUDGETS 6
#define MAX_CANDIDATES 20000   /* product of every free motif's options */

typedef struct {
    const Background *background;
    RateSource rates;         /* how chance-hit rates are measured */
    MutationModel model;      /* substitution may be overridden per clause */
    int samples;
    uint64_t seed;
    const FastaFile *sites;   /* known sites: if set, sensitivity is measured on them */
    const char *sites_label;
} SynthesisOptions;

/* Returns false if some contract cannot be met (errors reported). */
bool synthesize_budgets(Program *program, Diagnostics *diag, const SynthesisOptions *opts);

#endif
