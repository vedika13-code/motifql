#ifndef MOTIFQL_SEMA_H
#define MOTIFQL_SEMA_H

#include "common.h"
#include "symtab.h"
#include "diagnostics.h"

/*
 * Semantic analysis. Runs five distinct checks over the parsed
 * program, in an order forced by their dependencies:
 *
 *   1. Symbol resolution   -- bind every @name to a motif index, and
 *                             bind each search statement's target.
 *   2. Cycle detection     -- reject reference cycles (@a -> @b -> @a),
 *                             which describe infinitely large patterns.
 *                             Must precede (3) and (5), both of which
 *                             recurse through references and would not
 *                             terminate on a cyclic graph.
 *   3. Alphabet checking   -- resolve every residue character against
 *                             the declared alphabet, rejecting symbols
 *                             that are invalid in it (e.g. 'U' in DNA).
 *   4. Quantifier checking -- reject bounds where min > max, and
 *                             {0,0}, which can never match anything.
 *   5. Length inference    -- compute each pattern's min/max match
 *                             length bottom-up, then use it to reject
 *                             patterns that can match the empty string
 *                             and mismatch budgets that exceed the
 *                             motif's own minimum length.
 *
 * Returns true only when no errors were reported. Warnings (unused
 * motifs) do not fail the compile.
 */
/* Recomputes inferred match lengths, e.g. after budget synthesis has
 * changed some motifs' error budgets. */
void sema_recompute_lengths(Program *program);

/* Index of the one motif whose `require` clause governs the free
 * budget of motif `free_index`: -1 if none, -2 if more than one. */
int sema_governing_requirement(Program *program, int free_index);

/* Length of a motif that is a fixed run of residues, or -1. */
int sema_linear_length(const MotifDef *m);

bool analyze_program(Program *program, Diagnostics *diag);

#endif
