#ifndef MOTIFQL_NFA_H
#define MOTIFQL_NFA_H

#include "common.h"
#include "alphabet.h"
#include "symtab.h"
#include "ast.h"

/*
 * NFA representation for Thompson's construction. Three kinds of
 * state, exactly as in the classical construction:
 *
 *   NFA_SPLIT   -- an epsilon branch to up to two successor states.
 *                  Used for alternation and for the loop points in
 *                  repetition. Consumes no input.
 *   NFA_CONSUME -- consumes exactly one input residue, provided it is
 *                  in `mask`, and moves to `out1`.
 *   NFA_MATCH   -- accepting state. No outgoing edges.
 *
 * This is deliberately the same three-state-kind formulation used in
 * Ken Thompson's original construction (and popularised by Russ Cox's
 * writeup of it): every regular operator -- concatenation,
 * alternation, and repetition -- reduces to wiring these three
 * primitives together, never to a fourth kind of state.
 */

typedef enum {
    NFA_SPLIT,
    NFA_CONSUME,
    NFA_MATCH,
} NfaKind;

typedef struct {
    NfaKind kind;
    SymMask mask;    /* NFA_CONSUME only */
    int out1;        /* NFA_SPLIT: first branch. NFA_CONSUME: target after consuming. */
    int out2;        /* NFA_SPLIT: second branch. Unused otherwise. */
    int pos;         /* NFA_CONSUME: 0-based residue position in the motif's own
                      * sequence, in construction order; -1 for synthetic states.
                      * Used by `allow k ... at P` (position-restricted errors). */
} NfaState;

typedef struct {
    NfaState *states;
    int count;
    int capacity;
    int start;
    AlphabetKind alphabet;

    /* Set when this NFA was built with `allow k mismatch`: the state
     * space is k+1 replicated layers (see nfa.c for the construction),
     * and `layer_count` / `states_per_layer` let later stages (and
     * diagnostics) make sense of which layer a state index belongs
     * to. layer_count is 1 for an NFA with no mismatch tolerance. */
    int layer_count;
    int states_per_layer;
    int next_pos;    /* next residue position to assign while building */
} Nfa;

void nfa_init(Nfa *nfa, AlphabetKind alphabet);
void nfa_free(Nfa *nfa);

/* Builds the plain (zero-mismatch) NFA for one motif's pattern via
 * Thompson's construction, inlining any @references (the reference
 * graph is already known acyclic by Phase 1's semantic analysis, so
 * this recursion always terminates). */
void nfa_build_for_motif(Program *program, int motif_index, Nfa *out);

/* Given a plain NFA and a mismatch budget k > 0, produces a new NFA
 * accepting the same language augmented with up to k substitutions.
 * See nfa.c for the layered construction. */
void nfa_apply_mismatch_tolerance(const Nfa *plain, int k, Nfa *out);

/* Generalisation used by the rest of the compiler: layers the NFA for
 * a budget of k errors under the given model. ERR_SUBSTITUTION gives
 * exactly nfa_apply_mismatch_tolerance; ERR_EDIT additionally adds
 * insertion and deletion edges (Levenshtein distance), with every
 * inserted or substituted residue restricted to the alphabet. */
void nfa_apply_error_tolerance(const Nfa *plain, int k, ErrorModel model, Nfa *out);

/* As above, but errors may only occur at the residue positions whose
 * bits are set in `positions` (bit i = position i+1 as written).
 * positions == 0 means every position. */
void nfa_apply_error_tolerance_at(const Nfa *plain, int k, ErrorModel model, uint32_t positions, Nfa *out);

void nfa_print(const Nfa *nfa);

#endif
