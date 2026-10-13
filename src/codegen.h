#ifndef MOTIFQL_CODEGEN_H
#define MOTIFQL_CODEGEN_H

#include "common.h"
#include "dfa.h"
#include "symtab.h"

/*
 * Target code generation. A MatcherProgram is the compiled artifact
 * a motif reduces to: the minimized DFA's transition table (the
 * "instructions" the matcher executes, one table lookup per residue
 * consumed) plus the metadata the matcher and any future optimizer
 * need -- inferred length bounds from Phase 1, and the literal prefix
 * extracted from the DFA's own structure.
 */

typedef struct {
    int *transitions;    /* num_states * alphabet_size, flattened, DFA_DEAD sentinel */
    bool *accepting;
    int num_states;
    int alphabet_size;
    AlphabetKind alphabet;
    int start_state;

    int min_len;
    int max_len;          /* LEN_INF if unbounded */

    char literal_prefix[64];

    char motif_name[64];
} MatcherProgram;

/* Runs the whole Phase 2 pipeline for one motif: Thompson's
 * construction, optional mismatch-tolerance layering, subset
 * construction, dead-state elimination, minimization, and literal-
 * prefix extraction, then packages the result. `stats`, if non-NULL,
 * is filled with the state counts at each stage for reporting. */
typedef struct {
    int nfa_states;
    int dfa_states_raw;
    int dfa_states_after_dead_elim;
    int dead_states_removed;
    int dfa_states_final;
    int states_removed_by_minimization;
} CompileStats;

void compile_motif_to_bytecode(Program *program, int motif_index, MatcherProgram *out, CompileStats *stats);

void matcher_program_free(MatcherProgram *mp);
void matcher_program_print(const MatcherProgram *mp);

#endif
