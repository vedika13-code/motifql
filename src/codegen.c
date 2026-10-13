#include "codegen.h"
#include "nfa.h"

void compile_motif_to_bytecode(Program *program, int motif_index, MatcherProgram *out, CompileStats *stats) {
    MotifDef *m = &program->motifs[motif_index];

    Nfa plain;
    nfa_build_for_motif(program, motif_index, &plain);

    Nfa *final_nfa = &plain;
    Nfa layered;
    bool built_layered = false;
    if (m->mismatch_budget > 0) {
        nfa_apply_error_tolerance_at(&plain, m->mismatch_budget, m->error_model, m->error_positions, &layered);
        final_nfa = &layered;
        built_layered = true;
    }

    if (stats) stats->nfa_states = final_nfa->count;

    Dfa dfa;
    dfa_from_nfa(final_nfa, &dfa);
    if (stats) stats->dfa_states_raw = dfa.count;

    int dead_removed = dfa_eliminate_dead_states(&dfa);
    if (stats) {
        stats->dead_states_removed = dead_removed;
        stats->dfa_states_after_dead_elim = dfa.count;
    }

    int before_min = dfa.count;
    dfa_minimize(&dfa);
    if (stats) {
        stats->states_removed_by_minimization = before_min - dfa.count;
        stats->dfa_states_final = dfa.count;
    }

    memset(out, 0, sizeof(*out));
    out->num_states = dfa.count;
    out->alphabet_size = dfa.alphabet_size;
    out->alphabet = dfa.alphabet;
    out->start_state = dfa.start;
    out->min_len = m->min_len;
    out->max_len = m->max_len;
    /* Insertions can make a match up to k residues longer than the
     * motif's own longest form; the matcher uses max_len as its scan
     * cut-off, so it has to know. */
    if (m->error_model == ERR_EDIT && m->max_len != LEN_INF) out->max_len = m->max_len + m->mismatch_budget;
    if (m->error_model == ERR_EDIT) out->min_len = m->min_len > m->mismatch_budget ? m->min_len - m->mismatch_budget : 1;

    size_t trans_bytes = sizeof(int) * (size_t)dfa.count * (size_t)dfa.alphabet_size;
    out->transitions = malloc(trans_bytes);
    memcpy(out->transitions, dfa.transitions, trans_bytes);
    out->accepting = malloc(sizeof(bool) * (size_t)dfa.count);
    memcpy(out->accepting, dfa.accepting, sizeof(bool) * (size_t)dfa.count);

    dfa_extract_literal_prefix(&dfa, out->literal_prefix, sizeof(out->literal_prefix));

    int name_len = m->name_length < 63 ? m->name_length : 63;
    memcpy(out->motif_name, m->name, (size_t)name_len);
    out->motif_name[name_len] = '\0';

    dfa_free(&dfa);
    if (built_layered) nfa_free(&layered);
    nfa_free(&plain);
}

void matcher_program_free(MatcherProgram *mp) {
    free(mp->transitions);
    free(mp->accepting);
    memset(mp, 0, sizeof(*mp));
}

void matcher_program_print(const MatcherProgram *mp) {
    const char *residues = alphabet_residues(mp->alphabet);
    printf("== matcher bytecode: %s ==\n", mp->motif_name);
    printf("states: %d   start: %d   alphabet: %s\n",
           mp->num_states, mp->start_state, alphabet_name(mp->alphabet));
    if (mp->max_len == LEN_INF) {
        printf("length bounds: %d..inf\n", mp->min_len);
    } else {
        printf("length bounds: %d..%d\n", mp->min_len, mp->max_len);
    }
    printf("literal prefix: %s\n", mp->literal_prefix[0] ? mp->literal_prefix : "(none)");
    printf("\n     %-8s", "state");
    for (int r = 0; r < mp->alphabet_size; r++) printf(" %3c", residues[r]);
    printf("   accept\n");
    for (int s = 0; s < mp->num_states; s++) {
        printf("     %-8d", s);
        for (int r = 0; r < mp->alphabet_size; r++) {
            int t = mp->transitions[s * mp->alphabet_size + r];
            if (t == DFA_DEAD) printf("   .");
            else printf(" %3d", t);
        }
        printf("   %s\n", mp->accepting[s] ? "yes" : "");
    }
    printf("\n");
}
