#ifndef MOTIFQL_DFA_H
#define MOTIFQL_DFA_H

#include "common.h"
#include "alphabet.h"
#include "nfa.h"

#define DFA_DEAD (-1)

typedef struct {
    int *transitions;   /* count * alphabet_size, flattened; DFA_DEAD = no edge */
    bool *accepting;
    int count;
    int alphabet_size;
    AlphabetKind alphabet;
    int start;
} Dfa;

void dfa_init(Dfa *dfa, AlphabetKind alphabet);
void dfa_free(Dfa *dfa);

static inline int dfa_transition(const Dfa *dfa, int state, int residue) {
    return dfa->transitions[state * dfa->alphabet_size + residue];
}

/* Subset construction: builds a DFA equivalent to `nfa` via the
 * standard epsilon-closure / subset algorithm. A DFA state absent
 * from the result is simply never constructed -- there is no
 * physical trap state, an unlisted transition already means DFA_DEAD.
 * (See dfa_eliminate_dead_states for the separate case of a state
 * that IS reachable and constructed, but can never reach acceptance.) */
void dfa_from_nfa(const Nfa *nfa, Dfa *out);

/* Finds every constructed state from which no accepting state is
 * reachable, and redirects all transitions into them to DFA_DEAD,
 * then physically removes them and renumbers the rest. Returns the
 * number of states removed. */
int dfa_eliminate_dead_states(Dfa *dfa);

/* Minimises via iterative partition refinement (Moore's algorithm):
 * states start split only by accept/non-accept, and are repeatedly
 * split further whenever two states in the same block transition to
 * different blocks on some residue, until no split fires. This reaches
 * the same unique minimal DFA as Hopcroft's algorithm; it is simply
 * the O(n^2 * |alphabet|) formulation rather than Hopcroft's O(n log n)
 * worklist-driven one, which does not matter at the automaton sizes
 * this project deals with. Returns the number of states removed. */
int dfa_minimize(Dfa *dfa);

/* Walks forward from the start state as long as exactly one residue
 * has a live (non-dead) transition out of the current state, i.e. as
 * long as the DFA's structure leaves no ambiguity about what the next
 * residue must be. Writes up to buf_size-1 residue characters into
 * buf and NUL-terminates. An empty result is normal and correct for
 * a motif whose very first position is already ambiguous. */
void dfa_extract_literal_prefix(const Dfa *dfa, char *buf, size_t buf_size);

void dfa_print(const Dfa *dfa);

#endif
