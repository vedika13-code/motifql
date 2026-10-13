#include "dfa.h"

void dfa_init(Dfa *dfa, AlphabetKind alphabet) {
    memset(dfa, 0, sizeof(*dfa));
    dfa->alphabet = alphabet;
    dfa->alphabet_size = alphabet_size(alphabet);
}

void dfa_free(Dfa *dfa) {
    free(dfa->transitions);
    free(dfa->accepting);
    memset(dfa, 0, sizeof(*dfa));
}

/* ---- small integer-set helper, used to name a DFA state by the set
 * of NFA states it corresponds to (the textbook "subset" in subset
 * construction). Automata in this project are small enough that
 * linear search over discovered sets is entirely adequate; this is
 * not the asymptotic bottleneck anywhere in the pipeline. */

typedef struct {
    int *items;
    int count;
    int capacity;
} IntSet;

static void intset_free(IntSet *s) { free(s->items); s->items = NULL; s->count = s->capacity = 0; }

static bool intset_contains(const IntSet *s, int v) {
    for (int i = 0; i < s->count; i++) if (s->items[i] == v) return true;
    return false;
}

static void intset_add(IntSet *s, int v) {
    if (intset_contains(s, v)) return;
    if (s->count + 1 > s->capacity) {
        int new_cap = s->capacity == 0 ? 8 : s->capacity * 2;
        int *grown = realloc(s->items, sizeof(int) * (size_t)new_cap);
        if (!grown) { fprintf(stderr, "out of memory in subset construction\n"); exit(1); }
        s->items = grown;
        s->capacity = new_cap;
    }
    s->items[s->count++] = v;
}

static int int_cmp(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }

static void intset_sort(IntSet *s) { qsort(s->items, (size_t)s->count, sizeof(int), int_cmp); }

static bool intset_equal(const IntSet *a, const IntSet *b) {
    if (a->count != b->count) return false;
    for (int i = 0; i < a->count; i++) if (a->items[i] != b->items[i]) return false;
    return true;
}

/* ---- epsilon closure and move -------------------------------------------- */

static void epsilon_closure(const Nfa *nfa, IntSet *set) {
    /* Grow `set` in place by following NFA_SPLIT edges to a fixed
     * point. A simple worklist over the set's own growing item list
     * works because intset_add is a no-op on states already present,
     * so this terminates. */
    for (int i = 0; i < set->count; i++) {
        int s = set->items[i];
        if (nfa->states[s].kind != NFA_SPLIT) continue;
        int o1 = nfa->states[s].out1, o2 = nfa->states[s].out2;
        if (o1 >= 0) intset_add(set, o1);
        if (o2 >= 0) intset_add(set, o2);
    }
}

static IntSet move_on_residue(const Nfa *nfa, const IntSet *set, int residue) {
    IntSet result = {0};
    SymMask bit = (SymMask)1u << residue;
    for (int i = 0; i < set->count; i++) {
        const NfaState *s = &nfa->states[set->items[i]];
        if (s->kind == NFA_CONSUME && (s->mask & bit) && s->out1 >= 0) {
            intset_add(&result, s->out1);
        }
    }
    return result;
}

static bool set_has_match(const Nfa *nfa, const IntSet *set) {
    for (int i = 0; i < set->count; i++) {
        if (nfa->states[set->items[i]].kind == NFA_MATCH) return true;
    }
    return false;
}

/* ---- subset construction --------------------------------------------------- */

typedef struct {
    IntSet *keys;      /* one per discovered DFA state, parallel to below */
    int count;
    int capacity;
} DiscoveredStates;

static int find_or_add(DiscoveredStates *ds, IntSet key, bool *is_new) {
    for (int i = 0; i < ds->count; i++) {
        if (intset_equal(&ds->keys[i], &key)) {
            intset_free(&key);
            *is_new = false;
            return i;
        }
    }
    if (ds->count + 1 > ds->capacity) {
        int new_cap = ds->capacity == 0 ? 16 : ds->capacity * 2;
        IntSet *grown = realloc(ds->keys, sizeof(IntSet) * (size_t)new_cap);
        if (!grown) { fprintf(stderr, "out of memory in subset construction\n"); exit(1); }
        ds->keys = grown;
        ds->capacity = new_cap;
    }
    ds->keys[ds->count] = key;
    *is_new = true;
    return ds->count++;
}

void dfa_from_nfa(const Nfa *nfa, Dfa *out) {
    dfa_init(out, nfa->alphabet);
    int alpha = out->alphabet_size;

    DiscoveredStates ds = {0};

    IntSet start_set = {0};
    intset_add(&start_set, nfa->start);
    epsilon_closure(nfa, &start_set);
    intset_sort(&start_set);

    bool is_new;
    int start_index = find_or_add(&ds, start_set, &is_new);
    out->start = start_index;

    int processed = 0;
    int *transitions = NULL;
    int trans_capacity = 0;

    while (processed < ds.count) {
        int current = processed++;

        while ((current + 1) * alpha > trans_capacity) {
            int new_cap = trans_capacity == 0 ? alpha * 16 : trans_capacity * 2;
            int *grown = realloc(transitions, sizeof(int) * (size_t)new_cap);
            if (!grown) { fprintf(stderr, "out of memory in subset construction\n"); exit(1); }
            transitions = grown;
            trans_capacity = new_cap;
        }

        for (int r = 0; r < alpha; r++) {
            IntSet moved = move_on_residue(nfa, &ds.keys[current], r);
            if (moved.count == 0) {
                transitions[current * alpha + r] = DFA_DEAD;
                intset_free(&moved);
                continue;
            }
            epsilon_closure(nfa, &moved);
            intset_sort(&moved);
            bool added;
            int target = find_or_add(&ds, moved, &added);

            while (ds.count * alpha > trans_capacity) {
                int new_cap = trans_capacity == 0 ? alpha * 16 : trans_capacity * 2;
                while (new_cap < ds.count * alpha) new_cap *= 2;
                int *grown = realloc(transitions, sizeof(int) * (size_t)new_cap);
                if (!grown) { fprintf(stderr, "out of memory in subset construction\n"); exit(1); }
                transitions = grown;
                trans_capacity = new_cap;
            }
            transitions[current * alpha + r] = target;
            (void)added;
        }
    }

    out->count = ds.count;
    size_t nstates = ds.count > 0 ? (size_t)ds.count : 1;   /* count is never negative; this keeps -O2's range analysis quiet */
    out->transitions = realloc(transitions, sizeof(int) * nstates * (size_t)alpha);
    out->accepting = calloc(nstates, sizeof(bool));
    for (int i = 0; i < ds.count; i++) {
        out->accepting[i] = set_has_match(nfa, &ds.keys[i]);
        intset_free(&ds.keys[i]);
    }
    free(ds.keys);
}

/* ---- dead-state elimination ------------------------------------------------
 *
 * A state is "dead" if no accepting state is reachable from it at
 * all -- distinct from DFA_DEAD (which marks the ABSENCE of a
 * transition). Computed by backward reachability from the accepting
 * states over the transition relation, then anything not reached is
 * dead and gets pruned. */

int dfa_eliminate_dead_states(Dfa *dfa) {
    int n = dfa->count;
    bool *can_reach_accept = calloc((size_t)n, sizeof(bool));
    bool changed = true;

    while (changed) {
        changed = false;
        for (int s = 0; s < n; s++) {
            if (can_reach_accept[s]) continue;
            if (dfa->accepting[s]) { can_reach_accept[s] = true; changed = true; continue; }
            for (int r = 0; r < dfa->alphabet_size; r++) {
                int t = dfa_transition(dfa, s, r);
                if (t != DFA_DEAD && can_reach_accept[t]) { can_reach_accept[s] = true; changed = true; break; }
            }
        }
    }

    int removed = 0;
    int *remap = malloc(sizeof(int) * (size_t)n);
    int next = 0;
    for (int s = 0; s < n; s++) {
        if (can_reach_accept[s]) remap[s] = next++;
        else { remap[s] = DFA_DEAD; removed++; }
    }

    if (removed > 0) {
        int *new_trans = malloc(sizeof(int) * (size_t)next * (size_t)dfa->alphabet_size);
        bool *new_accept = malloc(sizeof(bool) * (size_t)next);
        for (int s = 0; s < n; s++) {
            if (remap[s] == DFA_DEAD) continue;
            new_accept[remap[s]] = dfa->accepting[s];
            for (int r = 0; r < dfa->alphabet_size; r++) {
                int t = dfa_transition(dfa, s, r);
                new_trans[remap[s] * dfa->alphabet_size + r] = (t == DFA_DEAD || !can_reach_accept[t]) ? DFA_DEAD : remap[t];
            }
        }
        free(dfa->transitions);
        free(dfa->accepting);
        dfa->transitions = new_trans;
        dfa->accepting = new_accept;
        dfa->count = next;
        dfa->start = (remap[dfa->start] == DFA_DEAD) ? DFA_DEAD : remap[dfa->start];
    }

    free(can_reach_accept);
    free(remap);
    return removed;
}

/* ---- minimization (Moore's partition-refinement algorithm) --------------- */

int dfa_minimize(Dfa *dfa) {
    int n = dfa->count;
    if (n == 0) return 0;

    int *block = malloc(sizeof(int) * (size_t)n);
    for (int s = 0; s < n; s++) block[s] = dfa->accepting[s] ? 1 : 0;
    int num_blocks = 2;
    {
        bool any0 = false, any1 = false;
        for (int s = 0; s < n; s++) { if (block[s] == 0) any0 = true; else any1 = true; }
        if (!any0 || !any1) num_blocks = 1;
    }

    bool changed = true;
    while (changed) {
        changed = false;

        typedef struct { int block; int *sig; } Row;
        Row *rows = malloc(sizeof(Row) * (size_t)n);
        for (int s = 0; s < n; s++) {
            rows[s].block = block[s];
            rows[s].sig = malloc(sizeof(int) * (size_t)dfa->alphabet_size);
            for (int r = 0; r < dfa->alphabet_size; r++) {
                int t = dfa_transition(dfa, s, r);
                rows[s].sig[r] = (t == DFA_DEAD) ? -1 : block[t];
            }
        }

        int *new_block = malloc(sizeof(int) * (size_t)n);
        for (int s = 0; s < n; s++) new_block[s] = -1;
        int next_block_id = 0;

        for (int old_b = 0; old_b < num_blocks; old_b++) {
            for (int s = 0; s < n; s++) {
                if (rows[s].block != old_b || new_block[s] != -1) continue;
                new_block[s] = next_block_id;
                for (int t = s + 1; t < n; t++) {
                    if (rows[t].block != old_b || new_block[t] != -1) continue;
                    bool same = true;
                    for (int r = 0; r < dfa->alphabet_size; r++) {
                        if (rows[s].sig[r] != rows[t].sig[r]) { same = false; break; }
                    }
                    if (same) new_block[t] = next_block_id;
                }
                next_block_id++;
            }
        }

        if (next_block_id != num_blocks) changed = true;
        else {
            for (int s = 0; s < n; s++) if (new_block[s] != block[s]) { changed = true; break; }
        }

        free(block);
        block = new_block;
        num_blocks = next_block_id;

        for (int s = 0; s < n; s++) free(rows[s].sig);
        free(rows);
    }

    int *new_trans = malloc(sizeof(int) * (size_t)num_blocks * (size_t)dfa->alphabet_size);
    bool *new_accept = calloc((size_t)num_blocks, sizeof(bool));
    bool *filled = calloc((size_t)num_blocks, sizeof(bool));

    for (int s = 0; s < n; s++) {
        int b = block[s];
        if (filled[b]) continue;
        filled[b] = true;
        new_accept[b] = dfa->accepting[s];
        for (int r = 0; r < dfa->alphabet_size; r++) {
            int t = dfa_transition(dfa, s, r);
            new_trans[b * dfa->alphabet_size + r] = (t == DFA_DEAD) ? DFA_DEAD : block[t];
        }
    }

    int removed = n - num_blocks;

    free(dfa->transitions);
    free(dfa->accepting);
    dfa->transitions = new_trans;
    dfa->accepting = new_accept;
    dfa->start = block[dfa->start];
    dfa->count = num_blocks;

    free(filled);
    free(block);
    return removed;
}

/* ---- literal-prefix extraction --------------------------------------------- */

void dfa_extract_literal_prefix(const Dfa *dfa, char *buf, size_t buf_size) {
    const char *residues = alphabet_residues(dfa->alphabet);
    size_t pos = 0;
    int current = dfa->start;
    bool *seen = calloc((size_t)(dfa->count > 0 ? dfa->count : 1), sizeof(bool));

    while (pos + 1 < buf_size && current != DFA_DEAD && !seen[current]) {
        seen[current] = true;
        int only_residue = -1;
        int live_count = 0;
        for (int r = 0; r < dfa->alphabet_size; r++) {
            if (dfa_transition(dfa, current, r) != DFA_DEAD) {
                live_count++;
                only_residue = r;
            }
        }
        if (live_count != 1) break;

        buf[pos++] = residues[only_residue];
        current = dfa_transition(dfa, current, only_residue);
    }
    buf[pos] = '\0';
    free(seen);
}

void dfa_print(const Dfa *dfa) {
    const char *residues = alphabet_residues(dfa->alphabet);
    printf("== DFA (%d states, start=%d) ==\n", dfa->count, dfa->start);
    printf("     %-8s", "state");
    for (int r = 0; r < dfa->alphabet_size; r++) printf(" %3c", residues[r]);
    printf("   accept\n");
    for (int s = 0; s < dfa->count; s++) {
        printf("     %-8d", s);
        for (int r = 0; r < dfa->alphabet_size; r++) {
            int t = dfa_transition(dfa, s, r);
            if (t == DFA_DEAD) printf("   .");
            else printf(" %3d", t);
        }
        printf("   %s\n", dfa->accepting[s] ? "yes" : "");
    }
}
