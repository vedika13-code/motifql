#include "nfa.h"

void nfa_init(Nfa *nfa, AlphabetKind alphabet) {
    memset(nfa, 0, sizeof(*nfa));
    nfa->alphabet = alphabet;
    nfa->layer_count = 1;
}

void nfa_free(Nfa *nfa) {
    free(nfa->states);
    memset(nfa, 0, sizeof(*nfa));
}

static int nfa_add_state(Nfa *nfa, NfaKind kind) {
    if (nfa->count + 1 > nfa->capacity) {
        int new_cap = nfa->capacity == 0 ? 16 : nfa->capacity * 2;
        NfaState *grown = realloc(nfa->states, sizeof(NfaState) * (size_t)new_cap);
        if (!grown) { fprintf(stderr, "out of memory building NFA\n"); exit(1); }
        nfa->states = grown;
        nfa->capacity = new_cap;
    }
    NfaState *s = &nfa->states[nfa->count];
    s->kind = kind;
    s->mask = 0;
    s->out1 = -1;
    s->out2 = -1;
    s->pos = -1;
    return nfa->count++;
}

/* ---- fragment / patch-list machinery --------------------------------------
 *
 * A Frag is a partially built piece of NFA: `start` is its entry
 * state, and `dangling` lists the (state, slot) output pointers that
 * are not yet wired to anything -- exactly the classical Thompson
 * construction's "dangling arrow list", threaded through recursive
 * descent over the AST instead of over a token stream.
 */

typedef struct { int state; int slot; } Dangling; /* slot: 1 or 2, matching out1/out2 */

typedef struct {
    Dangling *items;
    int count;
    int capacity;
} PatchList;

static void patchlist_push(PatchList *pl, int state, int slot) {
    if (pl->count + 1 > pl->capacity) {
        int new_cap = pl->capacity == 0 ? 4 : pl->capacity * 2;
        Dangling *grown = realloc(pl->items, sizeof(Dangling) * (size_t)new_cap);
        if (!grown) { fprintf(stderr, "out of memory building patch list\n"); exit(1); }
        pl->items = grown;
        pl->capacity = new_cap;
    }
    pl->items[pl->count].state = state;
    pl->items[pl->count].slot = slot;
    pl->count++;
}

static void patchlist_append(PatchList *dst, const PatchList *src) {
    for (int i = 0; i < src->count; i++) patchlist_push(dst, src->items[i].state, src->items[i].slot);
}

static void patchlist_free(PatchList *pl) {
    free(pl->items);
    pl->items = NULL;
    pl->count = 0;
    pl->capacity = 0;
}

static void patch(Nfa *nfa, PatchList *pl, int target) {
    for (int i = 0; i < pl->count; i++) {
        NfaState *s = &nfa->states[pl->items[i].state];
        if (pl->items[i].slot == 1) s->out1 = target;
        else s->out2 = target;
    }
    patchlist_free(pl);
}

typedef struct {
    int start;
    PatchList dangling;
} Frag;

/* An epsilon pass-through: a SPLIT whose both arrows are dangling.
 * Entering it and leaving it consumes no input, which is exactly what
 * a zero-repetition (the `{0,...}` case, or an empty concat base)
 * needs to represent. */
static Frag frag_epsilon(Nfa *nfa) {
    int s = nfa_add_state(nfa, NFA_SPLIT);
    Frag f;
    f.start = s;
    f.dangling = (PatchList){0};
    patchlist_push(&f.dangling, s, 1);
    patchlist_push(&f.dangling, s, 2);
    return f;
}

static Frag frag_symbol(Nfa *nfa, SymMask mask) {
    int s = nfa_add_state(nfa, NFA_CONSUME);
    nfa->states[s].mask = mask;
    nfa->states[s].pos = nfa->next_pos++;
    Frag f;
    f.start = s;
    f.dangling = (PatchList){0};
    patchlist_push(&f.dangling, s, 1);
    return f;
}

static Frag frag_concat(Nfa *nfa, Frag a, Frag b) {
    patch(nfa, &a.dangling, b.start);
    Frag f;
    f.start = a.start;
    f.dangling = b.dangling;
    return f;
}

static Frag frag_alt(Nfa *nfa, Frag a, Frag b) {
    int s = nfa_add_state(nfa, NFA_SPLIT);
    nfa->states[s].out1 = a.start;
    nfa->states[s].out2 = b.start;
    Frag f;
    f.start = s;
    f.dangling = (PatchList){0};
    patchlist_append(&f.dangling, &a.dangling);
    patchlist_append(&f.dangling, &b.dangling);
    patchlist_free(&a.dangling);
    patchlist_free(&b.dangling);
    return f;
}

/* Zero-or-more of an already-built fragment (used only for the
 * unbounded tail of a repetition -- see build_repeat). */
static Frag frag_star(Nfa *nfa, Frag a) {
    int s = nfa_add_state(nfa, NFA_SPLIT);
    nfa->states[s].out1 = a.start;
    patch(nfa, &a.dangling, s);
    Frag f;
    f.start = s;
    f.dangling = (PatchList){0};
    patchlist_push(&f.dangling, s, 2);
    return f;
}

/* Zero-or-one of an already-built fragment. */
static Frag frag_quest(Nfa *nfa, Frag a) {
    int s = nfa_add_state(nfa, NFA_SPLIT);
    nfa->states[s].out1 = a.start;
    Frag f;
    f.start = s;
    f.dangling = (PatchList){0};
    patchlist_push(&f.dangling, s, 2);
    patchlist_append(&f.dangling, &a.dangling);
    patchlist_free(&a.dangling);
    return f;
}

static Frag build_fragment(Nfa *nfa, Program *program, Node *node);

/* Copies a complete, separately built NFA into `nfa` as a fragment.
 * Every MATCH state of the copy becomes an epsilon exit (a SPLIT whose
 * first arrow is left dangling), so whatever follows the fragment is
 * reachable from every layer the copy finished in. Used for a
 * referenced motif that carries its own error budget: the reference
 * is built and layered on its own, then dropped into place, so its
 * budget applies to its own span and nowhere else. */
static Frag frag_embed(Nfa *nfa, const Nfa *sub) {
    int offset = nfa->count;
    Frag f;
    f.dangling = (PatchList){0};
    for (int i = 0; i < sub->count; i++) {
        const NfaState *src = &sub->states[i];
        int s = nfa_add_state(nfa, src->kind == NFA_MATCH ? NFA_SPLIT : src->kind);
        if (src->kind == NFA_MATCH) {
            patchlist_push(&f.dangling, s, 1);
        } else {
            nfa->states[s].mask = src->mask;
            nfa->states[s].out1 = src->out1 >= 0 ? src->out1 + offset : -1;
            nfa->states[s].out2 = src->out2 >= 0 ? src->out2 + offset : -1;
        }
    }
    f.start = sub->start + offset;
    return f;
}

/* Bounded/unbounded repetition, unrolled into `min` mandatory copies
 * followed either by one starred copy (unbounded) or (max-min)
 * quest-wrapped optional copies (bounded). Each copy is a FRESH call
 * into build_fragment, never a shared sub-NFA: Thompson's
 * construction gives every occurrence its own states, which is what
 * keeps the automaton acyclic-by-construction except for the
 * deliberate loop points frag_star introduces. */
static Frag build_repeat(Nfa *nfa, Program *program, Node *child, int mn, int mx) {
    Frag result = frag_epsilon(nfa);

    for (int i = 0; i < mn; i++) {
        result = frag_concat(nfa, result, build_fragment(nfa, program, child));
    }

    if (mx == REPEAT_INF) {
        Frag tail = frag_star(nfa, build_fragment(nfa, program, child));
        result = frag_concat(nfa, result, tail);
    } else {
        for (int i = 0; i < mx - mn; i++) {
            Frag opt = frag_quest(nfa, build_fragment(nfa, program, child));
            result = frag_concat(nfa, result, opt);
        }
    }
    return result;
}

static Frag build_fragment(Nfa *nfa, Program *program, Node *node) {
    switch (node->type) {
        case NODE_LITERAL: {
            if (node->mask_count == 0) return frag_epsilon(nfa);
            Frag result = frag_symbol(nfa, node->masks[0]);
            for (int i = 1; i < node->mask_count; i++) {
                result = frag_concat(nfa, result, frag_symbol(nfa, node->masks[i]));
            }
            return result;
        }
        case NODE_CLASS:
            return frag_symbol(nfa, node->class_mask);

        case NODE_CONCAT: {
            Frag result = build_fragment(nfa, program, node->children[0]);
            for (int i = 1; i < node->child_count; i++) {
                result = frag_concat(nfa, result, build_fragment(nfa, program, node->children[i]));
            }
            return result;
        }
        case NODE_ALTERNATE: {
            Frag result = build_fragment(nfa, program, node->children[0]);
            for (int i = 1; i < node->child_count; i++) {
                result = frag_alt(nfa, result, build_fragment(nfa, program, node->children[i]));
            }
            return result;
        }
        case NODE_REPEAT:
            return build_repeat(nfa, program, node->children[0], node->repeat_min, node->repeat_max);

        case NODE_REF: {
            /* Inlined expansion: safe because Phase 1 already proved
             * the reference graph acyclic. A referenced motif with its
             * own `allow k ...` keeps that budget inside the composite
             * (per-part budgets): it is built and layered separately,
             * then embedded. */
            MotifDef *target = &program->motifs[node->ref_index];
            if (target->mismatch_budget <= 0) return build_fragment(nfa, program, target->pattern);
            Nfa sub_plain, sub_layered;
            nfa_build_for_motif(program, node->ref_index, &sub_plain);
            nfa_apply_error_tolerance_at(&sub_plain, target->mismatch_budget, target->error_model,
                                         target->error_positions, &sub_layered);
            Frag f = frag_embed(nfa, &sub_layered);
            nfa_free(&sub_layered);
            nfa_free(&sub_plain);
            return f;
        }
    }
    return frag_epsilon(nfa);
}

void nfa_build_for_motif(Program *program, int motif_index, Nfa *out) {
    MotifDef *m = &program->motifs[motif_index];
    nfa_init(out, program->alphabet);

    Frag body = build_fragment(out, program, m->pattern);
    int accept = nfa_add_state(out, NFA_MATCH);
    patch(out, &body.dangling, accept);

    out->start = body.start;
    out->layer_count = 1;
    out->states_per_layer = out->count;
}

/* ---- mismatch-tolerance layering -------------------------------------------
 *
 * Standard construction for Hamming-distance (substitution-only)
 * approximate matching: replicate the whole automaton into k+1
 * layers, 0..k, where layer i means "i substitutions used so far".
 * A SPLIT state is copied per layer with same-layer targets (it
 * consumes nothing, so it costs nothing). A CONSUME state (mask, t)
 * becomes, at layer i:
 *
 *   - a normal edge on any residue IN mask, to (i, t)          [no cost]
 *   - if i < k, a substitution edge on any residue NOT in mask,
 *     to (i+1, t)                                              [+1 cost]
 *
 * Both edges are represented by branching through a SPLIT into two
 * synthetic CONSUME states, since NfaState only carries one mask.
 * MATCH is accepting at every layer: reaching the end having used
 * at most k substitutions is success regardless of which layer.
 */
void nfa_apply_mismatch_tolerance(const Nfa *plain, int k, Nfa *out) {
    nfa_apply_error_tolerance(plain, k, ERR_SUBSTITUTION, out);
}

/* ---- edit-distance layering ------------------------------------------------
 *
 * Same k+1 layers as above, but a CONSUME state (mask, t) at layer
 * i < k now gets FOUR ways forward instead of two:
 *
 *   match         residue in mask            -> (i,   t)   [free]
 *   substitution  residue in (alphabet-mask) -> (i+1, t)   [+1]
 *   deletion      epsilon, nothing consumed  -> (i+1, t)   [+1]
 *                 (the pattern position is skipped)
 *   insertion     any residue of the alphabet -> (i+1, j)  [+1]
 *                 (an extra residue in the sequence; the
 *                  pattern position is NOT advanced)
 *
 * "Any residue of the alphabet" is the alphabet-typing point: the
 * insertion edge's mask is the full residue set of the declared
 * alphabet, so a DNA motif can absorb an extra A/C/G/T but never a U
 * or an unknown symbol.
 *
 * One deliberate restriction: no insertion edges on the CONSUME
 * states the automaton can be in before it has consumed anything
 * (the epsilon closure of the start state). A residue "inserted"
 * before the first pattern position would just make every hit also
 * report a copy starting one residue earlier; for search, where the
 * start position is free anyway, that adds duplicate hits and no
 * information. This is the usual convention for approximate search
 * (semi-global alignment with a free start). */
static void mark_initial_consumes(const Nfa *nfa, bool *initial) {
    int *stack = malloc(sizeof(int) * (size_t)(nfa->count + 1));
    bool *seen = calloc((size_t)nfa->count, sizeof(bool));
    int top = 0;
    stack[top++] = nfa->start;
    seen[nfa->start] = true;
    while (top > 0) {
        int s = stack[--top];
        const NfaState *st = &nfa->states[s];
        if (st->kind == NFA_CONSUME) { initial[s] = true; continue; }
        if (st->kind != NFA_SPLIT) continue;
        int outs[2] = { st->out1, st->out2 };
        for (int o = 0; o < 2; o++) {
            if (outs[o] >= 0 && !seen[outs[o]]) { seen[outs[o]] = true; stack[top++] = outs[o]; }
        }
    }
    free(stack);
    free(seen);
}

void nfa_apply_error_tolerance(const Nfa *plain, int k, ErrorModel model, Nfa *out) {
    nfa_apply_error_tolerance_at(plain, k, model, 0, out);
}

void nfa_apply_error_tolerance_at(const Nfa *plain, int k, ErrorModel model, uint32_t positions, Nfa *out) {
    int n = plain->count;
    bool *initial = calloc((size_t)n, sizeof(bool));
    if (model == ERR_EDIT) mark_initial_consumes(plain, initial);
    int layers = k + 1;

    nfa_init(out, plain->alphabet);
    out->layer_count = layers;

    /* First pass: allocate one state per (layer, old-state) pair, in
     * lockstep, so index (layer*n + j) is reserved for old state j at
     * `layer` before any synthetic states are added. This keeps the
     * "which old state does this belong to" mapping trivial. */
    for (int layer = 0; layer < layers; layer++) {
        for (int j = 0; j < n; j++) {
            NfaKind kind = plain->states[j].kind;
            /* placeholder kind; fixed up in the second pass. CONSUME
             * states get replaced by a SPLIT to two synthetic
             * CONSUME states, added later, so reserve SPLIT here. */
            nfa_add_state(out, kind == NFA_CONSUME ? NFA_SPLIT : kind);
        }
    }
    out->states_per_layer = n;

    #define MAPPED(layer, j) ((layer) * n + (j))

    for (int layer = 0; layer < layers; layer++) {
        for (int j = 0; j < n; j++) {
            const NfaState *src = &plain->states[j];
            int self = MAPPED(layer, j);

            if (src->kind == NFA_SPLIT) {
                out->states[self].out1 = (src->out1 >= 0) ? MAPPED(layer, src->out1) : -1;
                out->states[self].out2 = (src->out2 >= 0) ? MAPPED(layer, src->out2) : -1;

            } else if (src->kind == NFA_MATCH) {
                /* nothing to wire; MATCH has no outgoing edges */

            } else { /* NFA_CONSUME */
                int target = MAPPED(layer, src->out1);

                int normal = nfa_add_state(out, NFA_CONSUME);
                out->states[normal].mask = src->mask;
                out->states[normal].out1 = target;
                out->states[normal].pos = src->pos;

                out->states[self].out1 = normal;

                /* `at P`: errors only at the listed motif positions. */
                bool error_here = positions == 0 ||
                                  (src->pos >= 0 && src->pos < 32 && ((positions >> src->pos) & 1u));

                if (layer < layers - 1 && error_here) {
                    SymMask full = ((SymMask)1u << alphabet_size(plain->alphabet)) - 1u;
                    SymMask complement = full & ~src->mask;
                    int mismatch_target = MAPPED(layer + 1, src->out1);

                    int mism = nfa_add_state(out, NFA_CONSUME);
                    out->states[mism].mask = complement;
                    out->states[mism].out1 = mismatch_target;

                    if (model == ERR_SUBSTITUTION) {
                        out->states[self].out2 = mism;
                    } else {
                        /* self -> {normal, aux}; aux -> {substitution, aux2};
                         * aux2 -> {deletion (epsilon), insertion}. */
                        int aux = nfa_add_state(out, NFA_SPLIT);
                        int aux2 = nfa_add_state(out, NFA_SPLIT);
                        out->states[self].out2 = aux;
                        out->states[aux].out1 = mism;
                        out->states[aux].out2 = aux2;
                        out->states[aux2].out1 = mismatch_target;   /* deletion */
                        if (!initial[j]) {
                            int ins = nfa_add_state(out, NFA_CONSUME);
                            out->states[ins].mask = full;
                            out->states[ins].out1 = MAPPED(layer + 1, j);
                            out->states[aux2].out2 = ins;               /* insertion */
                        } else {
                            out->states[aux2].out2 = -1;
                        }
                    }
                } else {
                    out->states[self].out2 = -1; /* no budget left: no mismatch branch */
                }
            }
        }
    }

    #undef MAPPED

    free(initial);
    out->start = plain->start; /* layer 0's mapping is the identity: 0*n + j == j */
}

static void print_state(const Nfa *nfa, int i) {
    const NfaState *s = &nfa->states[i];
    switch (s->kind) {
        case NFA_SPLIT:
            printf("  %4d  SPLIT    -> %d, %d\n", i, s->out1, s->out2);
            break;
        case NFA_CONSUME: {
            char buf[64];
            alphabet_format_mask(nfa->alphabet, s->mask, buf, sizeof(buf));
            printf("  %4d  CONSUME  [%s] -> %d\n", i, buf, s->out1);
            break;
        }
        case NFA_MATCH:
            printf("  %4d  MATCH\n", i);
            break;
    }
}

void nfa_print(const Nfa *nfa) {
    printf("== NFA (%d states, start=%d", nfa->count, nfa->start);
    if (nfa->layer_count > 1) {
        printf(", %d error layers x %d states", nfa->layer_count, nfa->states_per_layer);
    }
    printf(") ==\n");
    for (int i = 0; i < nfa->count; i++) print_state(nfa, i);
}
