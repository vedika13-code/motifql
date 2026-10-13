#include "sema.h"
#include <stdarg.h>

typedef struct {
    Program *program;
    Diagnostics *diag;
    bool had_error;
} Sema;

static void sema_error(Sema *s, int line, int column, const char *fmt, ...) {
    char buffer[512];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buffer, sizeof(buffer), fmt, args);
    va_end(args);
    diag_report(s->diag, PHASE_SEMANTIC, SEV_ERROR, line, column, buffer);
    s->had_error = true;
}

static void sema_warn(Sema *s, int line, int column, const char *fmt, ...) {
    char buffer[512];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buffer, sizeof(buffer), fmt, args);
    va_end(args);
    diag_report(s->diag, PHASE_SEMANTIC, SEV_WARNING, line, column, buffer);
}

/* ---- pass 1: symbol resolution ------------------------------------------- */

static void resolve_refs(Sema *s, Node *node) {
    if (!node) return;

    if (node->type == NODE_REF && node->text) {
        int index = program_find_motif(s->program, node->text, node->text_length);
        if (index < 0) {
            sema_error(s, node->line, node->column,
                       "reference to undefined motif '@%.*s'",
                       node->text_length, node->text);
        } else {
            node->ref_index = index;
            s->program->motifs[index].used = true;
        }
    }

    for (int i = 0; i < node->child_count; i++) resolve_refs(s, node->children[i]);
}

static void resolve_all_symbols(Sema *s) {
    for (int i = 0; i < s->program->motif_count; i++) {
        resolve_refs(s, s->program->motifs[i].pattern);
    }

    for (int i = 0; i < s->program->search_count; i++) {
        SearchStmt *stmt = &s->program->searches[i];
        int index = program_find_motif(s->program, stmt->motif_name, stmt->motif_name_length);
        if (index < 0) {
            sema_error(s, stmt->line, stmt->column,
                       "search targets undefined motif '%.*s'",
                       stmt->motif_name_length, stmt->motif_name);
        } else {
            stmt->motif_index = index;
            s->program->motifs[index].used = true;
        }
    }
}

/* ---- pass 2: reference cycle detection -----------------------------------
 *
 * A motif referencing itself, directly or through a chain, describes a
 * pattern of unbounded nesting depth. A regular language cannot express
 * that, so it cannot be compiled to a finite automaton in Phase 2 --
 * this is exactly the boundary between regular and context-free
 * grammars, and rejecting it here is what keeps Phase 2's Thompson
 * construction guaranteed to terminate.
 *
 * Standard three-colour DFS: 0 = unvisited, 1 = on the current stack,
 * 2 = fully explored. Reaching a node coloured 1 means a back edge,
 * hence a cycle.
 */

static bool detect_cycle_in(Sema *s, int motif_index, Node *node) {
    if (!node) return false;

    if (node->type == NODE_REF && node->ref_index >= 0) {
        MotifDef *target = &s->program->motifs[node->ref_index];
        if (target->visit_state == 1) {
            sema_error(s, node->line, node->column,
                       "circular motif reference: '@%.*s' is already being defined "
                       "(a motif cannot refer to itself, directly or indirectly)",
                       node->text_length, node->text);
            return true;
        }
        if (target->visit_state == 0) {
            target->visit_state = 1;
            bool cyclic = detect_cycle_in(s, node->ref_index, target->pattern);
            target->visit_state = 2;
            if (cyclic) return true;
        }
    }

    for (int i = 0; i < node->child_count; i++) {
        if (detect_cycle_in(s, motif_index, node->children[i])) return true;
    }
    return false;
}

static bool detect_cycles(Sema *s) {
    for (int i = 0; i < s->program->motif_count; i++) {
        s->program->motifs[i].visit_state = 0;
    }
    bool any = false;
    for (int i = 0; i < s->program->motif_count; i++) {
        MotifDef *m = &s->program->motifs[i];
        if (m->visit_state != 0) continue;
        m->visit_state = 1;
        if (detect_cycle_in(s, i, m->pattern)) any = true;
        m->visit_state = 2;
    }
    return !any;
}

/* ---- pass 3: alphabet checking -------------------------------------------- */

static void check_alphabet(Sema *s, Node *node) {
    if (!node) return;
    AlphabetKind kind = s->program->alphabet;

    if (node->type == NODE_LITERAL && node->text) {
        node->mask_count = node->text_length;
        node->masks = calloc((size_t)(node->text_length > 0 ? node->text_length : 1), sizeof(SymMask));
        if (!node->masks) { fprintf(stderr, "out of memory\n"); exit(1); }

        for (int i = 0; i < node->text_length; i++) {
            char c = node->text[i];
            SymMask mask = alphabet_symbol_mask(kind, c);
            if (mask == 0) {
                sema_error(s, node->line, node->column + 1 + i,
                           "'%c' is not a valid symbol in alphabet %s",
                           c, alphabet_name(kind));
            }
            node->masks[i] = mask;
        }
    }

    if (node->type == NODE_CLASS && node->text) {
        SymMask mask = 0;
        for (int i = 0; i < node->text_length; i++) {
            char c = node->text[i];
            SymMask bit = alphabet_symbol_mask(kind, c);
            if (bit == 0) {
                sema_error(s, node->line, node->column + 1 + i,
                           "'%c' is not a valid symbol in alphabet %s",
                           c, alphabet_name(kind));
            }
            mask |= bit;
        }
        if (node->negated) {
            int n = alphabet_size(kind);
            SymMask full = (n >= 32) ? (SymMask)0xFFFFFFFFu : (((SymMask)1u << n) - 1u);
            mask = full & ~mask;
            if (mask == 0) {
                sema_error(s, node->line, node->column,
                           "negated class excludes every symbol in alphabet %s "
                           "and can never match", alphabet_name(kind));
            }
        }
        node->class_mask = mask;
    }

    for (int i = 0; i < node->child_count; i++) check_alphabet(s, node->children[i]);
}

/* ---- pass 4: quantifier bound checking ------------------------------------ */

static void check_quantifiers(Sema *s, Node *node) {
    if (!node) return;

    if (node->type == NODE_REPEAT) {
        if (node->repeat_max != REPEAT_INF && node->repeat_min > node->repeat_max) {
            sema_error(s, node->line, node->column,
                       "quantifier lower bound %d exceeds upper bound %d",
                       node->repeat_min, node->repeat_max);
            /* Normalise so length inference downstream doesn't produce
             * an impossible range like 15..6 in the symbol table dump.
             * The program is already rejected; this only keeps the
             * diagnostic output honest. */
            node->repeat_max = node->repeat_min;
        }
        if (node->repeat_max == 0) {
            sema_error(s, node->line, node->column,
                       "quantifier {0,0} makes this element unmatchable");
        }
    }

    for (int i = 0; i < node->child_count; i++) check_quantifiers(s, node->children[i]);
}

/* ---- pass 5: match-length inference ---------------------------------------
 *
 * Computed bottom-up over the AST. Beyond diagnostics, Phase 2 needs
 * these numbers: the minimum length bounds how far a scan window can
 * be advanced on a failed match, and a finite maximum tells the
 * matcher how much lookahead to buffer over a FASTA stream.
 *
 * Saturating arithmetic keeps unbounded ranges (from '*' and '+')
 * from overflowing into meaningless finite numbers.
 */

static int add_len(int a, int b) {
    if (a == LEN_INF || b == LEN_INF) return LEN_INF;
    return a + b;
}

static int mul_len(int a, int count) {
    if (a == LEN_INF) return (count == 0) ? 0 : LEN_INF;
    if (count == REPEAT_INF) return (a == 0) ? 0 : LEN_INF;
    return a * count;
}

static void infer_motif_lengths(Sema *s, int motif_index);

static void infer_lengths(Sema *s, Node *node) {
    if (!node) return;
    for (int i = 0; i < node->child_count; i++) infer_lengths(s, node->children[i]);

    switch (node->type) {
        case NODE_LITERAL:
            node->min_len = node->mask_count;
            node->max_len = node->mask_count;
            break;

        case NODE_CLASS:
            node->min_len = 1;
            node->max_len = 1;
            break;

        case NODE_CONCAT: {
            int mn = 0, mx = 0;
            for (int i = 0; i < node->child_count; i++) {
                mn = add_len(mn, node->children[i]->min_len);
                mx = add_len(mx, node->children[i]->max_len);
            }
            node->min_len = mn;
            node->max_len = mx;
            break;
        }

        case NODE_ALTERNATE: {
            if (node->child_count == 0) { node->min_len = 0; node->max_len = 0; break; }
            int mn = node->children[0]->min_len;
            int mx = node->children[0]->max_len;
            for (int i = 1; i < node->child_count; i++) {
                int cmin = node->children[i]->min_len;
                int cmax = node->children[i]->max_len;
                if (cmin < mn) mn = cmin;
                if (mx != LEN_INF && (cmax == LEN_INF || cmax > mx)) mx = cmax;
            }
            node->min_len = mn;
            node->max_len = mx;
            break;
        }

        case NODE_REPEAT: {
            Node *child = node->child_count > 0 ? node->children[0] : NULL;
            int cmin = child ? child->min_len : 0;
            int cmax = child ? child->max_len : 0;
            node->min_len = mul_len(cmin, node->repeat_min);
            node->max_len = mul_len(cmax, node->repeat_max);
            break;
        }

        case NODE_REF: {
            if (node->ref_index < 0) { node->min_len = 0; node->max_len = 0; break; }
            infer_motif_lengths(s, node->ref_index);
            MotifDef *target = &s->program->motifs[node->ref_index];
            node->min_len = target->min_len;
            node->max_len = target->max_len;
            /* A referenced motif's own edit budget lets its span shrink
             * or grow by up to k inside the composite. The composite's
             * length bounds must include that, or the matcher's length
             * cut-off would stop reading before an insertion finished. */
            if (target->error_model == ERR_EDIT && target->mismatch_budget > 0) {
                int k = target->mismatch_budget;
                node->min_len = node->min_len > k ? node->min_len - k : 0;
                if (node->max_len != LEN_INF) node->max_len += k;
            }
            break;
        }
    }
}

/* Memoised via visit_state so a motif referenced from several places
 * is only measured once. Safe to recurse because pass 2 already
 * proved the reference graph is acyclic. */
static void infer_motif_lengths(Sema *s, int motif_index) {
    MotifDef *m = &s->program->motifs[motif_index];
    if (m->visit_state == 2) return;
    m->visit_state = 2;

    infer_lengths(s, m->pattern);
    m->min_len = m->pattern ? m->pattern->min_len : 0;
    m->max_len = m->pattern ? m->pattern->max_len : 0;
}

static void infer_all_lengths(Sema *s) {
    for (int i = 0; i < s->program->motif_count; i++) {
        s->program->motifs[i].visit_state = 0;
    }
    for (int i = 0; i < s->program->motif_count; i++) {
        infer_motif_lengths(s, i);
    }
}

/* Re-runs length inference after budgets change (budget synthesis
 * tries many assignments). No diagnostics are produced. */
void sema_recompute_lengths(Program *program) {
    Sema s;
    s.program = program;
    s.diag = NULL;
    s.had_error = false;
    infer_all_lengths(&s);
}

/* ---- pass 5b: checks that depend on inferred lengths ---------------------- */

static void check_length_derived_rules(Sema *s) {
    for (int i = 0; i < s->program->motif_count; i++) {
        MotifDef *m = &s->program->motifs[i];

        /* A motif that can match zero residues matches at every
         * position in the sequence, which is never what anyone means
         * and would make the Phase 2 scan degenerate. */
        if (m->min_len == 0) {
            sema_error(s, m->line, m->column,
                       "motif '%.*s' can match the empty sequence "
                       "(every quantifier on its top level allows zero repetitions)",
                       m->name_length, m->name);
        }

        /* Allowing as many mismatches as the motif is long means every
         * window matches. Comparing against min_len is the strict
         * reading: if the shortest alternative can be fully mismatched,
         * the budget is too large. */
        /* With `at P`, at most |P| positions can be in error, so the
         * effective budget for the "too loose" checks is the smaller. */
        int saved_budget = m->mismatch_budget;
        if (m->error_positions && m->mismatch_budget > 0) {
            int npos = 0;
            for (uint32_t b = m->error_positions; b; b >>= 1) npos += (int)(b & 1u);
            if (npos < m->min_len) m->mismatch_budget = 0;   /* some positions stay exact: never degenerate */
        }
        if (m->mismatch_budget < 0 && m->mismatch_budget != BUDGET_FREE) {
            sema_error(s, m->line, m->column,
                       "mismatch allowance cannot be negative");
        } else if (m->error_model == ERR_SUBSTITUTION && m->min_len != LEN_INF &&
                   m->mismatch_budget >= m->min_len && m->min_len > 0) {
            sema_error(s, m->line, m->column,
                       "motif '%.*s' allows %d mismatch(es) but its shortest match is "
                       "only %d residue(s) long, so every position would match",
                       m->name_length, m->name, m->mismatch_budget, m->min_len);
        } else if (m->error_model == ERR_EDIT && m->min_len != LEN_INF &&
                   m->mismatch_budget >= m->min_len && m->min_len > 0) {
            /* Under edit distance the budget can be spent on deletions,
             * and deleting every position of the shortest alternative
             * leaves the empty string -- the same degenerate "matches
             * everywhere" motif the empty-match rule above rejects. */
            sema_error(s, m->line, m->column,
                       "motif '%.*s' allows %d edit(s) but its shortest match is only "
                       "%d residue(s) long, so deletions alone could match the empty sequence",
                       m->name_length, m->name, m->mismatch_budget, m->min_len);
        } else if (m->error_model == ERR_EDIT && m->min_len != LEN_INF &&
                   m->min_len > 0 && 2 * m->mismatch_budget >= m->min_len) {
            /* Not wrong, but after k deletions fewer than k residues of
             * the motif are left to pin the match down. */
            sema_warn(s, m->line, m->column,
                      "motif '%.*s' allows %d edit(s) against a shortest match of %d residue(s); "
                      "after deletions only %d residue(s) would remain, so expect many chance hits",
                      m->name_length, m->name, m->mismatch_budget, m->min_len,
                      m->min_len - m->mismatch_budget);
        }

        m->mismatch_budget = saved_budget;

        if (!m->used) {
            sema_warn(s, m->line, m->column,
                      "motif '%.*s' is never referenced or searched",
                      m->name_length, m->name);
        }
    }
}

/* ---- pass 5c: position-restricted errors ---------------------------------- */

/* `at P` names residue positions, which only makes sense for a motif
 * that is a fixed run of residues: literals and classes, concatenated.
 * Alternation, repetition or references would make "position 3"
 * ambiguous. Returns the motif's length, or -1 if it is not linear. */
static int linear_length(const Node *node) {
    if (!node) return -1;
    switch (node->type) {
        case NODE_LITERAL: return node->mask_count;
        case NODE_CLASS:   return 1;
        case NODE_CONCAT: {
            int total = 0;
            for (int i = 0; i < node->child_count; i++) {
                int l = linear_length(node->children[i]);
                if (l < 0) return -1;
                total += l;
            }
            return total;
        }
        default: return -1;
    }
}

int sema_linear_length(const MotifDef *m) { return linear_length(m->pattern); }

static void check_positions(Sema *s) {
    for (int i = 0; i < s->program->motif_count; i++) {
        MotifDef *m = &s->program->motifs[i];
        if (!m->error_positions && !m->positions_free) continue;
        int len = linear_length(m->pattern);
        if (len < 0) {
            sema_error(s, m->line, m->column,
                       "'at' needs motif '%.*s' to be a fixed run of residues (no |, repeats or @references), "
                       "so that positions are well defined", m->name_length, m->name);
        } else if (len > 32) {
            sema_error(s, m->line, m->column, "'at' supports motifs of up to 32 residues");
        } else if (m->error_positions >> len) {
            sema_error(s, m->line, m->column,
                       "motif '%.*s' has %d residues, but 'at' names a position beyond that",
                       m->name_length, m->name, len);
        }
        if (m->positions_free && m->mismatch_budget != BUDGET_FREE) {
            sema_error(s, m->line, m->column, "'at ?' needs a free budget too: write 'allow ? ... at ?'");
        }
    }
}

/* ---- pass 6: budget-synthesis contracts ----------------------------------- */

/* Marks every motif reachable from `index` through @references,
 * including `index` itself. The reference graph is acyclic by now. */
static void mark_reachable_node(Program *program, Node *node, bool *mark);
static void mark_reachable(Program *program, int index, bool *mark) {
    if (mark[index]) return;
    mark[index] = true;
    mark_reachable_node(program, program->motifs[index].pattern, mark);
}
static void mark_reachable_node(Program *program, Node *node, bool *mark) {
    if (!node) return;
    if (node->type == NODE_REF && node->ref_index >= 0) mark_reachable(program, node->ref_index, mark);
    for (int i = 0; i < node->child_count; i++) mark_reachable_node(program, node->children[i], mark);
}

int sema_governing_requirement(Program *program, int free_index) {
    int found = -1;
    for (int r = 0; r < program->motif_count; r++) {
        if (!program->motifs[r].has_requirement) continue;
        bool mark[MAX_MOTIFS] = {false};
        mark_reachable(program, r, mark);
        if (mark[free_index]) {
            if (found >= 0) return -2;   /* governed by two contracts: ambiguous */
            found = r;
        }
    }
    return found;
}

/* A free budget (`allow ?`) must be settled by exactly one `require`
 * clause: the motif's own, or that of the one motif that uses it.
 * With none, nothing says what to aim for; with two, the contracts
 * could ask for different budgets for the same motif. */
static void check_requirements(Sema *s) {
    Program *program = s->program;
    for (int i = 0; i < program->motif_count; i++) {
        MotifDef *m = &program->motifs[i];
        if (m->has_requirement) {
            if (m->required_rate <= 0) {
                sema_error(s, m->require_line, m->require_column,
                           "a required rate of 0 can never be met; use a positive number of hits");
            }
            if (m->divergence >= 1.0) {
                sema_error(s, m->require_line, m->require_column,
                           "divergence must be below 100%%");
            }
        }
        if (m->mismatch_budget != BUDGET_FREE) continue;
        int gov = sema_governing_requirement(program, i);
        if (gov == -1) {
            sema_error(s, m->line, m->column,
                       "motif '%.*s' has 'allow ?' but no 'require' clause uses it, so there is "
                       "nothing to choose its budget against",
                       m->name_length, m->name);
        } else if (gov == -2) {
            sema_error(s, m->line, m->column,
                       "motif '%.*s' has 'allow ?' and is used by more than one motif with a "
                       "'require' clause; its budget would be ambiguous",
                       m->name_length, m->name);
        }
    }
}

/* ---- driver ---------------------------------------------------------------- */

bool analyze_program(Program *program, Diagnostics *diag) {
    Sema s;
    s.program = program;
    s.diag = diag;
    s.had_error = false;

    if (!program->alphabet_declared) {
        sema_warn(&s, 1, 1,
                  "no 'alphabet' declaration; defaulting to DNA");
    }

    resolve_all_symbols(&s);

    /* Cycle detection must succeed before anything recurses through
     * references, or length inference would not terminate. */
    if (!detect_cycles(&s)) return false;

    for (int i = 0; i < program->motif_count; i++) {
        check_alphabet(&s, program->motifs[i].pattern);
        check_quantifiers(&s, program->motifs[i].pattern);
    }

    infer_all_lengths(&s);
    check_length_derived_rules(&s);
    check_positions(&s);
    check_requirements(&s);

    return !s.had_error;
}
