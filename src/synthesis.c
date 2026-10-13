#include "synthesis.h"
#include "sema.h"

/* One choice for one free motif: a budget and (for `at ?`) positions.
 * mask 0 means "every position". */
typedef struct { int k; uint32_t mask; } Option;

typedef struct {
    int opt[MAX_FREE_BUDGETS];    /* index into each motif's option list */
    int ksum, psum;
    int state;                    /* 0 pruned, 1 feasible, 2 infeasible (compiled) */
    double rate;                  /* chance hits per position, if compiled */
    bool have_sens;
    SensitivityResult sens;
} Candidate;

typedef struct {
    int motif;
    int length;                   /* residues, for `at ?` motifs; 0 otherwise */
    Option *opts;
    int count;
} FreeMotif;

static int popcount32(uint32_t x) { int c = 0; for (; x; x >>= 1) c += (int)(x & 1u); return c; }

static uint32_t full_mask(int length) { return length >= 32 ? 0xFFFFFFFFu : (((uint32_t)1u << length) - 1u); }

/* a >= b for one motif's options. With k = 0 the positions do not matter. */
static bool opt_above(const FreeMotif *f, Option a, Option b) {
    if (a.k < b.k) return false;
    if (b.k == 0) return true;
    uint32_t ma = a.mask ? a.mask : full_mask(f->length ? f->length : 32);
    uint32_t mb = b.mask ? b.mask : full_mask(f->length ? f->length : 32);
    return (ma & mb) == mb;
}

static bool dominates(const FreeMotif *fm, int n, const Candidate *a, const Candidate *b) {
    for (int i = 0; i < n; i++)
        if (!opt_above(&fm[i], fm[i].opts[a->opt[i]], fm[i].opts[b->opt[i]])) return false;
    return true;
}

static int cand_cmp(const void *pa, const void *pb) {
    const Candidate *a = pa, *b = pb;
    if (a->ksum != b->ksum) return a->ksum - b->ksum;
    return a->psum - b->psum;
}

static void apply(Program *program, const FreeMotif *fm, int n, const Candidate *c) {
    for (int i = 0; i < n; i++) {
        MotifDef *m = &program->motifs[fm[i].motif];
        Option o = fm[i].opts[c->opt[i]];
        m->mismatch_budget = o.k;
        m->error_positions = o.mask;
    }
    sema_recompute_lengths(program);
}

static void apply_zero(Program *program, const FreeMotif *fm, int n) {
    for (int i = 0; i < n; i++) {
        program->motifs[fm[i].motif].mismatch_budget = 0;
        program->motifs[fm[i].motif].error_positions = 0;
    }
    sema_recompute_lengths(program);
}

static double evaluate_rate(Program *program, int motif, const RateSource *src) {
    MatcherProgram mp;
    compile_motif_to_bytecode(program, motif, &mp, NULL);
    double q = rate_of(&mp, src);
    matcher_program_free(&mp);
    return q;
}

static SensitivityResult evaluate_sens(Program *program, int motif, const MutationModel *model,
                                       const SynthesisOptions *opts) {
    MatcherProgram mp;
    compile_motif_to_bytecode(program, motif, &mp, NULL);
    SensitivityResult s = opts->sites ? sensitivity_on_sites(&mp, opts->sites)
                                      : sensitivity_estimate(program, motif, &mp, model, opts->samples, opts->seed);
    matcher_program_free(&mp);
    return s;
}

static void option_text(const MotifDef *m, Option o, char *buf, size_t n) {
    int len = snprintf(buf, n, "%d %s", o.k,
                       m->error_model == ERR_EDIT ? (o.k == 1 ? "edit" : "edits")
                                                  : (o.k == 1 ? "mismatch" : "mismatches"));
    if (o.k > 0 && o.mask) {
        len += snprintf(buf + len, n - (size_t)len, " at ");
        bool first = true;
        for (int p = 0; p < 32 && (size_t)len < n; p++) {
            if (!((o.mask >> p) & 1u)) continue;
            len += snprintf(buf + len, n - (size_t)len, first ? "%d" : ",%d", p + 1);
            first = false;
        }
    }
}

static bool build_options(Program *program, int i, FreeMotif *f) {
    MotifDef *m = &program->motifs[i];
    f->motif = i;
    f->length = 0;
    int base = m->min_len;
    int cap = (m->error_model == ERR_EDIT) ? (base - 1) / 2 : base - 1;
    if (cap > MAX_SYNTH_BUDGET) cap = MAX_SYNTH_BUDGET;
    if (cap < 0) cap = 0;

    if (!m->positions_free) {
        f->opts = calloc((size_t)cap + 1, sizeof(Option));
        f->count = 0;
        for (int k = 0; k <= cap; k++) f->opts[f->count++] = (Option){ k, 0 };
        return true;
    }

    int L = sema_linear_length(m);
    if (L <= 0 || L > 16) return false;            /* 2^16 subsets is the practical limit */
    f->length = L;
    /* A position subset smaller than the motif can never make it degenerate,
     * so the budget may go up to the number of positions chosen. */
    int kcap = MAX_SYNTH_BUDGET < L ? MAX_SYNTH_BUDGET : L;
    size_t max_opts = 1 + (size_t)kcap * ((size_t)1 << L);
    f->opts = calloc(max_opts, sizeof(Option));
    f->count = 0;
    f->opts[f->count++] = (Option){ 0, 0 };
    for (int k = 1; k <= kcap; k++) {
        for (uint32_t mask = 1; mask <= full_mask(L); mask++) {
            int pc = popcount32(mask);
            if (pc < k) continue;                                   /* k beyond |positions| adds nothing */
            if (pc == L && k > cap) continue;                       /* all positions: usual cap applies */
            f->opts[f->count++] = (Option){ k, mask };
        }
    }
    return true;
}

static bool synthesize_one(Program *program, Diagnostics *diag, int r, const SynthesisOptions *opts) {
    MotifDef *req = &program->motifs[r];
    MutationModel model = opts->model;
    if (req->divergence >= 0) model.substitution = req->divergence;
    double unit_len = req->required_rate_shown / req->required_rate;   /* positions per contract unit */

    FreeMotif fm[MAX_FREE_BUDGETS];
    int n = 0;
    for (int i = 0; i < program->motif_count; i++) {
        if (program->motifs[i].mismatch_budget != BUDGET_FREE) continue;
        if (sema_governing_requirement(program, i) != r) continue;
        if (n == MAX_FREE_BUDGETS) {
            diag_reportf(diag, PHASE_SPECIFICITY, SEV_ERROR, req->require_line, req->require_column,
                         "more than %d free budgets under one 'require' clause", MAX_FREE_BUDGETS);
            return false;
        }
        if (!build_options(program, i, &fm[n])) {
            diag_reportf(diag, PHASE_SPECIFICITY, SEV_ERROR, program->motifs[i].line, program->motifs[i].column,
                         "'at ?' supports fixed motifs of up to 16 residues");
            return false;
        }
        n++;
    }

    char srcbuf[256];
    printf("== budget synthesis: %.*s (require rate <= %g per %s",
           req->name_length, req->name, req->required_rate_shown, req->required_unit);
    if (req->divergence >= 0 && !opts->sites) printf(", divergence %g%%", req->divergence * 100);
    printf(") ==\n");
    printf("chance hits: %s\n", rate_source_text(&opts->rates, srcbuf, sizeof(srcbuf)));
    if (opts->sites)
        printf("sensitivity: share of %d known sites in %s that are found\n", opts->sites->count, opts->sites_label);
    else
        printf("sensitivity: sites mutated with %g%% substitution, %g%% indel; %d sampled sites, seed %llu\n",
               model.substitution * 100, model.indel * 100, opts->samples, (unsigned long long)opts->seed);

    if (n == 0) {
        double q = evaluate_rate(program, r, &opts->rates);
        SensitivityResult s = evaluate_sens(program, r, &model, opts);
        printf("no free budgets: checking the contract as written\n");
        printf("  expected %.3g per %s (limit %g), sensitivity %.2f [%.2f, %.2f]\n\n",
               q * unit_len, req->required_unit, req->required_rate_shown, s.sensitivity, s.ci_low, s.ci_high);
        if (q > req->required_rate) {
            diag_reportf(diag, PHASE_SPECIFICITY, SEV_ERROR, req->require_line, req->require_column,
                         "motif '%.*s' expects %.3g chance hits per %s, above the required %g",
                         req->name_length, req->name, q * unit_len, req->required_unit, req->required_rate_shown);
            return false;
        }
        return true;
    }

    /* Every combination of options, smallest first. */
    long total = 1;
    for (int i = 0; i < n; i++) {
        total *= fm[i].count;
        if (total > MAX_CANDIDATES) {
            diag_reportf(diag, PHASE_SPECIFICITY, SEV_ERROR, req->require_line, req->require_column,
                         "too many budget combinations to search (more than %d); fix some budgets or positions",
                         MAX_CANDIDATES);
            for (int j = 0; j < n; j++) free(fm[j].opts);
            return false;
        }
    }
    Candidate *c = calloc((size_t)total, sizeof(Candidate));
    for (long idx = 0; idx < total; idx++) {
        long rem = idx;
        for (int i = 0; i < n; i++) {
            c[idx].opt[i] = (int)(rem % fm[i].count);
            rem /= fm[i].count;
            Option o = fm[i].opts[c[idx].opt[i]];
            c[idx].ksum += o.k;
            c[idx].psum += o.k == 0 ? 0 : (o.mask ? popcount32(o.mask) : 32);
        }
    }
    qsort(c, (size_t)total, sizeof(Candidate), cand_cmp);

    int compiled = 0, pruned = 0;
    int *infeasible = malloc(sizeof(int) * (size_t)total);
    int ninf = 0;
    for (long idx = 0; idx < total; idx++) {
        bool skip = false;
        for (int j = 0; j < ninf && !skip; j++)
            if (dominates(fm, n, &c[idx], &c[infeasible[j]])) skip = true;
        if (skip) { c[idx].state = 0; pruned++; continue; }
        apply(program, fm, n, &c[idx]);
        c[idx].rate = evaluate_rate(program, r, &opts->rates);
        if (c[idx].rate <= req->required_rate) c[idx].state = 1;
        else { c[idx].state = 2; infeasible[ninf++] = (int)idx; }
        compiled++;
    }

    /* Maximal feasible: no other feasible candidate lies above it. */
    int best = -1, maximal = 0;
    for (long idx = 0; idx < total; idx++) {
        if (c[idx].state != 1) continue;
        bool is_max = true;
        for (long j = 0; j < total && is_max; j++)
            if (j != idx && c[j].state == 1 && c[j].ksum + c[j].psum > c[idx].ksum + c[idx].psum &&
                dominates(fm, n, &c[j], &c[idx])) is_max = false;
        if (!is_max) continue;
        maximal++;
        apply(program, fm, n, &c[idx]);
        c[idx].sens = evaluate_sens(program, r, &model, opts);
        c[idx].have_sens = true;
        if (best < 0 || c[idx].sens.sensitivity > c[best].sens.sensitivity ||
            (c[idx].sens.sensitivity == c[best].sens.sensitivity && c[idx].rate < c[best].rate))
            best = (int)idx;
    }
    /* Sensitivity is the first goal, the chance-hit rate the second:
     * among feasible assignments as sensitive as the best maximal one,
     * keep the one with the fewest chance hits. With a finite set of
     * known sites sensitivity saturates (all sites found), and a much
     * tighter pattern often finds them all too. */
    if (best >= 0) {
        double target_sens = c[best].sens.sensitivity;
        for (;;) {
            long next = -1;
            for (long idx = 0; idx < total; idx++)
                if (c[idx].state == 1 && c[idx].rate < c[best].rate && !c[idx].have_sens &&
                    (next < 0 || c[idx].rate < c[next].rate)) next = idx;
            if (next < 0) break;
            apply(program, fm, n, &c[next]);
            c[next].sens = evaluate_sens(program, r, &model, opts);
            c[next].have_sens = true;
            if (c[next].sens.sensitivity >= target_sens) { best = (int)next; break; }
        }
    }

    for (int j = 0; j < ninf; j++) {                 /* minimal infeasible: what one more step costs */
        apply(program, fm, n, &c[infeasible[j]]);
        c[infeasible[j]].sens = evaluate_sens(program, r, &model, opts);
        c[infeasible[j]].have_sens = true;
    }

    printf("free budgets:");
    for (int i = 0; i < n; i++) {
        MotifDef *m = &program->motifs[fm[i].motif];
        printf(" %.*s (%s%s, %d options)%s", m->name_length, m->name,
               m->error_model == ERR_EDIT ? "edits" : "mismatches",
               m->positions_free ? ", positions free" : "", fm[i].count, i + 1 < n ? "," : "");
    }
    printf("\nassignments: %ld possible, %d compiled, %d skipped by monotone pruning, %d maximal feasible\n",
           total, compiled, pruned, maximal);

    /* Show the chosen assignment, the best other maximal feasible ones,
     * and the most sensitive infeasible ones (at most 5 of each). */
    long *order = malloc(sizeof(long) * (size_t)total);
    printf("  %-44s %-12s %-22s %s\n", "budgets", "chance/unit", "sensitivity [95% CI]", "status");
    for (int pass = 0; pass < 2; pass++) {
        long m = 0;
        for (long idx = 0; idx < total; idx++)
            if (c[idx].have_sens && (pass == 0) == (c[idx].state == 1)) order[m++] = idx;
        /* chosen first, then by sensitivity, highest first (simple selection; m is small) */
        for (long a2 = 0; a2 < m; a2++) {
            long top = a2;
            for (long b2 = a2 + 1; b2 < m; b2++) {
                bool better = (order[b2] == best) ||
                              (order[top] != best && c[order[b2]].sens.sensitivity > c[order[top]].sens.sensitivity);
                if (better) top = b2;
            }
            long t = order[a2]; order[a2] = order[top]; order[top] = t;
        }
        for (long a2 = 0; a2 < m && a2 < 5; a2++) {
            long pick = order[a2];
            char opts_txt[256] = "";
            size_t len = 0;
            for (int i = 0; i < n; i++) {
                MotifDef *md = &program->motifs[fm[i].motif];
                char one[96];
                option_text(md, fm[i].opts[c[pick].opt[i]], one, sizeof(one));
                len += (size_t)snprintf(opts_txt + len, sizeof(opts_txt) - len, "%s%.*s %s",
                                        i ? ", " : "", md->name_length > 8 ? 8 : md->name_length, md->name, one);
            }
            char sens[48];
            snprintf(sens, sizeof(sens), "%.3f [%.3f,%.3f]", c[pick].sens.sensitivity, c[pick].sens.ci_low, c[pick].sens.ci_high);
            printf("  %-44s %-12.3g %-22s %s\n", opts_txt, c[pick].rate * unit_len, sens,
                   pick == best ? "CHOSEN" : (c[pick].state == 1 ? "feasible, maximal" : "exceeds the limit"));
        }
        if (m > 5) printf("  ... %ld more %s\n", m - 5, pass == 0 ? "feasible" : "minimal infeasible");
    }
    free(order);
    printf("(chance/unit = expected chance hits per %s)\n", req->required_unit);

    bool ok = best >= 0;
    if (ok) {
        apply(program, fm, n, &c[best]);
        printf("note: chose");
        for (int i = 0; i < n; i++) {
            MotifDef *m = &program->motifs[fm[i].motif];
            char one[96];
            option_text(m, fm[i].opts[c[best].opt[i]], one, sizeof(one));
            printf(" %.*s = %s%s", m->name_length, m->name, one, i + 1 < n ? "," : "");
        }
        printf("\n\n");
    } else {
        apply_zero(program, fm, n);
        diag_reportf(diag, PHASE_SPECIFICITY, SEV_ERROR, req->require_line, req->require_column,
                     "no budget assignment meets 'rate <= %g per %s': even with no errors allowed "
                     "motif '%.*s' expects %.3g",
                     req->required_rate_shown, req->required_unit, req->name_length, req->name,
                     evaluate_rate(program, r, &opts->rates) * unit_len);
        printf("\n");
    }
    for (int i = 0; i < n; i++) free(fm[i].opts);
    free(infeasible);
    free(c);
    return ok;
}

bool synthesize_budgets(Program *program, Diagnostics *diag, const SynthesisOptions *opts) {
    bool ok = true;
    for (int r = 0; r < program->motif_count; r++) {
        if (!program->motifs[r].has_requirement) continue;
        if (!synthesize_one(program, diag, r, opts)) ok = false;
    }
    return ok;
}
