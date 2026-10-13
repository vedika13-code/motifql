#include "specificity.h"
#include "matcher.h"

static int residue_index(AlphabetKind kind, char c) {
    const char *residues = alphabet_residues(kind);
    const char *p = strchr(residues, c);
    return (p && c) ? (int)(p - residues) : -1;
}

void background_free(Background *bg) {
    free(bg->init);
    free(bg->trans);
    bg->init = bg->trans = NULL;
}

void background_uniform(Background *bg, AlphabetKind alphabet) {
    memset(bg, 0, sizeof(*bg));
    bg->alphabet = alphabet;
    bg->size = alphabet_size(alphabet);
    bg->order = 0;
    bg->contexts = 1;
    bg->init = calloc(1, sizeof(double));
    bg->trans = calloc((size_t)bg->size, sizeof(double));
    bg->init[0] = 1.0;
    for (int r = 0; r < bg->size; r++) bg->p0[r] = bg->trans[r] = 1.0 / bg->size;
    snprintf(bg->source, sizeof(bg->source), "uniform %s", alphabet_name(alphabet));
}

bool background_from_fasta(Background *bg, AlphabetKind alphabet, const FastaFile *f, int order,
                           const char *label) {
    memset(bg, 0, sizeof(*bg));
    bg->alphabet = alphabet;
    bg->size = alphabet_size(alphabet);
    int max_order = (bg->size > 4) ? 1 : BG_MAX_ORDER;
    if (order < 0) order = 0;
    if (order > max_order) order = max_order;
    bg->order = order;
    int A = bg->size, C = 1;
    for (int i = 0; i < order; i++) C *= A;
    bg->contexts = C;

    double *ctx_count = calloc((size_t)C, sizeof(double));
    double *pair_count = calloc((size_t)C * A, sizeof(double));
    double c0[BG_MAX] = {0}, total = 0;

    for (int i = 0; i < f->count; i++) {
        int ctx = 0, have = 0;   /* rolling context of the last `have` (<= order) residues */
        for (long j = 0; j < f->records[i].length; j++) {
            int r = residue_index(alphabet, f->records[i].sequence[j]);
            if (r < 0) { ctx = 0; have = 0; continue; }
            c0[r] += 1; total += 1;
            if (have == order) {
                ctx_count[ctx] += 1;
                pair_count[(size_t)ctx * A + r] += 1;
            }
            if (C > 1) ctx = (ctx * A + r) % C;
            if (have < order) have++;
        }
    }
    if (total == 0) { free(ctx_count); free(pair_count); return false; }

    /* Pseudocounts of 1 keep every transition possible, so a context
     * that happens to be absent from a short file is rare, not impossible. */
    bg->init = calloc((size_t)C, sizeof(double));
    bg->trans = calloc((size_t)C * A, sizeof(double));
    double ctx_total = 0;
    for (int c = 0; c < C; c++) ctx_total += ctx_count[c] + 1.0;
    for (int c = 0; c < C; c++) {
        bg->init[c] = (ctx_count[c] + 1.0) / ctx_total;
        double row = 0;
        for (int r = 0; r < A; r++) row += pair_count[(size_t)c * A + r] + 1.0;
        for (int r = 0; r < A; r++) bg->trans[(size_t)c * A + r] = (pair_count[(size_t)c * A + r] + 1.0) / row;
    }
    for (int r = 0; r < A; r++) bg->p0[r] = (c0[r] + 1.0) / (total + A);
    free(ctx_count);
    free(pair_count);
    snprintf(bg->source, sizeof(bg->source), "order-%d background of %s", order, label);
    return true;
}

long background_valid_positions(const FastaFile *f, AlphabetKind alphabet) {
    long n = 0;
    for (int i = 0; i < f->count; i++)
        for (long j = 0; j < f->records[i].length; j++)
            if (residue_index(alphabet, f->records[i].sequence[j]) >= 0) n++;
    return n;
}

/* Upper limit on how many residues we push mass through. Motifs with a
 * bounded maximum length stop at that length (the matcher stops there
 * too). Unbounded motifs stop once the remaining mass is negligible. */
#define MAX_STEPS 1000000
#define MASS_EPS 1e-15

double specificity_hit_probability(const MatcherProgram *mp, const Background *bg) {
    int n = mp->num_states, A = mp->alphabet_size;
    int C = bg->contexts;                     /* mass is tracked per (DFA state, last k residues) */

    double *mass = calloc((size_t)n * C, sizeof(double));
    double *next = calloc((size_t)n * C, sizeof(double));

    /* Start: DFA at its start state; the k residues before the start
     * position follow the background's own context frequencies. */
    for (int c = 0; c < C; c++) mass[(size_t)mp->start_state * C + c] = bg->init[c];

    double hit = 0.0;
    long limit = (mp->max_len == LEN_INF) ? MAX_STEPS : mp->max_len;

    for (long step = 0; step < limit; step++) {
        memset(next, 0, sizeof(double) * (size_t)n * C);
        double alive = 0.0;
        for (int s = 0; s < n; s++) {
            for (int c = 0; c < C; c++) {
                double m = mass[(size_t)s * C + c];
                if (m == 0.0) continue;
                const double *pr = &bg->trans[(size_t)c * A];
                for (int r = 0; r < A; r++) {
                    int t = mp->transitions[s * A + r];
                    if (t == DFA_DEAD) continue;
                    double flow = m * pr[r];
                    if (mp->accepting[t]) {
                        hit += flow;              /* first acceptance: this start position is a hit */
                    } else {
                        int nc = (C > 1) ? (c * A + r) % C : 0;
                        next[(size_t)t * C + nc] += flow;
                        alive += flow;
                    }
                }
            }
        }
        double *tmp = mass; mass = next; next = tmp;
        if (alive < MASS_EPS) break;
    }

    free(mass);
    free(next);
    return hit;
}

double specificity_empirical_rate(const MatcherProgram *mp, const FastaFile *reference) {
    long positions = background_valid_positions(reference, mp->alphabet);
    if (positions == 0) return 0.0;
    long hits = 0;
    for (int i = 0; i < reference->count; i++) {
        HitList list;
        matcher_scan(mp, reference->records[i].sequence, reference->records[i].length, &list);
        hits += list.count;
        hitlist_free(&list);
    }
    return (double)hits / (double)positions;
}

double rate_of(const MatcherProgram *mp, const RateSource *src) {
    if (src->reference) return specificity_empirical_rate(mp, src->reference);
    return specificity_hit_probability(mp, src->background);
}

const char *rate_source_text(const RateSource *src, char *buf, size_t n) {
    if (src->reference) snprintf(buf, n, "empirical rate on %s", src->reference_label);
    else snprintf(buf, n, "%s", src->background->source);
    return buf;
}
