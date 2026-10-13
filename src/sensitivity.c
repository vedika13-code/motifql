#include "sensitivity.h"
#include "matcher.h"
#include <math.h>

/* xorshift64*: small, fast, and identical on every platform, so a
 * given seed gives the same sensitivity figure on any machine. */
static uint64_t rng_next(uint64_t *s) {
    uint64_t x = *s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    *s = x;
    return x * 2685821657736338717ULL;
}
static double rng_unit(uint64_t *s) { return (double)(rng_next(s) >> 11) / 9007199254740992.0; }
static int rng_below(uint64_t *s, int n) { return (int)(rng_unit(s) * n); }

static char pick_from_mask(AlphabetKind kind, SymMask mask, uint64_t *rng) {
    const char *residues = alphabet_residues(kind);
    int choices[32], count = 0;
    for (int r = 0; r < alphabet_size(kind); r++) if (mask & ((SymMask)1u << r)) choices[count++] = r;
    if (count == 0) return residues[0];
    return residues[choices[rng_below(rng, count)]];
}

typedef struct { char *buf; int len; int cap; bool overflow; } Out;

static void emit(Out *o, char c) {
    if (o->len + 1 >= o->cap) { o->overflow = true; return; }
    o->buf[o->len++] = c;
}

/* Unbounded repeats ({n,}, *, +) draw extra copies geometrically with
 * mean 1, capped, so instances stay short and finite. */
#define UNBOUNDED_EXTRA_CAP 20

static void sample_node(const Program *program, const Node *node, Out *o, uint64_t *rng) {
    if (!node || o->overflow) return;
    AlphabetKind kind = program->alphabet;
    switch (node->type) {
        case NODE_LITERAL:
            for (int i = 0; i < node->mask_count; i++) emit(o, pick_from_mask(kind, node->masks[i], rng));
            break;
        case NODE_CLASS:
            emit(o, pick_from_mask(kind, node->class_mask, rng));
            break;
        case NODE_CONCAT:
            for (int i = 0; i < node->child_count; i++) sample_node(program, node->children[i], o, rng);
            break;
        case NODE_ALTERNATE:
            if (node->child_count > 0) sample_node(program, node->children[rng_below(rng, node->child_count)], o, rng);
            break;
        case NODE_REPEAT: {
            int count = node->repeat_min;
            if (node->repeat_max == REPEAT_INF) {
                int extra = 0;
                while (extra < UNBOUNDED_EXTRA_CAP && rng_unit(rng) < 0.5) extra++;
                count += extra;
            } else if (node->repeat_max > node->repeat_min) {
                count += rng_below(rng, node->repeat_max - node->repeat_min + 1);
            }
            for (int i = 0; i < count; i++) sample_node(program, node->children[0], o, rng);
            break;
        }
        case NODE_REF:
            if (node->ref_index >= 0) sample_node(program, program->motifs[node->ref_index].pattern, o, rng);
            break;
    }
}

int sensitivity_sample_instance(const Program *program, int motif_index, char *buf, int cap, uint64_t *rng) {
    Out o = { buf, 0, cap, false };
    sample_node(program, program->motifs[motif_index].pattern, &o, rng);
    if (o.overflow) return -1;
    buf[o.len] = '\0';
    return o.len;
}

static int mutate(const char *in, int len, char *out, int cap, AlphabetKind kind,
                  const MutationModel *m, uint64_t *rng) {
    const char *residues = alphabet_residues(kind);
    int A = alphabet_size(kind), n = 0;
    for (int i = 0; i < len; i++) {
        double u = rng_unit(rng);
        if (u < m->indel / 2) continue;                       /* deletion */
        if (u < m->indel && n + 1 < cap) out[n++] = residues[rng_below(rng, A)];   /* insertion */
        char c = in[i];
        if (rng_unit(rng) < m->substitution) {
            char d;
            do { d = residues[rng_below(rng, A)]; } while (d == c);
            c = d;
        }
        if (n + 1 < cap) out[n++] = c;
    }
    out[n] = '\0';
    return n;
}

/* Wilson score interval: well behaved even when the estimate is near 0 or 1. */
static void wilson(int found, int n, double *lo, double *hi) {
    if (n == 0) { *lo = 0; *hi = 1; return; }
    double z = 1.96, p = (double)found / n, z2 = z * z;
    double centre = (p + z2 / (2.0 * n)) / (1 + z2 / n);
    double half = z * sqrt(p * (1 - p) / n + z2 / (4.0 * n * n)) / (1 + z2 / n);
    *lo = centre - half < 0 ? 0 : centre - half;
    *hi = centre + half > 1 ? 1 : centre + half;
}

#define SITE_CAP 4096

SensitivityResult sensitivity_estimate(const Program *program, int motif_index, const MatcherProgram *mp,
                                       const MutationModel *model, int samples, uint64_t seed) {
    SensitivityResult r = {0};
    uint64_t rng = seed ? seed : 0x9E3779B97F4A7C15ULL;
    char site[SITE_CAP], mutated[SITE_CAP * 2];
    for (int t = 0; t < samples; t++) {
        int len = sensitivity_sample_instance(program, motif_index, site, SITE_CAP, &rng);
        if (len < 0) continue;
        int mlen = mutate(site, len, mutated, (int)sizeof(mutated), program->alphabet, model, &rng);
        HitList hits;
        matcher_scan(mp, mutated, mlen, &hits);
        if (hits.count > 0) r.found++;
        hitlist_free(&hits);
        r.samples++;
    }
    r.sensitivity = r.samples ? (double)r.found / r.samples : 0;
    wilson(r.found, r.samples, &r.ci_low, &r.ci_high);
    return r;
}

SensitivityResult sensitivity_on_sites(const MatcherProgram *mp, const FastaFile *sites) {
    SensitivityResult r = {0};
    for (int i = 0; i < sites->count; i++) {
        HitList hits;
        matcher_scan(mp, sites->records[i].sequence, sites->records[i].length, &hits);
        if (hits.count > 0) r.found++;
        hitlist_free(&hits);
        r.samples++;
    }
    r.sensitivity = r.samples ? (double)r.found / r.samples : 0;
    wilson(r.found, r.samples, &r.ci_low, &r.ci_high);
    return r;
}
