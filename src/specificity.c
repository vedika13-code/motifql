#include "specificity.h"

static int residue_index(AlphabetKind kind, char c) {
    const char *residues = alphabet_residues(kind);
    const char *p = strchr(residues, c);
    return (p && c) ? (int)(p - residues) : -1;
}

void background_uniform(Background *bg, AlphabetKind alphabet) {
    memset(bg, 0, sizeof(*bg));
    bg->alphabet = alphabet;
    bg->size = alphabet_size(alphabet);
    bg->order = 0;
    for (int r = 0; r < bg->size; r++) bg->p0[r] = 1.0 / bg->size;
    snprintf(bg->source, sizeof(bg->source), "uniform %s", alphabet_name(alphabet));
}

bool background_from_fasta(Background *bg, AlphabetKind alphabet, const FastaFile *f, int order,
                           const char *label) {
    memset(bg, 0, sizeof(*bg));
    bg->alphabet = alphabet;
    bg->size = alphabet_size(alphabet);
    bg->order = order;

    double c0[BG_MAX] = {0}, c1[BG_MAX][BG_MAX] = {{0}};
    double total = 0;
    for (int i = 0; i < f->count; i++) {
        int prev = -1;
        for (long j = 0; j < f->records[i].length; j++) {
            int r = residue_index(alphabet, f->records[i].sequence[j]);
            if (r < 0) { prev = -1; continue; }
            c0[r] += 1; total += 1;
            if (prev >= 0) c1[prev][r] += 1;
            prev = r;
        }
    }
    if (total == 0) return false;

    /* A pseudocount of 1 keeps every transition possible, so a residue
     * pair that happens to be absent from a short file is treated as
     * rare rather than impossible. */
    for (int r = 0; r < bg->size; r++) bg->p0[r] = (c0[r] + 1.0) / (total + bg->size);
    for (int a = 0; a < bg->size; a++) {
        double row = 0;
        for (int b = 0; b < bg->size; b++) row += c1[a][b] + 1.0;
        for (int b = 0; b < bg->size; b++) bg->p1[a][b] = (c1[a][b] + 1.0) / row;
    }
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
    int ctx = (bg->order == 1) ? A : 1;           /* order 1 tracks the previous residue */

    double *mass = calloc((size_t)n * ctx, sizeof(double));
    double *next = calloc((size_t)n * ctx, sizeof(double));

    /* Start: DFA at its start state. For order 1 the residue before
     * the start position follows the background's own composition. */
    for (int c = 0; c < ctx; c++) mass[mp->start_state * ctx + c] = (ctx == 1) ? 1.0 : bg->p0[c];

    double hit = 0.0;
    long limit = (mp->max_len == LEN_INF) ? MAX_STEPS : mp->max_len;

    for (long step = 0; step < limit; step++) {
        memset(next, 0, sizeof(double) * (size_t)n * ctx);
        double alive = 0.0;
        for (int s = 0; s < n; s++) {
            for (int c = 0; c < ctx; c++) {
                double m = mass[s * ctx + c];
                if (m == 0.0) continue;
                for (int r = 0; r < A; r++) {
                    double pr = (ctx == 1) ? bg->p0[r] : bg->p1[c][r];
                    int t = mp->transitions[s * A + r];
                    if (t == DFA_DEAD) continue;
                    double flow = m * pr;
                    if (mp->accepting[t]) {
                        hit += flow;              /* first acceptance: this start position is a hit */
                    } else {
                        next[t * ctx + (ctx == 1 ? 0 : r)] += flow;
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
