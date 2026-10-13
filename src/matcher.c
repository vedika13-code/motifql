#include "matcher.h"

void hitlist_free(HitList *list) {
    free(list->hits);
    memset(list, 0, sizeof(*list));
}

static void hitlist_push(HitList *list, long start, long length) {
    if (list->count + 1 > list->capacity) {
        int new_cap = list->capacity == 0 ? 16 : list->capacity * 2;
        Hit *grown = realloc(list->hits, sizeof(Hit) * (size_t)new_cap);
        if (!grown) { fprintf(stderr, "out of memory recording hits\n"); exit(1); }
        list->hits = grown;
        list->capacity = new_cap;
    }
    list->hits[list->count].start = start;
    list->hits[list->count].length = length;
    list->count++;
}

static int residue_index(AlphabetKind kind, char c) {
    const char *residues = alphabet_residues(kind);
    const char *p = strchr(residues, c);
    return p ? (int)(p - residues) : -1;
}

/* Runs the DFA starting at sequence[from], for as long as it stays
 * alive or the sequence ends. Returns the longest length (>=1) at
 * which an accepting state was reached, or 0 if none was. */
static long run_from(const MatcherProgram *mp, const char *sequence, long len, long from) {
    int state = mp->start_state;
    long best = 0;

    for (long i = from; i < len; i++) {
        int r = residue_index(mp->alphabet, sequence[i]);
        if (r < 0) break; /* a character outside the alphabet ends this attempt */

        int next = mp->transitions[state * mp->alphabet_size + r];
        if (next == DFA_DEAD) break;

        state = next;
        long length_here = i - from + 1;
        if (mp->accepting[state]) best = length_here;

        /* No point continuing once we've exceeded the motif's own
         * inferred maximum length -- Phase 1's length inference
         * paying for itself a second time, this time as a scan bound
         * rather than a diagnostic. */
        if (mp->max_len != LEN_INF && length_here >= mp->max_len) break;
    }
    return best;
}

void matcher_scan(const MatcherProgram *mp, const char *sequence, long len, HitList *hits) {
    memset(hits, 0, sizeof(*hits));

    long prefix_len = (long)strlen(mp->literal_prefix);

    for (long p = 0; p < len; p++) {
        if (prefix_len > 0) {
            if (p + prefix_len > len) break; /* not enough room left for the prefix at all */
            if (strncmp(sequence + p, mp->literal_prefix, (size_t)prefix_len) != 0) continue;
        }

        long matched_len = run_from(mp, sequence, len, p);
        if (matched_len > 0) hitlist_push(hits, p, matched_len);
    }
}
