#ifndef MOTIFQL_MATCHER_H
#define MOTIFQL_MATCHER_H

#include "common.h"
#include "codegen.h"
#include "fasta.h"

typedef struct {
    long start;    /* 0-based offset into the record's sequence */
    long length;
} Hit;

typedef struct {
    Hit *hits;
    int count;
    int capacity;
} HitList;

void hitlist_free(HitList *list);

/* Scans one sequence buffer of length `len` against a compiled
 * matcher program, appending every hit found to `hits`. Overlapping
 * hits are all reported (a match starting at position p does not
 * suppress a match starting at p+1) -- see docs/PHASE2.md for why
 * that is the appropriate default for motif search.
 *
 * Policy for one start position: run the DFA from that position for
 * as long as it stays alive; if it passes through one or more
 * accepting states, report the LONGEST accepting length reached
 * (greedy), then move on to the next start position regardless.
 *
 * When the compiled program has a non-empty literal prefix, start
 * positions where the sequence cannot possibly match that prefix are
 * skipped without running the automaton at all -- the first payoff
 * from the length/prefix analysis Phase 1 and this phase compute. */
void matcher_scan(const MatcherProgram *mp, const char *sequence, long len, HitList *hits);

#endif
