#ifndef MOTIFQL_DIAGNOSTICS_H
#define MOTIFQL_DIAGNOSTICS_H

#include "common.h"

/*
 * One diagnostics sink shared by all three error-producing phases, so
 * a lexical error, a syntax error, and a semantic error all print in
 * the same format with the same source-line caret. The phase label is
 * carried on the message so the output makes it obvious which stage
 * of the compiler rejected the program.
 */

typedef enum {
    PHASE_LEXICAL,
    PHASE_SYNTAX,
    PHASE_SEMANTIC,
    PHASE_SPECIFICITY,  /* compile-time statistical analysis of the compiled automaton */
} CompilerPhase;

typedef enum {
    SEV_ERROR,
    SEV_WARNING,
} Severity;

typedef struct {
    const char *source;   /* whole source buffer, for echoing the offending line */
    const char *filename;
    int error_count;
    int warning_count;
} Diagnostics;

void diag_init(Diagnostics *d, const char *source, const char *filename);

void diag_report(Diagnostics *d, CompilerPhase phase, Severity severity,
                 int line, int column, const char *message);

/* Convenience wrapper that formats the message first. */
void diag_reportf(Diagnostics *d, CompilerPhase phase, Severity severity,
                  int line, int column, const char *fmt, ...);

const char *phase_name(CompilerPhase phase);

#endif
