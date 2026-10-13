#include "diagnostics.h"
#include <stdarg.h>

void diag_init(Diagnostics *d, const char *source, const char *filename) {
    d->source = source;
    d->filename = filename;
    d->error_count = 0;
    d->warning_count = 0;
}

const char *phase_name(CompilerPhase phase) {
    switch (phase) {
        case PHASE_LEXICAL:  return "lexical";
        case PHASE_SYNTAX:   return "syntax";
        case PHASE_SEMANTIC: return "semantic";
        case PHASE_SPECIFICITY: return "specificity";
    }
    return "?";
}

/* Echoes the source line containing the error, with a caret under the
 * offending column. Costs a scan through the source per diagnostic,
 * which is irrelevant at the file sizes this language deals with and
 * makes errors far easier to act on. */
static void echo_source_line(const char *source, int line, int column) {
    if (!source || line <= 0) return;

    const char *p = source;
    int current = 1;
    while (current < line && *p) {
        if (*p == '\n') current++;
        p++;
    }
    if (current != line) return;

    const char *end = p;
    while (*end && *end != '\n') end++;

    fprintf(stderr, "  %5d | %.*s\n", line, (int)(end - p), p);
    fprintf(stderr, "        | ");
    for (int i = 1; i < column; i++) fputc(' ', stderr);
    fprintf(stderr, "^\n");
}

void diag_report(Diagnostics *d, CompilerPhase phase, Severity severity,
                 int line, int column, const char *message) {
    const char *sev = (severity == SEV_ERROR) ? "error" : "warning";
    if (severity == SEV_ERROR) d->error_count++;
    else d->warning_count++;

    fprintf(stderr, "%s:%d:%d: %s %s: %s\n",
            d->filename ? d->filename : "<input>",
            line, column, phase_name(phase), sev, message);
    echo_source_line(d->source, line, column);
}

void diag_reportf(Diagnostics *d, CompilerPhase phase, Severity severity,
                  int line, int column, const char *fmt, ...) {
    char buffer[512];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buffer, sizeof(buffer), fmt, args);
    va_end(args);
    diag_report(d, phase, severity, line, column, buffer);
}
