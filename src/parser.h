#ifndef MOTIFQL_PARSER_H
#define MOTIFQL_PARSER_H

#include "common.h"
#include "symtab.h"
#include "diagnostics.h"

/* Runs lexical + syntax analysis over `source`, filling `program`
 * with motif definitions and their ASTs. Returns true when no errors
 * were reported. Parsing continues past errors via panic-mode
 * recovery, so a single run reports as many distinct problems as it
 * can rather than stopping at the first. */
bool parse_program(const char *source, Program *program, Diagnostics *diag);

#endif
