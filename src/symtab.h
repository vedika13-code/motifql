#ifndef MOTIFQL_SYMTAB_H
#define MOTIFQL_SYMTAB_H

#include "common.h"
#include "alphabet.h"
#include "ast.h"

/*
 * The symbol table for MotifQL holds motif definitions. Unlike a
 * general-purpose language there are no nested lexical scopes: motifs
 * live in one flat namespace, because a motif is a top-level named
 * pattern and nothing in the grammar introduces a new binding region.
 *
 * That does NOT make the table trivial. It still has to support
 * forward references (a motif may reference one defined later),
 * detect duplicate definitions, resolve @name references to indices,
 * track which motifs are actually used, and -- the interesting one --
 * detect reference CYCLES, since @a referring to @b referring back to
 * @a describes an infinitely large pattern that cannot be compiled.
 */

/*
 * What one unit of the error budget buys.
 *
 *   ERR_SUBSTITUTION  -- `allow k mismatch`: Hamming distance. A match
 *                        is always exactly as long as the motif.
 *   ERR_EDIT          -- `allow k edits`: Levenshtein distance. An edit
 *                        is a substitution, an inserted residue, or a
 *                        deleted pattern position, so match length can
 *                        shrink or grow by up to k.
 *
 * In both models every residue an error introduces is drawn from the
 * declared alphabet, never "any character" -- the alphabet typing that
 * the rest of the compiler enforces carries through to approximate
 * matching.
 */
typedef enum {
    ERR_SUBSTITUTION,
    ERR_EDIT,
} ErrorModel;

typedef struct {
    const char *name;      /* points into the source buffer */
    int name_length;
    int line;
    int column;

    Node *pattern;         /* root of this motif's AST */
    int mismatch_budget;   /* k from `allow k mismatch` or `allow k edits`, 0 if absent */
    ErrorModel error_model;/* which kind of error the budget counts (see below) */

    /* filled by semantic analysis */
    int min_len;
    int max_len;           /* LEN_INF if unbounded */
    bool used;             /* referenced by another motif or a search */
    int visit_state;       /* cycle-detection colouring: 0 unvisited, 1 in-progress, 2 done */
} MotifDef;

typedef struct {
    const char *motif_name;
    int motif_name_length;
    int motif_index;       /* resolved during semantic analysis, -1 if unresolved */
    const char *path;      /* FASTA path from `search X in "path"` */
    int path_length;
    int line;
    int column;
} SearchStmt;

typedef struct {
    bool alphabet_declared;
    AlphabetKind alphabet;
    int alphabet_line;

    MotifDef motifs[MAX_MOTIFS];
    int motif_count;

    SearchStmt searches[MAX_SEARCHES];
    int search_count;
} Program;

void program_init(Program *program);
void program_free(Program *program);

/* Returns the index of an existing motif with this name, or -1. */
int program_find_motif(Program *program, const char *name, int name_length);

/* Adds a motif and returns its index, or -1 if the table is full. */
int program_add_motif(Program *program, const char *name, int name_length, int line, int column);

void symtab_print(Program *program);

#endif
