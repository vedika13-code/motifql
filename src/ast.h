#ifndef MOTIFQL_AST_H
#define MOTIFQL_AST_H

#include "common.h"
#include "alphabet.h"
#include "lexer.h"

/*
 * The AST is the output of syntax analysis and the input to semantic
 * analysis. Note what it deliberately does NOT contain: resolved
 * symbol masks, resolved motif references, or inferred match lengths.
 * The parser records only what the source text SAYS; deciding what it
 * MEANS (is 'U' a legal residue here? does @tata refer to anything?)
 * is semantic analysis's job, and those fields are filled in during
 * that pass. Keeping the split strict is what makes it possible to
 * report a clean syntax tree even for a program that fails semantics.
 */

typedef enum {
    NODE_LITERAL,    /* "ATGN"  -- a run of residue symbols        */
    NODE_CLASS,      /* [ACG] or [^T]                              */
    NODE_CONCAT,     /* juxtaposition: A B                          */
    NODE_ALTERNATE,  /* A | B                                       */
    NODE_REPEAT,     /* A*, A+, A?, A{n}, A{n,m}, A{n,}             */
    NODE_REF,        /* @othermotif                                 */
} NodeType;

typedef struct Node Node;

struct Node {
    NodeType type;
    int line;
    int column;

    /* --- filled by the parser (syntax) --- */
    const char *text;      /* NODE_LITERAL/NODE_CLASS/NODE_REF: raw source span */
    int text_length;
    bool negated;          /* NODE_CLASS: [^...] */
    int repeat_min;
    int repeat_max;        /* REPEAT_INF for unbounded */
    Node **children;
    int child_count;
    int child_capacity;

    /* --- filled by semantic analysis --- */
    SymMask *masks;        /* NODE_LITERAL: one mask per residue position */
    int mask_count;
    SymMask class_mask;    /* NODE_CLASS: resolved (and negation applied) */
    int ref_index;         /* NODE_REF: index into the motif table, -1 if unresolved */
    int min_len;           /* inferred minimum match length */
    int max_len;           /* inferred maximum match length, LEN_INF if unbounded */
};

Node *node_new(NodeType type, int line, int column);
void node_add_child(Node *parent, Node *child);
void node_free(Node *node);

const char *node_type_name(NodeType type);

/* Pretty-prints the tree with indentation. `alphabet_known` controls
 * whether resolved masks are shown (they only exist after semantic
 * analysis has run successfully). */
void ast_print(Node *node, int indent, AlphabetKind kind, bool resolved);

#endif
