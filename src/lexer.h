#ifndef MOTIFQL_LEXER_H
#define MOTIFQL_LEXER_H

#include "common.h"

typedef enum {
    /* keywords */
    TK_ALPHABET, TK_MOTIF, TK_SEARCH, TK_IN, TK_ALLOW, TK_MISMATCH,
    /* literals and names */
    TK_IDENT, TK_NUMBER, TK_STRING,
    /* structure */
    TK_LBRACE, TK_RBRACE, TK_LBRACKET, TK_RBRACKET, TK_LPAREN, TK_RPAREN,
    /* pattern operators */
    TK_PIPE, TK_STAR, TK_PLUS, TK_QUESTION, TK_CARET, TK_AT,
    /* punctuation */
    TK_COMMA, TK_SEMICOLON, TK_EQUAL,
    /* control */
    TK_EOF, TK_ERROR,
} TokenType;

typedef struct {
    TokenType type;
    const char *start;  /* points into the source buffer; not NUL-terminated */
    int length;
    int line;
    int column;
    long number;        /* valid when type == TK_NUMBER */
} Token;

typedef struct {
    const char *source;
    const char *start;
    const char *current;
    int line;
    int line_start_offset;
} Lexer;

void lexer_init(Lexer *lexer, const char *source);
Token lexer_next_token(Lexer *lexer);

/* Human-readable token name, used by the --tokens dump mode so the
 * lexical analysis phase produces inspectable evidence of its own. */
const char *token_type_name(TokenType type);

#endif
