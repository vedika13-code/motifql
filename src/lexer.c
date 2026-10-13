#include "lexer.h"

void lexer_init(Lexer *lx, const char *source) {
    lx->source = source;
    lx->start = source;
    lx->current = source;
    lx->line = 1;
    lx->line_start_offset = 0;
}

const char *token_type_name(TokenType type) {
    switch (type) {
        case TK_ALPHABET:  return "ALPHABET";
        case TK_MOTIF:     return "MOTIF";
        case TK_SEARCH:    return "SEARCH";
        case TK_IN:        return "IN";
        case TK_ALLOW:     return "ALLOW";
        case TK_MISMATCH:  return "MISMATCH";
        case TK_EDITS:     return "EDITS";
        case TK_REQUIRE:   return "REQUIRE";
        case TK_LE:        return "LE";
        case TK_PERCENT:   return "PERCENT";
        case TK_IDENT:     return "IDENT";
        case TK_NUMBER:    return "NUMBER";
        case TK_STRING:    return "STRING";
        case TK_LBRACE:    return "LBRACE";
        case TK_RBRACE:    return "RBRACE";
        case TK_LBRACKET:  return "LBRACKET";
        case TK_RBRACKET:  return "RBRACKET";
        case TK_LPAREN:    return "LPAREN";
        case TK_RPAREN:    return "RPAREN";
        case TK_PIPE:      return "PIPE";
        case TK_STAR:      return "STAR";
        case TK_PLUS:      return "PLUS";
        case TK_QUESTION:  return "QUESTION";
        case TK_CARET:     return "CARET";
        case TK_AT:        return "AT";
        case TK_COMMA:     return "COMMA";
        case TK_SEMICOLON: return "SEMICOLON";
        case TK_EQUAL:     return "EQUAL";
        case TK_EOF:       return "EOF";
        case TK_ERROR:     return "ERROR";
    }
    return "?";
}

static bool is_at_end(Lexer *lx) { return *lx->current == '\0'; }
static char peek(Lexer *lx) { return *lx->current; }
static char peek_next(Lexer *lx) { return is_at_end(lx) ? '\0' : lx->current[1]; }
static char advance_char(Lexer *lx) { return *lx->current++; }

static bool is_digit(char c) { return c >= '0' && c <= '9'; }
static bool is_alpha(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_';
}

static int current_column(Lexer *lx, const char *at) {
    return (int)(at - lx->source) - lx->line_start_offset + 1;
}

static Token make_token(Lexer *lx, TokenType type) {
    Token t;
    t.type = type;
    t.start = lx->start;
    t.length = (int)(lx->current - lx->start);
    t.line = lx->line;
    t.column = current_column(lx, lx->start);
    t.number = 0;
    return t;
}

static Token error_token(Lexer *lx, const char *message) {
    Token t;
    t.type = TK_ERROR;
    t.start = message;   /* for TK_ERROR, start/length carry the message */
    t.length = (int)strlen(message);
    t.line = lx->line;
    t.column = current_column(lx, lx->start);
    t.number = 0;
    return t;
}

static void skip_whitespace_and_comments(Lexer *lx) {
    for (;;) {
        char c = peek(lx);
        if (c == ' ' || c == '\t' || c == '\r') {
            advance_char(lx);
        } else if (c == '\n') {
            advance_char(lx);
            lx->line++;
            lx->line_start_offset = (int)(lx->current - lx->source);
        } else if (c == '/' && peek_next(lx) == '/') {
            while (peek(lx) != '\n' && !is_at_end(lx)) advance_char(lx);
        } else if (c == '#') {
            /* '#' comments too: FASTA-adjacent tooling conventionally
             * uses them, and users coming from bed/gff files expect it. */
            while (peek(lx) != '\n' && !is_at_end(lx)) advance_char(lx);
        } else {
            return;
        }
    }
}

static TokenType keyword_type(const char *start, int length) {
    struct { const char *word; int len; TokenType type; } keywords[] = {
        {"alphabet", 8, TK_ALPHABET},
        {"motif",    5, TK_MOTIF},
        {"search",   6, TK_SEARCH},
        {"in",       2, TK_IN},
        {"allow",    5, TK_ALLOW},
        {"mismatch", 8, TK_MISMATCH},
        {"edits",    5, TK_EDITS},
        {"edit",     4, TK_EDITS},
        {"require",  7, TK_REQUIRE},
        {NULL, 0, TK_IDENT},
    };
    for (int i = 0; keywords[i].word; i++) {
        if (length == keywords[i].len && memcmp(start, keywords[i].word, (size_t)length) == 0) {
            return keywords[i].type;
        }
    }
    return TK_IDENT;
}

static Token identifier(Lexer *lx) {
    while (is_alpha(peek(lx)) || is_digit(peek(lx))) advance_char(lx);
    int length = (int)(lx->current - lx->start);
    return make_token(lx, keyword_type(lx->start, length));
}

static Token number(Lexer *lx) {
    while (is_digit(peek(lx))) advance_char(lx);
    Token t = make_token(lx, TK_NUMBER);
    char buf[32];
    int len = t.length < 31 ? t.length : 31;
    memcpy(buf, t.start, (size_t)len);
    buf[len] = '\0';
    t.number = strtol(buf, NULL, 10);
    return t;
}

/* Sequence literals are quoted: "ATG", "TATAWAWR".
 *
 * Requiring quotes is a deliberate design decision (see LANGUAGE.md):
 * without them, a bare motif reference like `N` would be
 * indistinguishable from an identifier at the token level, forcing
 * the lexer to know whether it sits inside a pattern. Quotes keep the
 * lexer context-free. */
static Token string_literal(Lexer *lx) {
    while (peek(lx) != '"' && !is_at_end(lx)) {
        if (peek(lx) == '\n') {
            return error_token(lx, "unterminated sequence literal (newline inside quotes)");
        }
        advance_char(lx);
    }
    if (is_at_end(lx)) {
        return error_token(lx, "unterminated sequence literal (reached end of file)");
    }
    advance_char(lx); /* closing quote */

    /* Trim the surrounding quotes from the reported span so the parser
     * sees only the residue characters. */
    Token t = make_token(lx, TK_STRING);
    t.start += 1;
    t.length -= 2;
    return t;
}

Token lexer_next_token(Lexer *lx) {
    skip_whitespace_and_comments(lx);
    lx->start = lx->current;

    if (is_at_end(lx)) return make_token(lx, TK_EOF);

    char c = advance_char(lx);

    if (is_alpha(c)) return identifier(lx);
    if (is_digit(c)) return number(lx);

    switch (c) {
        case '"': return string_literal(lx);
        case '{': return make_token(lx, TK_LBRACE);
        case '}': return make_token(lx, TK_RBRACE);
        case '[': return make_token(lx, TK_LBRACKET);
        case ']': return make_token(lx, TK_RBRACKET);
        case '(': return make_token(lx, TK_LPAREN);
        case ')': return make_token(lx, TK_RPAREN);
        case '|': return make_token(lx, TK_PIPE);
        case '*': return make_token(lx, TK_STAR);
        case '+': return make_token(lx, TK_PLUS);
        case '?': return make_token(lx, TK_QUESTION);
        case '^': return make_token(lx, TK_CARET);
        case '@': return make_token(lx, TK_AT);
        case ',': return make_token(lx, TK_COMMA);
        case ';': return make_token(lx, TK_SEMICOLON);
        case '=': return make_token(lx, TK_EQUAL);
        case '%': return make_token(lx, TK_PERCENT);
        case '<':
            if (peek(lx) == '=') { advance_char(lx); return make_token(lx, TK_LE); }
            return error_token(lx, "expect '<=' (MotifQL has no plain '<')");
        default:  return error_token(lx, "unexpected character");
    }
}
