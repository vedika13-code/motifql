#include "parser.h"
#include "lexer.h"
#include "ast.h"

typedef struct {
    Lexer lexer;
    Token current;
    Token previous;
    Program *program;
    Diagnostics *diag;
    bool had_error;
    bool panic_mode;
} Parser;

/* ---- error handling ------------------------------------------------------
 * Panic mode: after the first error inside a construct, further errors
 * are suppressed until synchronize() reaches a token that can safely
 * start a new declaration. Without this, one stray token produces a
 * cascade of meaningless follow-on errors. */

static void error_at(Parser *p, Token *token, const char *message) {
    if (p->panic_mode) return;
    p->panic_mode = true;
    p->had_error = true;
    diag_report(p->diag, PHASE_SYNTAX, SEV_ERROR, token->line, token->column, message);
}

static void error_at_current(Parser *p, const char *message) {
    error_at(p, &p->current, message);
}

static void advance(Parser *p) {
    p->previous = p->current;
    for (;;) {
        p->current = lexer_next_token(&p->lexer);
        if (p->current.type != TK_ERROR) break;
        /* TK_ERROR carries its message in start/length; this is the
         * lexical phase reporting through the same sink. */
        if (!p->panic_mode) {
            p->panic_mode = true;
            p->had_error = true;
            diag_reportf(p->diag, PHASE_LEXICAL, SEV_ERROR,
                         p->current.line, p->current.column,
                         "%.*s", p->current.length, p->current.start);
        }
    }
}

static bool check(Parser *p, TokenType type) { return p->current.type == type; }

static bool match(Parser *p, TokenType type) {
    if (!check(p, type)) return false;
    advance(p);
    return true;
}

static void consume(Parser *p, TokenType type, const char *message) {
    if (check(p, type)) { advance(p); return; }
    error_at_current(p, message);
}

static void synchronize(Parser *p) {
    p->panic_mode = false;
    while (!check(p, TK_EOF)) {
        if (p->previous.type == TK_SEMICOLON) return;
        switch (p->current.type) {
            case TK_ALPHABET:
            case TK_MOTIF:
            case TK_SEARCH:
                return;
            default:
                break;
        }
        advance(p);
    }
}

/* ---- pattern grammar -----------------------------------------------------
 * alternation := concat ( "|" concat )*
 * concat      := repeat+
 * repeat      := atom quantifier?
 * quantifier  := "*" | "+" | "?" | "{" NUMBER ( "," NUMBER? )? "}"
 * atom        := STRING | class | "(" alternation ")" | "@" IDENT
 * class       := "[" "^"? IDENT "]"
 */

static Node *parse_alternation(Parser *p);

static Node *parse_class(Parser *p) {
    Token bracket = p->previous;
    Node *node = node_new(NODE_CLASS, bracket.line, bracket.column);
    node->negated = match(p, TK_CARET);

    if (!check(p, TK_IDENT)) {
        error_at_current(p, "expect residue characters inside '[...]'");
        return node;
    }
    advance(p);
    node->text = p->previous.start;
    node->text_length = p->previous.length;

    consume(p, TK_RBRACKET, "expect ']' to close character class");
    return node;
}

static Node *parse_atom(Parser *p) {
    if (match(p, TK_STRING)) {
        Node *node = node_new(NODE_LITERAL, p->previous.line, p->previous.column);
        node->text = p->previous.start;
        node->text_length = p->previous.length;
        if (node->text_length == 0) {
            error_at(p, &p->previous, "empty sequence literal \"\" is not a valid pattern");
        }
        return node;
    }
    if (match(p, TK_LBRACKET)) {
        return parse_class(p);
    }
    if (match(p, TK_LPAREN)) {
        Node *inner = parse_alternation(p);
        consume(p, TK_RPAREN, "expect ')' after group");
        return inner;
    }
    if (match(p, TK_AT)) {
        Token at = p->previous;
        if (!check(p, TK_IDENT)) {
            error_at_current(p, "expect motif name after '@'");
            return node_new(NODE_REF, at.line, at.column);
        }
        advance(p);
        Node *node = node_new(NODE_REF, at.line, at.column);
        node->text = p->previous.start;
        node->text_length = p->previous.length;
        return node;
    }

    error_at_current(p, "expect a pattern element (sequence literal, '[', '(' or '@')");
    return NULL;
}

/* Parses a trailing quantifier, if present, and wraps `atom` in a
 * NODE_REPEAT. Bound *validity* (min <= max) is checked in semantic
 * analysis, not here: syntactically `{5,2}` is a well-formed
 * quantifier, and reporting it as a syntax error would misattribute
 * the problem. */
static Node *parse_quantifier(Parser *p, Node *atom) {
    if (!atom) return NULL;

    int min, max;
    Token op;

    if (match(p, TK_STAR)) {
        op = p->previous; min = 0; max = REPEAT_INF;
    } else if (match(p, TK_PLUS)) {
        op = p->previous; min = 1; max = REPEAT_INF;
    } else if (match(p, TK_QUESTION)) {
        op = p->previous; min = 0; max = 1;
    } else if (match(p, TK_LBRACE)) {
        op = p->previous;
        if (!check(p, TK_NUMBER)) {
            error_at_current(p, "expect a number after '{' in quantifier");
            return atom;
        }
        advance(p);
        min = (int)p->previous.number;
        max = min;

        if (match(p, TK_COMMA)) {
            if (check(p, TK_NUMBER)) {
                advance(p);
                max = (int)p->previous.number;
            } else {
                max = REPEAT_INF;   /* {n,} */
            }
        }
        consume(p, TK_RBRACE, "expect '}' to close quantifier");
    } else {
        return atom;
    }

    Node *repeat = node_new(NODE_REPEAT, op.line, op.column);
    repeat->repeat_min = min;
    repeat->repeat_max = max;
    node_add_child(repeat, atom);
    return repeat;
}

static bool starts_atom(Parser *p) {
    switch (p->current.type) {
        case TK_STRING: case TK_LBRACKET: case TK_LPAREN: case TK_AT:
            return true;
        default:
            return false;
    }
}

static Node *parse_concat(Parser *p) {
    Token start = p->current;
    Node *first = parse_quantifier(p, parse_atom(p));
    if (!first) return NULL;

    if (!starts_atom(p)) return first;

    Node *concat = node_new(NODE_CONCAT, start.line, start.column);
    node_add_child(concat, first);
    while (starts_atom(p)) {
        Node *next = parse_quantifier(p, parse_atom(p));
        if (!next) break;
        node_add_child(concat, next);
    }
    return concat;
}

static Node *parse_alternation(Parser *p) {
    Token start = p->current;
    Node *first = parse_concat(p);
    if (!first) return NULL;
    if (!check(p, TK_PIPE)) return first;

    Node *alt = node_new(NODE_ALTERNATE, start.line, start.column);
    node_add_child(alt, first);
    while (match(p, TK_PIPE)) {
        Node *next = parse_concat(p);
        if (!next) break;
        node_add_child(alt, next);
    }
    return alt;
}

/* ---- declarations --------------------------------------------------------- */

static void parse_alphabet_decl(Parser *p) {
    Token kw = p->previous;
    if (!check(p, TK_IDENT)) {
        error_at_current(p, "expect an alphabet name (DNA, RNA or PROTEIN)");
        return;
    }
    advance(p);
    Token name = p->previous;

    AlphabetKind kind;
    if (!alphabet_from_name(name.start, name.length, &kind)) {
        diag_reportf(p->diag, PHASE_SEMANTIC, SEV_ERROR, name.line, name.column,
                     "unknown alphabet '%.*s' (expected DNA, RNA or PROTEIN)",
                     name.length, name.start);
        p->had_error = true;
    } else if (p->program->alphabet_declared) {
        diag_reportf(p->diag, PHASE_SEMANTIC, SEV_ERROR, name.line, name.column,
                     "alphabet already declared on line %d", p->program->alphabet_line);
        p->had_error = true;
    } else {
        p->program->alphabet = kind;
        p->program->alphabet_declared = true;
        p->program->alphabet_line = kw.line;
    }

    consume(p, TK_SEMICOLON, "expect ';' after alphabet declaration");
}

static bool check_word(Parser *p, const char *word) {
    int n = (int)strlen(word);
    return p->current.type == TK_IDENT && p->current.length == n && memcmp(p->current.start, word, (size_t)n) == 0;
}

static void expect_word(Parser *p, const char *word, const char *message) {
    if (check_word(p, word)) { advance(p); return; }
    error_at_current(p, message);
}

static void parse_mismatch_modifier(Parser *p, MotifDef *motif) {
    if (!match(p, TK_ALLOW)) return;
    if (match(p, TK_QUESTION)) {
        /* `allow ? edits`: the budget is left for the compiler to choose
         * (budget synthesis, driven by a `require` clause). */
        motif->mismatch_budget = BUDGET_FREE;
    } else if (!check(p, TK_NUMBER)) {
        error_at_current(p, "expect a number or '?' after 'allow'");
        return;
    } else {
        advance(p);
        motif->mismatch_budget = (int)p->previous.number;
    }
    if (match(p, TK_MISMATCH)) {
        motif->error_model = ERR_SUBSTITUTION;
    } else if (match(p, TK_EDITS)) {
        motif->error_model = ERR_EDIT;
    } else {
        error_at_current(p, "expect 'mismatch' or 'edits' after the allowance count");
        motif->mismatch_budget = 0;   /* don't carry a half-parsed budget into later dumps */
        return;
    }

    /* Optional position restriction:  at 1,2,5   or   at ?
     * Positions are 1-based residue positions of the motif as written.
     * `at ?` lets budget synthesis choose them. */
    if (check_word(p, "at")) {
        advance(p);
        if (match(p, TK_QUESTION)) {
            motif->positions_free = true;
            return;
        }
        do {
            if (!check(p, TK_NUMBER)) { error_at_current(p, "expect a position number after 'at'"); return; }
            advance(p);
            long pos = p->previous.number;
            if (pos < 1 || pos > 32) { error_at(p, &p->previous, "positions run from 1 to 32"); return; }
            motif->error_positions |= (uint32_t)1u << (pos - 1);
        } while (match(p, TK_COMMA));
    }
}

/* require_clause := "require" "rate" "<=" NUMBER "per" [NUMBER] UNIT
 *                   [ "," "divergence" NUMBER "%" ]
 * UNIT           := "bp" | "kb" | "Mb"
 *
 * `rate`, `per`, `divergence` and the units are contextual words, not
 * reserved keywords, so they stay usable as motif names elsewhere. */
static void parse_require_clause(Parser *p, MotifDef *motif) {
    if (!match(p, TK_REQUIRE)) return;
    Token at = p->previous;
    motif->has_requirement = true;
    motif->divergence = -1.0;
    motif->require_line = at.line;
    motif->require_column = at.column;

    expect_word(p, "rate", "expect 'rate' after 'require'");
    consume(p, TK_LE, "expect '<=' after 'require rate'");
    if (!check(p, TK_NUMBER)) { error_at_current(p, "expect a number of hits after '<='"); return; }
    advance(p);
    double hits = (double)p->previous.number;
    expect_word(p, "per", "expect 'per' after the number of hits");

    double scale = 1.0;
    char scale_txt[16] = "";
    if (check(p, TK_NUMBER)) {
        advance(p);
        scale = (double)p->previous.number;
        snprintf(scale_txt, sizeof(scale_txt), "%ld", p->previous.number);
    }
    double unit = 0;
    if (check_word(p, "bp")) unit = 1;
    else if (check_word(p, "kb")) unit = 1e3;
    else if (check_word(p, "Mb")) unit = 1e6;
    if (unit == 0) { error_at_current(p, "expect a unit: 'bp', 'kb' or 'Mb'"); return; }
    snprintf(motif->required_unit, sizeof(motif->required_unit), "%.15s%.*s", scale_txt,
             p->current.length < 8 ? p->current.length : 8, p->current.start);
    advance(p);

    if (scale <= 0) { error_at(p, &at, "the length in 'per N unit' must be positive"); return; }
    motif->required_rate_shown = hits;
    motif->required_rate = hits / (scale * unit);

    if (match(p, TK_COMMA)) {
        expect_word(p, "divergence", "expect 'divergence' after ','");
        if (!check(p, TK_NUMBER)) { error_at_current(p, "expect a percentage after 'divergence'"); return; }
        advance(p);
        motif->divergence = (double)p->previous.number / 100.0;
        consume(p, TK_PERCENT, "expect '%' after the divergence");
    }
}

static void parse_motif_decl(Parser *p) {
    if (!check(p, TK_IDENT)) {
        error_at_current(p, "expect a motif name after 'motif'");
        return;
    }
    advance(p);
    Token name = p->previous;

    /* Duplicate detection happens here rather than in semantic
     * analysis so the second definition's AST is never built and can
     * never be silently used. */
    int existing = program_find_motif(p->program, name.start, name.length);
    if (existing >= 0) {
        diag_reportf(p->diag, PHASE_SEMANTIC, SEV_ERROR, name.line, name.column,
                     "motif '%.*s' is already defined on line %d",
                     name.length, name.start, p->program->motifs[existing].line);
        p->had_error = true;
    }

    int index = program_add_motif(p->program, name.start, name.length, name.line, name.column);
    if (index < 0) {
        error_at(p, &name, "too many motifs in one program");
        return;
    }
    MotifDef *motif = &p->program->motifs[index];

    Node *pattern = NULL;
    if (match(p, TK_LBRACE)) {
        pattern = parse_alternation(p);
        consume(p, TK_RBRACE, "expect '}' to close the motif body");
        parse_mismatch_modifier(p, motif);
        parse_require_clause(p, motif);
        match(p, TK_SEMICOLON);           /* optional after a braced body */
    } else if (match(p, TK_EQUAL)) {
        pattern = parse_alternation(p);
        parse_mismatch_modifier(p, motif);
        parse_require_clause(p, motif);
        consume(p, TK_SEMICOLON, "expect ';' after a '=' motif definition");
    } else {
        error_at_current(p, "expect '{' or '=' after the motif name");
        return;
    }

    motif->pattern = pattern;
}

static void parse_search_stmt(Parser *p) {
    Token kw = p->previous;
    if (!check(p, TK_IDENT)) {
        error_at_current(p, "expect a motif name after 'search'");
        return;
    }
    advance(p);
    Token name = p->previous;

    consume(p, TK_IN, "expect 'in' after the motif name");

    if (!check(p, TK_STRING)) {
        error_at_current(p, "expect a quoted file path after 'in'");
        return;
    }
    advance(p);
    Token path = p->previous;

    consume(p, TK_SEMICOLON, "expect ';' after a search statement");

    if (p->program->search_count >= MAX_SEARCHES) {
        error_at(p, &kw, "too many search statements in one program");
        return;
    }
    SearchStmt *s = &p->program->searches[p->program->search_count++];
    s->motif_name = name.start;
    s->motif_name_length = name.length;
    s->motif_index = -1;   /* resolved in semantic analysis */
    s->path = path.start;
    s->path_length = path.length;
    s->line = name.line;
    s->column = name.column;
}

static void parse_declaration(Parser *p) {
    if (match(p, TK_ALPHABET)) {
        parse_alphabet_decl(p);
    } else if (match(p, TK_MOTIF)) {
        parse_motif_decl(p);
    } else if (match(p, TK_SEARCH)) {
        parse_search_stmt(p);
    } else {
        error_at_current(p, "expect a declaration ('alphabet', 'motif' or 'search')");
        advance(p);
    }

    if (p->panic_mode) synchronize(p);
}

bool parse_program(const char *source, Program *program, Diagnostics *diag) {
    Parser p;
    memset(&p, 0, sizeof(p));
    lexer_init(&p.lexer, source);
    p.program = program;
    p.diag = diag;

    advance(&p);   /* prime the lookahead token */

    while (!check(&p, TK_EOF)) {
        parse_declaration(&p);
    }

    if (!p.had_error && program->motif_count == 0) {
        diag_report(diag, PHASE_SYNTAX, SEV_ERROR, 1, 1,
                    "program defines no motifs");
        p.had_error = true;
    }

    return !p.had_error;
}
