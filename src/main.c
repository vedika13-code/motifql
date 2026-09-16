#include "common.h"
#include "lexer.h"
#include "parser.h"
#include "sema.h"
#include "symtab.h"
#include "diagnostics.h"

/*
 * Phase 1 driver: front end only.
 *
 * Runs lexical analysis, syntax analysis, and semantic analysis over a
 * .motif source file, and dumps the intermediate products of each
 * stage so every phase produces inspectable evidence of its own work:
 *
 *   --tokens   token stream from the lexer
 *   --ast      annotated syntax tree from the parser + semantic analyser
 *   --symbols  symbol table with resolved references and inferred lengths
 *
 * With no flags, all three are printed. Phase 2 adds NFA/DFA
 * construction and bytecode generation on top of this same front end.
 */

static char *read_file(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "motifql: could not open '%s'\n", path);
        exit(1);
    }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    rewind(f);

    char *buffer = malloc((size_t)size + 1);
    if (!buffer) {
        fprintf(stderr, "motifql: out of memory reading '%s'\n", path);
        fclose(f);
        exit(1);
    }
    size_t read = fread(buffer, 1, (size_t)size, f);
    buffer[read] = '\0';
    fclose(f);
    return buffer;
}

static void dump_tokens(const char *source) {
    printf("== lexical analysis: token stream ==\n");
    Lexer lx;
    lexer_init(&lx, source);
    int count = 0;
    for (;;) {
        Token t = lexer_next_token(&lx);
        printf("%3d  %-10s line %-3d col %-3d  '%.*s'",
               count, token_type_name(t.type), t.line, t.column, t.length, t.start);
        if (t.type == TK_NUMBER) printf("   (value %ld)", t.number);
        printf("\n");
        count++;
        if (t.type == TK_EOF) break;
        if (count > 10000) { printf("  ... token limit reached\n"); break; }
    }
    printf("\n");
}

static void dump_ast(Program *program, bool resolved) {
    printf("== syntax analysis: annotated syntax trees ==\n");
    for (int i = 0; i < program->motif_count; i++) {
        MotifDef *m = &program->motifs[i];
        printf("motif %.*s", m->name_length, m->name);
        if (m->mismatch_budget > 0) printf("  (allow %d mismatch)", m->mismatch_budget);
        printf("\n");
        if (m->pattern) {
            ast_print(m->pattern, 1, program->alphabet, resolved);
        } else {
            printf("  (no pattern -- parse failed for this motif)\n");
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    const char *path = NULL;
    bool want_tokens = false, want_ast = false, want_symbols = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--tokens") == 0)       want_tokens = true;
        else if (strcmp(argv[i], "--ast") == 0)     want_ast = true;
        else if (strcmp(argv[i], "--symbols") == 0) want_symbols = true;
        else if (argv[i][0] == '-') {
            fprintf(stderr, "motifql: unknown option '%s'\n", argv[i]);
            return 1;
        } else {
            path = argv[i];
        }
    }

    if (!path) {
        fprintf(stderr, "usage: motifql [--tokens] [--ast] [--symbols] <file.motif>\n");
        return 1;
    }

    if (!want_tokens && !want_ast && !want_symbols) {
        want_tokens = want_ast = want_symbols = true;
    }

    char *source = read_file(path);

    Diagnostics diag;
    diag_init(&diag, source, path);

    if (want_tokens) dump_tokens(source);

    Program program;
    program_init(&program);

    bool parsed = parse_program(source, &program, &diag);
    bool analyzed = false;
    if (parsed) {
        analyzed = analyze_program(&program, &diag);
    }

    if (want_ast) dump_ast(&program, analyzed);
    if (want_symbols && parsed) symtab_print(&program);

    printf("== front end summary ==\n");
    printf("motifs: %d, searches: %d, errors: %d, warnings: %d\n",
           program.motif_count, program.search_count,
           diag.error_count, diag.warning_count);

    if (!parsed) {
        printf("result: REJECTED during syntax analysis\n");
    } else if (!analyzed) {
        printf("result: REJECTED during semantic analysis\n");
    } else {
        printf("result: ACCEPTED -- ready for Phase 2 (NFA construction)\n");
    }

    int exit_code = (parsed && analyzed) ? 0 : 1;
    program_free(&program);
    free(source);
    return exit_code;
}
