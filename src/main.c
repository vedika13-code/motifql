#include "common.h"
#include "lexer.h"
#include "parser.h"
#include "sema.h"
#include "symtab.h"
#include "diagnostics.h"
#include "nfa.h"
#include "dfa.h"
#include "codegen.h"
#include "fasta.h"
#include "matcher.h"
#include "specificity.h"

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
 * With no flags, all three are printed. Phase 2 adds, behind the same
 * front end:
 *
 *   --nfa / --dfa / --bytecode   the automaton at each stage
 *   --specificity                expected chance hits for every motif,
 *                                computed from the compiled automaton
 *   --background FILE            estimate the background from a FASTA file
 *   --markov 0|1                 order of that background (default 1)
 *   --max-rate R                 warn above R expected chance hits per Mb
 *   --merge-contained            drop hits lying wholly inside the previous hit
 *   --no-search                  compile only
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
        if (m->mismatch_budget > 0)
            printf("  (allow %d %s)", m->mismatch_budget,
                   m->error_model == ERR_EDIT ? "edits" : "mismatch");
        printf("\n");
        if (m->pattern) {
            ast_print(m->pattern, 1, program->alphabet, resolved);
        } else {
            printf("  (no pattern -- parse failed for this motif)\n");
        }
        printf("\n");
    }
}

static void dump_nfa_dfa_bytecode(Program *program, bool want_nfa, bool want_dfa, bool want_bytecode) {
    for (int i = 0; i < program->motif_count; i++) {
        MotifDef *m = &program->motifs[i];

        Nfa plain;
        nfa_build_for_motif(program, i, &plain);

        Nfa *final_nfa = &plain;
        Nfa layered;
        bool built_layered = false;
        if (m->mismatch_budget > 0) {
            nfa_apply_error_tolerance(&plain, m->mismatch_budget, m->error_model, &layered);
            final_nfa = &layered;
            built_layered = true;
        }

        if (want_nfa) {
            printf("motif %.*s:\n", m->name_length, m->name);
            nfa_print(final_nfa);
            printf("\n");
        }

        if (want_dfa) {
            Dfa dfa;
            dfa_from_nfa(final_nfa, &dfa);
            int raw = dfa.count;
            int dead = dfa_eliminate_dead_states(&dfa);
            int before_min = dfa.count;
            dfa_minimize(&dfa);
            printf("motif %.*s:\n", m->name_length, m->name);
            printf("subset construction: %d states\n", raw);
            printf("dead-state elimination: removed %d, %d states remain\n", dead, before_min);
            printf("minimization: removed %d, %d states remain\n", before_min - dfa.count, dfa.count);
            dfa_print(&dfa);
            printf("\n");
            dfa_free(&dfa);
        }

        if (want_bytecode) {
            MatcherProgram mp;
            CompileStats stats = {0};
            compile_motif_to_bytecode(program, i, &mp, &stats);
            printf("motif %.*s: NFA %d -> DFA %d -> (dead-elim) %d -> (minimized) %d states\n",
                   m->name_length, m->name, stats.nfa_states, stats.dfa_states_raw,
                   stats.dfa_states_after_dead_elim, stats.dfa_states_final);
            matcher_program_print(&mp);
            matcher_program_free(&mp);
        }

        if (built_layered) nfa_free(&layered);
        nfa_free(&plain);
    }
}

static void motif_tolerance_text(const MotifDef *m, char *buf, size_t n) {
    if (m->mismatch_budget == 0) snprintf(buf, n, "exact");
    else snprintf(buf, n, "%d %s", m->mismatch_budget,
                  m->error_model == ERR_EDIT ? (m->mismatch_budget == 1 ? "edit" : "edits")
                                             : (m->mismatch_budget == 1 ? "mismatch" : "mismatches"));
}

/* Compile-time specificity report: no sequence is scanned. */
static void run_specificity(Program *program, Diagnostics *diag, const Background *bg, double max_rate) {
    int flagged = 0;
    printf("== specificity analysis (%s) ==\n", bg->source);
    printf("%-16s %-12s %-6s %-14s %-14s %s\n",
           "motif", "tolerance", "DFA", "P(hit)/pos", "chance/Mb", "chance/100kb");
    for (int i = 0; i < program->motif_count; i++) {
        MotifDef *m = &program->motifs[i];
        MatcherProgram mp;
        CompileStats stats = {0};
        compile_motif_to_bytecode(program, i, &mp, &stats);
        double q = specificity_hit_probability(&mp, bg);
        char tol[32];
        motif_tolerance_text(m, tol, sizeof(tol));
        printf("%-16.*s %-12s %-6d %-14.3e %-14.1f %.2f\n",
               m->name_length, m->name, tol, stats.dfa_states_final, q, q * 1e6, q * 1e5);
        if (q * 1e6 > max_rate) {
            diag_reportf(diag, PHASE_SPECIFICITY, SEV_WARNING, m->line, m->column,
                         "motif '%.*s' (%s) expects about %.0f chance hits per Mb under the %s "
                         "(threshold %.0f); real sites may be swamped",
                         m->name_length, m->name, tol, q * 1e6, bg->source, max_rate);
            flagged++;
        }
        matcher_program_free(&mp);
    }
    printf("motifs above %.0f chance hits per Mb: %d\n\n", max_rate, flagged);
}

static void run_searches(Program *program, bool merge_contained) {
    if (program->search_count == 0) return;

    printf("== running searches ==\n");
    for (int i = 0; i < program->search_count; i++) {
        SearchStmt *s = &program->searches[i];
        if (s->motif_index < 0) continue;

        char path[512];
        int plen = s->path_length < 511 ? s->path_length : 511;
        memcpy(path, s->path, (size_t)plen);
        path[plen] = '\0';

        FastaFile fasta;
        if (!fasta_read(path, &fasta)) continue;

        MatcherProgram mp;
        CompileStats stats = {0};
        compile_motif_to_bytecode(program, s->motif_index, &mp, &stats);

        printf("%.*s in \"%s\": %d record(s), automaton %d NFA -> %d DFA states, prefix \"%s\"\n",
               s->motif_name_length, s->motif_name, path, fasta.count,
               stats.nfa_states, stats.dfa_states_final, mp.literal_prefix);

        int total_hits = 0;
        for (int r = 0; r < fasta.count; r++) {
            FastaRecord *rec = &fasta.records[r];
            HitList hits;
            matcher_scan(&mp, rec->sequence, rec->length, &hits);
            long kept_start = -1, kept_end = -1;
            for (int h = 0; h < hits.count; h++) {
                long hs = hits.hits[h].start, he = hs + hits.hits[h].length;
                /* --merge-contained: with edits, one real site also
                 * produces shifted copies whose span lies inside it
                 * (e.g. the same site minus its first residue, paid
                 * for by a deletion). Drop a hit wholly contained in
                 * the last one kept. */
                if (merge_contained && kept_start >= 0 && hs >= kept_start && he <= kept_end) continue;
                kept_start = hs; kept_end = he;
                printf("  %s:%ld-%ld (len %ld)  %.*s\n",
                       rec->name, hits.hits[h].start + 1, hits.hits[h].start + hits.hits[h].length,
                       hits.hits[h].length,
                       (int)hits.hits[h].length, rec->sequence + hits.hits[h].start);
                total_hits++;
            }
            hitlist_free(&hits);
        }
        printf("  -> %d hit(s)\n", total_hits);

        /* Put the observed count next to what chance alone would give
         * on a random sequence with this file's own composition. */
        Background bg;
        if (background_from_fasta(&bg, program->alphabet, &fasta, 1, path)) {
            long positions = background_valid_positions(&fasta, program->alphabet);
            double q = specificity_hit_probability(&mp, &bg);
            printf("  expected by chance: %.2f (%s, %ld positions)\n", q * (double)positions, bg.source, positions);
        }

        matcher_program_free(&mp);
        fasta_free(&fasta);
    }
    printf("\n");
}

int main(int argc, char **argv) {
    const char *path = NULL;
    bool want_tokens = false, want_ast = false, want_symbols = false;
    bool want_nfa = false, want_dfa = false, want_bytecode = false;
    bool want_search = true;
    bool merge_contained = false;
    bool want_specificity = false;
    const char *background_path = NULL;
    int markov_order = 1;
    double max_rate = 1000.0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--tokens") == 0)        want_tokens = true;
        else if (strcmp(argv[i], "--ast") == 0)      want_ast = true;
        else if (strcmp(argv[i], "--symbols") == 0)  want_symbols = true;
        else if (strcmp(argv[i], "--nfa") == 0)      want_nfa = true;
        else if (strcmp(argv[i], "--dfa") == 0)      want_dfa = true;
        else if (strcmp(argv[i], "--bytecode") == 0) want_bytecode = true;
        else if (strcmp(argv[i], "--no-search") == 0) want_search = false;
        else if (strcmp(argv[i], "--merge-contained") == 0) merge_contained = true;
        else if (strcmp(argv[i], "--specificity") == 0) want_specificity = true;
        else if (strcmp(argv[i], "--background") == 0 && i + 1 < argc) background_path = argv[++i];
        else if (strcmp(argv[i], "--markov") == 0 && i + 1 < argc) markov_order = atoi(argv[++i]) ? 1 : 0;
        else if (strcmp(argv[i], "--max-rate") == 0 && i + 1 < argc) max_rate = atof(argv[++i]);
        else if (argv[i][0] == '-') {
            fprintf(stderr, "motifql: unknown option '%s'\n", argv[i]);
            return 1;
        } else {
            path = argv[i];
        }
    }

    if (!path) {
        fprintf(stderr, "usage: motifql [--tokens] [--ast] [--symbols] [--nfa] [--dfa] [--bytecode]\n"
                        "               [--specificity [--background FILE] [--markov 0|1] [--max-rate R]]\n"
                        "               [--merge-contained] [--no-search] <file.motif>\n");
        return 1;
    }

    if (!want_tokens && !want_ast && !want_symbols && !want_nfa && !want_dfa && !want_bytecode && !want_specificity) {
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
        printf("result: ACCEPTED\n\n");
        if (want_nfa || want_dfa || want_bytecode) {
            dump_nfa_dfa_bytecode(&program, want_nfa, want_dfa, want_bytecode);
        }
        if (want_specificity) {
            Background bg;
            background_uniform(&bg, program.alphabet);
            if (background_path) {
                FastaFile bgf;
                if (fasta_read(background_path, &bgf)) {
                    if (!background_from_fasta(&bg, program.alphabet, &bgf, markov_order, background_path))
                        fprintf(stderr, "motifql: no usable residues in '%s'; using a uniform background\n", background_path);
                    fasta_free(&bgf);
                }
            }
            run_specificity(&program, &diag, &bg, max_rate);
        }
        if (want_search) {
            run_searches(&program, merge_contained);
        }
    }

    int exit_code = (parsed && analyzed) ? 0 : 1;
    program_free(&program);
    free(source);
    return exit_code;
}
