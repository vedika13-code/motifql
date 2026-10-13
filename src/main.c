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
#include "sensitivity.h"
#include "synthesis.h"

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
 *   --markov K                   order of that background, 0..6 (default 1)
 *   --max-rate R                 warn above R expected chance hits per Mb
 *   --merge-contained            drop hits lying wholly inside the previous hit
 *   --count-only                 print hit counts, not every hit
 *
 * Phase 3:
 *   --sensitivity                chance hits AND sensitivity for every motif
 *   --divergence D               substitution rate of diverged sites (default 0.10)
 *   --indel G                    insertion + deletion rate (default 0.02)
 *   --samples N --seed S         Monte Carlo settings (default 5000, 1)
 *   --sites FILE                 known sites (one per record): sensitivity is the
 *                                share found; synthesis learns from them
 *   --empirical FILE             measure chance-hit rates by scanning a reference
 *                                sequence (not the one being searched)
 *   `require` clauses trigger budget synthesis automatically.
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
            nfa_apply_error_tolerance_at(&plain, m->mismatch_budget, m->error_model, m->error_positions, &layered);
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
    if (m->mismatch_budget == BUDGET_FREE) snprintf(buf, n, "? %s", m->error_model == ERR_EDIT ? "edits" : "mismatch");
    else if (m->mismatch_budget == 0) snprintf(buf, n, "exact");
    else snprintf(buf, n, "%d %s", m->mismatch_budget,
                  m->error_model == ERR_EDIT ? (m->mismatch_budget == 1 ? "edit" : "edits")
                                             : (m->mismatch_budget == 1 ? "mismatch" : "mismatches"));
}

/* Compile-time specificity report: no sequence is scanned. */
static void run_specificity(Program *program, Diagnostics *diag, const RateSource *rates, double max_rate) {
    int flagged = 0;
    char srcbuf[256];
    const char *src_text = rate_source_text(rates, srcbuf, sizeof(srcbuf));
    printf("== specificity analysis (%s) ==\n", src_text);
    printf("%-16s %-12s %-6s %-14s %-14s %s\n",
           "motif", "tolerance", "DFA", "P(hit)/pos", "chance/Mb", "chance/100kb");
    for (int i = 0; i < program->motif_count; i++) {
        MotifDef *m = &program->motifs[i];
        MatcherProgram mp;
        CompileStats stats = {0};
        compile_motif_to_bytecode(program, i, &mp, &stats);
        double q = rate_of(&mp, rates);
        char tol[32];
        motif_tolerance_text(m, tol, sizeof(tol));
        printf("%-16.*s %-12s %-6d %-14.3e %-14.1f %.2f\n",
               m->name_length, m->name, tol, stats.dfa_states_final, q, q * 1e6, q * 1e5);
        if (q * 1e6 > max_rate) {
            diag_reportf(diag, PHASE_SPECIFICITY, SEV_WARNING, m->line, m->column,
                         "motif '%.*s' (%s) expects about %.0f chance hits per Mb under the %s "
                         "(threshold %.0f); real sites may be swamped",
                         m->name_length, m->name, tol, q * 1e6, src_text, max_rate);
            flagged++;
        }
        matcher_program_free(&mp);
    }
    printf("motifs above %.0f chance hits per Mb: %d\n\n", max_rate, flagged);
}

/* Specificity and sensitivity side by side: what a motif costs in
 * chance hits, and what it buys in diverged sites found. */
static void run_sensitivity(Program *program, const RateSource *rates, const MutationModel *model,
                            int samples, uint64_t seed, const FastaFile *sites) {
    char srcbuf[256];
    printf("== sensitivity analysis (%s; sites mutated with %g%% substitution, %g%% indel; %d samples, seed %llu) ==\n",
           rate_source_text(rates, srcbuf, sizeof(srcbuf)), model->substitution * 100, model->indel * 100, samples, (unsigned long long)seed);
    printf("%-16s %-12s %-14s %s\n", "motif", "tolerance", "chance/Mb", "sensitivity [95% CI]");
    for (int i = 0; i < program->motif_count; i++) {
        MotifDef *m = &program->motifs[i];
        MatcherProgram mp;
        compile_motif_to_bytecode(program, i, &mp, NULL);
        double q = rate_of(&mp, rates);
        SensitivityResult s = sites ? sensitivity_on_sites(&mp, sites)
                                    : sensitivity_estimate(program, i, &mp, model, samples, seed);
        char tol[32];
        motif_tolerance_text(m, tol, sizeof(tol));
        printf("%-16.*s %-12s %-14.1f %.3f [%.3f, %.3f]\n", m->name_length, m->name, tol, q * 1e6,
               s.sensitivity, s.ci_low, s.ci_high);
        matcher_program_free(&mp);
    }
    printf("\n");
}

static void run_searches(Program *program, bool merge_contained, bool count_only, int search_order) {
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
                if (count_only) { total_hits++; continue; }
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
        if (background_from_fasta(&bg, program->alphabet, &fasta, search_order, path)) {
            long positions = background_valid_positions(&fasta, program->alphabet);
            double q = specificity_hit_probability(&mp, &bg);
            printf("  expected by chance: %.2f (%s, %ld positions)\n", q * (double)positions, bg.source, positions);
            background_free(&bg);
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
    bool count_only = false;
    bool want_specificity = false;
    const char *background_path = NULL;
    const char *empirical_path = NULL;
    const char *sites_path = NULL;
    int markov_order = 1;
    double max_rate = 1000.0;
    bool want_sensitivity = false;
    MutationModel model = { 0.10, 0.02 };
    int samples = 5000;
    uint64_t seed = 1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--tokens") == 0)        want_tokens = true;
        else if (strcmp(argv[i], "--ast") == 0)      want_ast = true;
        else if (strcmp(argv[i], "--symbols") == 0)  want_symbols = true;
        else if (strcmp(argv[i], "--nfa") == 0)      want_nfa = true;
        else if (strcmp(argv[i], "--dfa") == 0)      want_dfa = true;
        else if (strcmp(argv[i], "--bytecode") == 0) want_bytecode = true;
        else if (strcmp(argv[i], "--no-search") == 0) want_search = false;
        else if (strcmp(argv[i], "--merge-contained") == 0) merge_contained = true;
        else if (strcmp(argv[i], "--count-only") == 0) count_only = true;
        else if (strcmp(argv[i], "--specificity") == 0) want_specificity = true;
        else if (strcmp(argv[i], "--background") == 0 && i + 1 < argc) background_path = argv[++i];
        else if (strcmp(argv[i], "--empirical") == 0 && i + 1 < argc) empirical_path = argv[++i];
        else if (strcmp(argv[i], "--sites") == 0 && i + 1 < argc) sites_path = argv[++i];
        else if (strcmp(argv[i], "--markov") == 0 && i + 1 < argc) markov_order = atoi(argv[++i]);
        else if (strcmp(argv[i], "--max-rate") == 0 && i + 1 < argc) max_rate = atof(argv[++i]);
        else if (strcmp(argv[i], "--sensitivity") == 0) want_sensitivity = true;
        else if (strcmp(argv[i], "--divergence") == 0 && i + 1 < argc) model.substitution = atof(argv[++i]);
        else if (strcmp(argv[i], "--indel") == 0 && i + 1 < argc) model.indel = atof(argv[++i]);
        else if (strcmp(argv[i], "--samples") == 0 && i + 1 < argc) samples = atoi(argv[++i]);
        else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) seed = (uint64_t)strtoull(argv[++i], NULL, 10);
        else if (argv[i][0] == '-') {
            fprintf(stderr, "motifql: unknown option '%s'\n", argv[i]);
            return 1;
        } else {
            path = argv[i];
        }
    }

    if (!path) {
        fprintf(stderr, "usage: motifql [--tokens] [--ast] [--symbols] [--nfa] [--dfa] [--bytecode]\n"
                        "               [--specificity [--background FILE] [--markov K] [--max-rate R]]\n"
                        "               [--sensitivity [--divergence D] [--indel G] [--samples N] [--seed S]]\n"
                        "               [--empirical REFERENCE.fa] [--sites KNOWN.fa] [--count-only]\n"
                        "               [--merge-contained] [--no-search] <file.motif>\n");
        return 1;
    }

    if (!want_tokens && !want_ast && !want_symbols && !want_nfa && !want_dfa && !want_bytecode &&
        !want_specificity && !want_sensitivity) {
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

        /* One background for every statistical pass: uniform, or
         * estimated from --background FILE. */
        Background bg;
        background_uniform(&bg, program.alphabet);
        if (background_path) {
            FastaFile bgf;
            if (fasta_read(background_path, &bgf)) {
                Background fitted;
                if (background_from_fasta(&fitted, program.alphabet, &bgf, markov_order, background_path)) {
                    background_free(&bg);
                    bg = fitted;
                } else {
                    fprintf(stderr, "motifql: no usable residues in '%s'; using a uniform background\n", background_path);
                }
                fasta_free(&bgf);
            }
        }

        /* --empirical FILE: measure chance-hit rates by scanning a
         * reference sequence instead of using the background model. */
        FastaFile empirical_fasta = {0};
        RateSource rates = { &bg, NULL, NULL };
        if (empirical_path && fasta_read(empirical_path, &empirical_fasta)) {
            rates.reference = &empirical_fasta;
            rates.reference_label = empirical_path;
        }

        /* --sites FILE: known sites, one per record; sensitivity becomes
         * the share of them the motif finds. */
        FastaFile sites_fasta = {0};
        bool sites_loaded = sites_path && fasta_read(sites_path, &sites_fasta);

        /* Budget synthesis settles every `allow ?` before anything is
         * compiled for output or search. */
        bool synth_ok = true;
        for (int i = 0; i < program.motif_count; i++) {
            if (program.motifs[i].has_requirement) {
                SynthesisOptions so = { &bg, rates, model, samples, seed,
                                        sites_loaded ? &sites_fasta : NULL, sites_path };
                synth_ok = synthesize_budgets(&program, &diag, &so);
                break;
            }
        }
        if (!synth_ok) {
            printf("result: REJECTED during budget synthesis\n");
            background_free(&bg);
            if (rates.reference) fasta_free(&empirical_fasta);
            if (sites_loaded) fasta_free(&sites_fasta);
            program_free(&program);
            free(source);
            return 1;
        }

        if (want_nfa || want_dfa || want_bytecode) {
            dump_nfa_dfa_bytecode(&program, want_nfa, want_dfa, want_bytecode);
        }
        if (want_specificity) run_specificity(&program, &diag, &rates, max_rate);
        if (want_sensitivity) run_sensitivity(&program, &rates, &model, samples, seed,
                                              sites_loaded ? &sites_fasta : NULL);
        if (want_search) {
            run_searches(&program, merge_contained, count_only, markov_order);
        }
        background_free(&bg);
        if (rates.reference) fasta_free(&empirical_fasta);
        if (sites_loaded) fasta_free(&sites_fasta);
    }

    int exit_code = (parsed && analyzed) ? 0 : 1;
    program_free(&program);
    free(source);
    return exit_code;
}
