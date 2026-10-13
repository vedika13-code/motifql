#include "symtab.h"

void program_init(Program *program) {
    memset(program, 0, sizeof(*program));
    program->alphabet = ALPHA_DNA;
}

void program_free(Program *program) {
    for (int i = 0; i < program->motif_count; i++) {
        node_free(program->motifs[i].pattern);
        program->motifs[i].pattern = NULL;
    }
    program->motif_count = 0;
    program->search_count = 0;
}

int program_find_motif(Program *program, const char *name, int name_length) {
    for (int i = 0; i < program->motif_count; i++) {
        MotifDef *m = &program->motifs[i];
        if (m->name_length == name_length && memcmp(m->name, name, (size_t)name_length) == 0) {
            return i;
        }
    }
    return -1;
}

int program_add_motif(Program *program, const char *name, int name_length, int line, int column) {
    if (program->motif_count >= MAX_MOTIFS) return -1;
    MotifDef *m = &program->motifs[program->motif_count];
    memset(m, 0, sizeof(*m));
    m->name = name;
    m->name_length = name_length;
    m->line = line;
    m->column = column;
    m->min_len = 0;
    m->max_len = 0;
    return program->motif_count++;
}

void symtab_print(Program *program) {
    printf("== symbol table ==\n");
    printf("alphabet: %s\n", program->alphabet_declared
        ? alphabet_name(program->alphabet)
        : "(not declared, defaulted to DNA)");
    printf("%-4s %-16s %-10s %-10s %-9s %s\n",
           "#", "motif", "tolerance", "len", "used", "defined at");
    for (int i = 0; i < program->motif_count; i++) {
        MotifDef *m = &program->motifs[i];
        char lenbuf[32];
        if (m->max_len == LEN_INF) {
            snprintf(lenbuf, sizeof(lenbuf), "%d..inf", m->min_len);
        } else {
            snprintf(lenbuf, sizeof(lenbuf), "%d..%d", m->min_len, m->max_len);
        }
        char tolbuf[24];
        if (m->mismatch_budget == BUDGET_FREE) snprintf(tolbuf, sizeof(tolbuf), "? %s",
                      m->error_model == ERR_EDIT ? "edit" : "mism");
        else if (m->mismatch_budget == 0) snprintf(tolbuf, sizeof(tolbuf), "exact");
        else snprintf(tolbuf, sizeof(tolbuf), "%d %s", m->mismatch_budget,
                      m->error_model == ERR_EDIT ? "edit" : "mism");
        printf("%-4d %-16.*s %-10s %-10s %-9s line %d\n",
               i, m->name_length, m->name, tolbuf,
               lenbuf, m->used ? "yes" : "no", m->line);
    }

    if (program->search_count > 0) {
        printf("\nsearches:\n");
        for (int i = 0; i < program->search_count; i++) {
            SearchStmt *s = &program->searches[i];
            printf("  %.*s in \"%.*s\"  (motif #%d, line %d)\n",
                   s->motif_name_length, s->motif_name,
                   s->path_length, s->path,
                   s->motif_index, s->line);
        }
    }
    printf("\n");
}
