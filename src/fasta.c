#include "fasta.h"

typedef struct {
    char *data;
    long length;
    long capacity;
} DynString;

static void dynstring_push(DynString *ds, char c) {
    if (ds->length + 1 > ds->capacity) {
        long new_cap = ds->capacity == 0 ? 256 : ds->capacity * 2;
        char *grown = realloc(ds->data, (size_t)new_cap);
        if (!grown) { fprintf(stderr, "out of memory reading FASTA\n"); exit(1); }
        ds->data = grown;
        ds->capacity = new_cap;
    }
    ds->data[ds->length++] = c;
}

bool fasta_read(const char *path, FastaFile *out) {
    out->records = NULL;
    out->count = 0;

    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "motifql: could not open FASTA file '%s'\n", path);
        return false;
    }

    int capacity = 0;
    DynString name = {0};
    DynString seq = {0};
    bool have_record = false;

    int c;
    bool at_line_start = true;

    #define FLUSH_RECORD() do {                                         \
        if (have_record) {                                              \
            if (out->count + 1 > capacity) {                            \
                capacity = capacity == 0 ? 4 : capacity * 2;             \
                out->records = realloc(out->records, sizeof(FastaRecord) * (size_t)capacity); \
            }                                                           \
            dynstring_push(&name, '\0');                                \
            dynstring_push(&seq, '\0');                                 \
            out->records[out->count].name = name.data;                  \
            out->records[out->count].sequence = seq.data;                \
            out->records[out->count].length = seq.length - 1;           \
            out->count++;                                              \
            name = (DynString){0};                                      \
            seq = (DynString){0};                                       \
        }                                                                \
    } while (0)

    while ((c = fgetc(f)) != EOF) {
        if (at_line_start && c == '>') {
            FLUSH_RECORD();
            have_record = true;
            /* consume the rest of the header line into `name` */
            while ((c = fgetc(f)) != EOF && c != '\n') {
                if (c != '\r') dynstring_push(&name, (char)c);
            }
            at_line_start = true;
            continue;
        }
        if (c == '\n') { at_line_start = true; continue; }
        at_line_start = false;
        if (c == '\r' || c == ' ' || c == '\t') continue;

        char up = (char)((c >= 'a' && c <= 'z') ? c - 'a' + 'A' : c);
        dynstring_push(&seq, up);
    }
    FLUSH_RECORD();

    #undef FLUSH_RECORD

    free(name.data);
    free(seq.data);
    fclose(f);
    return true;
}

void fasta_free(FastaFile *f) {
    for (int i = 0; i < f->count; i++) {
        free(f->records[i].name);
        free(f->records[i].sequence);
    }
    free(f->records);
    f->records = NULL;
    f->count = 0;
}
