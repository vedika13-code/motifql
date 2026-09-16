#include "alphabet.h"

#define DNA_RESIDUES     "ACGT"
#define RNA_RESIDUES     "ACGU"
#define PROTEIN_RESIDUES "ACDEFGHIKLMNPQRSTVWY"

int alphabet_size(AlphabetKind kind) {
    return (int)strlen(alphabet_residues(kind));
}

const char *alphabet_name(AlphabetKind kind) {
    switch (kind) {
        case ALPHA_DNA:     return "DNA";
        case ALPHA_RNA:     return "RNA";
        case ALPHA_PROTEIN: return "PROTEIN";
    }
    return "?";
}

const char *alphabet_residues(AlphabetKind kind) {
    switch (kind) {
        case ALPHA_DNA:     return DNA_RESIDUES;
        case ALPHA_RNA:     return RNA_RESIDUES;
        case ALPHA_PROTEIN: return PROTEIN_RESIDUES;
    }
    return "";
}

bool alphabet_from_name(const char *name, int length, AlphabetKind *out) {
    if (length == 3 && memcmp(name, "DNA", 3) == 0)         { *out = ALPHA_DNA;     return true; }
    if (length == 3 && memcmp(name, "RNA", 3) == 0)         { *out = ALPHA_RNA;     return true; }
    if (length == 7 && memcmp(name, "PROTEIN", 7) == 0)     { *out = ALPHA_PROTEIN; return true; }
    return false;
}

static SymMask bit_for(const char *residues, char c) {
    const char *p = strchr(residues, c);
    if (!p) return 0;
    return (SymMask)1u << (int)(p - residues);
}

/* Nucleotide IUPAC ambiguity codes, shared between DNA and RNA. The
 * table is written in terms of A/C/G/T and remapped for RNA below,
 * because the codes are defined identically apart from T<->U. */
typedef struct { char code; const char *expands_to; } IupacEntry;

static const IupacEntry NUC_IUPAC[] = {
    {'R', "AG"},   {'Y', "CT"},   {'S', "GC"},   {'W', "AT"},
    {'K', "GT"},   {'M', "AC"},   {'B', "CGT"},  {'D', "AGT"},
    {'H', "ACT"},  {'V', "ACG"},  {'N', "ACGT"},
    {0, NULL}
};

/* Protein ambiguity codes. 'X' is any residue; B and Z are the two
 * classical ambiguity pairs. Note there is deliberately no 'N' entry
 * here -- in protein, N is asparagine, a concrete residue. */
static const IupacEntry PROT_IUPAC[] = {
    {'X', PROTEIN_RESIDUES},
    {'B', "DN"},
    {'Z', "EQ"},
    {0, NULL}
};

static SymMask nucleotide_mask(AlphabetKind kind, char c) {
    const char *residues = alphabet_residues(kind);
    char thymine = (kind == ALPHA_RNA) ? 'U' : 'T';

    SymMask direct = bit_for(residues, c);
    if (direct) return direct;

    for (int i = 0; NUC_IUPAC[i].code; i++) {
        if (NUC_IUPAC[i].code != c) continue;
        SymMask mask = 0;
        for (const char *e = NUC_IUPAC[i].expands_to; *e; e++) {
            char residue = (*e == 'T') ? thymine : *e;
            mask |= bit_for(residues, residue);
        }
        return mask;
    }
    return 0;
}

static SymMask protein_mask(char c) {
    SymMask direct = bit_for(PROTEIN_RESIDUES, c);
    if (direct) return direct;

    for (int i = 0; PROT_IUPAC[i].code; i++) {
        if (PROT_IUPAC[i].code != c) continue;
        SymMask mask = 0;
        for (const char *e = PROT_IUPAC[i].expands_to; *e; e++) {
            mask |= bit_for(PROTEIN_RESIDUES, *e);
        }
        return mask;
    }
    return 0;
}

SymMask alphabet_symbol_mask(AlphabetKind kind, char c) {
    if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
    if (kind == ALPHA_PROTEIN) return protein_mask(c);
    return nucleotide_mask(kind, c);
}

void alphabet_format_mask(AlphabetKind kind, SymMask mask, char *buf, size_t buf_size) {
    const char *residues = alphabet_residues(kind);
    int n = alphabet_size(kind);

    int count = 0;
    for (int i = 0; i < n; i++) if (mask & ((SymMask)1u << i)) count++;

    if (count == 1) {
        for (int i = 0; i < n; i++) {
            if (mask & ((SymMask)1u << i)) {
                snprintf(buf, buf_size, "%c", residues[i]);
                return;
            }
        }
    }

    /* Prefer the IUPAC letter when the mask corresponds to one exactly. */
    const IupacEntry *table = (kind == ALPHA_PROTEIN) ? PROT_IUPAC : NUC_IUPAC;
    char thymine = (kind == ALPHA_RNA) ? 'U' : 'T';
    for (int i = 0; table[i].code; i++) {
        SymMask candidate = 0;
        for (const char *e = table[i].expands_to; *e; e++) {
            char residue = *e;
            if (kind != ALPHA_PROTEIN && residue == 'T') residue = thymine;
            candidate |= bit_for(residues, residue);
        }
        if (candidate == mask) {
            snprintf(buf, buf_size, "%c", table[i].code);
            return;
        }
    }

    /* Fall back to an explicit set. */
    size_t pos = 0;
    if (pos + 1 < buf_size) buf[pos++] = '[';
    for (int i = 0; i < n && pos + 1 < buf_size; i++) {
        if (mask & ((SymMask)1u << i)) buf[pos++] = residues[i];
    }
    if (pos + 1 < buf_size) buf[pos++] = ']';
    buf[pos] = '\0';
}
