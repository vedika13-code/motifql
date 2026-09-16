#ifndef MOTIFQL_ALPHABET_H
#define MOTIFQL_ALPHABET_H

#include "common.h"

/*
 * A symbol in a motif is stored as a BITMASK over the declared
 * alphabet's residues, not as a character. This is what makes IUPAC
 * ambiguity codes a first-class concept rather than a preprocessing
 * step: 'N' in DNA is simply the mask {A,C,G,T}, and 'R' is {A,G}.
 * Matching a residue against a symbol is then a single AND, which
 * matters once Phase 2 builds an NFA out of these.
 */

typedef enum {
    ALPHA_DNA,
    ALPHA_RNA,
    ALPHA_PROTEIN,
} AlphabetKind;

typedef uint32_t SymMask;

/* Number of distinct residues in the alphabet (4 for DNA/RNA, 20 for
 * protein). Determines how wide a DFA transition row is in Phase 2. */
int alphabet_size(AlphabetKind kind);

const char *alphabet_name(AlphabetKind kind);

/* Ordered residue characters, e.g. "ACGT". Index i corresponds to bit i. */
const char *alphabet_residues(AlphabetKind kind);

/* Resolves one source character to its residue mask under `kind`.
 * Returns 0 if the character is not valid in this alphabet -- note
 * that validity is alphabet-dependent by design: 'U' is legal in RNA
 * but not DNA, and 'N' means "any base" in DNA but is asparagine (a
 * single specific residue) in protein. Catching that difference is
 * the point of doing this in semantic analysis rather than the lexer. */
SymMask alphabet_symbol_mask(AlphabetKind kind, char c);

/* Renders a mask back to readable form for AST dumps: a single
 * residue prints bare ("A"), an ambiguity code prints as its IUPAC
 * letter when one exists ("N"), otherwise as a set ("[ACG]"). */
void alphabet_format_mask(AlphabetKind kind, SymMask mask, char *buf, size_t buf_size);

/* Parses an alphabet name from a declaration. Returns false if unknown. */
bool alphabet_from_name(const char *name, int length, AlphabetKind *out);

#endif
