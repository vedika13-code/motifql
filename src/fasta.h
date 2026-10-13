#ifndef MOTIFQL_FASTA_H
#define MOTIFQL_FASTA_H

#include "common.h"

typedef struct {
    char *name;       /* header line, without '>' */
    char *sequence;   /* concatenated, whitespace-stripped residues, uppercased */
    long length;
} FastaRecord;

typedef struct {
    FastaRecord *records;
    int count;
} FastaFile;

/* Reads and parses a FASTA file. Returns false (with a message on
 * stderr) if the file cannot be opened; a file that opens but
 * contains no ">" header still yields a zero-record FastaFile, which
 * is reported as-is rather than treated as an error, since an empty
 * input is not a malformed one. */
bool fasta_read(const char *path, FastaFile *out);
void fasta_free(FastaFile *f);

#endif
