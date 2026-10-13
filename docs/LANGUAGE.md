# MotifQL Language Specification

A domain-specific language for describing biological sequence motifs,
compiled to a finite automaton and executed as matcher bytecode.

## 1. Example

```
alphabet DNA;

motif tata     = "TATAWAWR";
motif kozak    { ("A" | "G") "CC" "ATG" "G" }
motif polyA    { "AATAAA" } allow 1 mismatch;
motif promoter { @tata "N"{20,40} @kozak }

search promoter in "chr1.fa";
```

## 2. Grammar (EBNF)

```
program      := declaration+

declaration  := alphabetDecl | motifDecl | searchStmt

alphabetDecl := "alphabet" IDENT ";"

motifDecl    := "motif" IDENT ( "{" pattern "}" modifiers? require? ";"?
                              | "=" pattern modifiers? require? ";" )

modifiers    := "allow" ( NUMBER | "?" ) ( "mismatch" | "edit" | "edits" )
                ( "at" ( "?" | NUMBER ( "," NUMBER )* ) )?

require      := "require" "rate" "<=" NUMBER "per" NUMBER? UNIT
                ( "," "divergence" NUMBER "%" )?
UNIT         := "bp" | "kb" | "Mb"     (* contextual words, not reserved *)

searchStmt   := "search" IDENT "in" STRING ";"

pattern      := alternation
alternation  := concat ( "|" concat )*
concat       := repeat+
repeat       := atom quantifier?
quantifier   := "*" | "+" | "?"
              | "{" NUMBER ( "," NUMBER? )? "}"
atom         := STRING            (* "ATGN" -- a run of residues   *)
              | class             (* [ACG] or [^T]                 *)
              | "(" pattern ")"
              | "@" IDENT         (* reference to another motif    *)
class        := "[" "^"? IDENT "]"
```

Comments run to end of line and may start with `//` or `#`.

## 3. Alphabets and IUPAC symbols

One `alphabet` declaration per program selects the residue set:

| Alphabet  | Residues               | Wildcard |
|-----------|------------------------|----------|
| `DNA`     | A C G T                | `N`      |
| `RNA`     | A C G U                | `N`      |
| `PROTEIN` | the 20 standard amino acids | `X` |

Nucleotide ambiguity codes are native symbols, not preprocessing:
`R`=AG, `Y`=CT, `S`=GC, `W`=AT, `K`=GT, `M`=AC, `B`=CGT, `D`=AGT,
`H`=ACT, `V`=ACG, `N`=any. Protein supports `X`=any, `B`=DN, `Z`=EQ.

Internally each symbol is a **bitmask over the alphabet's residues**,
so matching a residue against a symbol is a single AND. This is why
ambiguity codes cost nothing at match time and why `[^AG]` in DNA
resolves cleanly to the pyrimidine code `Y`.

**Symbol meaning is alphabet-dependent, and this is the point.** `N`
in DNA is the wildcard {A,C,G,T}; `N` in protein is asparagine, one
specific residue. `U` is legal in RNA and a compile error in DNA.
Resolving symbols during semantic analysis rather than lexing is what
makes that distinction expressible at all.

## 4. Design decisions

### 4.1 Sequence literals are quoted

`"ATG"`, not bare `ATG`. Without quotes, a single-letter IUPAC code
like `N` would be lexically identical to an identifier, forcing the
lexer to track whether it is inside a pattern. Quoting keeps the lexer
context-free and the grammar unambiguous, at the cost of slightly more
typing.

### 4.2 `{` is context-disambiguated, not tokenised separately

`{` opens a motif body after a motif name, and opens a quantifier
after a pattern atom. The parser resolves this by position, so one
`TK_LBRACE` token serves both roles and no lexer mode is needed.

### 4.3 Quantifier bounds are checked semantically, not syntactically

`"ATG"{5,2}` is a *syntactically valid* quantifier with an invalid
bound relationship. Reporting it as a syntax error would misattribute
the problem, so the parser accepts it and semantic analysis rejects
it. The same reasoning applies to residues: the parser records the
characters in a literal without judging them, because whether `U` is
valid depends on a declaration the parser has no business
interpreting.

### 4.4 Motif references may point forward, but never cycle

`@name` resolves against the whole program, so a motif may reference
one defined later in the file. Reference *cycles* are rejected:
`@a -> @b -> @a` describes a pattern of unbounded nesting depth, which
is not a regular language and therefore has no finite automaton. This
check is what guarantees Phase 2's Thompson construction terminates.

### 4.5 A motif may not match the empty sequence

If every quantifier on a motif's top level permits zero repetitions,
the motif matches at every position in the sequence. This is never the
user's intent and would make the scan degenerate, so it is a compile
error rather than a warning.

### 4.6 Mismatch budgets are bounded by the motif's own length

`"ATG" allow 3 mismatch` permits every one of its three residues to
differ, so it matches everywhere. The budget must be strictly less
than the motif's inferred minimum match length.

### 4.7 Two error models: mismatches and edits

`allow k mismatch` counts substitutions only (Hamming distance), so a
hit is always exactly as long as the motif. `allow k edits` counts
substitutions, inserted residues and deleted motif positions
(Levenshtein distance), so a hit may be up to k residues shorter or
longer. In both models an error can only introduce residues of the
declared alphabet.

The empty-match rule extends naturally: under edits, a budget of k
deletions can erase a motif whose shortest form has k residues or
fewer, so `"ATG" allow 3 edits` is an error. When k is at least half
the shortest length the compiler warns instead: the motif is legal,
but after deletions too little of it is left to be specific.

### 4.8 Budgets travel with references

A referenced motif keeps its own budget inside a composite. With
`motif tata = "TATAWAWR" allow 1 edits;`, the composite
`{ @tata "N"{20,40} @kozak }` tolerates one edit inside the TATA box
and none inside the Kozak sequence. Edits count only within a part's
own span: no residue is inserted before a part's first position or
after its last (the spacer already absorbs such residues).

### 4.9 Contracts and chosen budgets

`require rate <= R per N unit` is a statistical contract: the motif may
expect at most R chance hits per N bases. `allow ? edits` leaves a
budget for the compiler to choose. Every free budget must be governed
by exactly one contract, either its own motif's or that of the single
motif that uses it. The compiler then picks the budgets with the
highest sensitivity at the stated divergence that still meet the rate,
or reports that none can. A contract with no free budgets is checked
as written, and breaking it is a compile-time error.

### 4.10 Position-restricted errors

`allow 1 mismatch at 1,2,5` lets errors fall only at the listed
residue positions (1-based, as written). It needs a motif that is a
fixed run of residues, so that "position 3" is unambiguous. With
`allow ? mismatch at ?` the compiler chooses both the budget and the
positions; with `--sites FILE` it scores each choice by how many of the
given known sites it finds, so the positions are learned from examples.

## 5. Inferred match lengths

Semantic analysis computes a `[min, max]` match length for every
node bottom-up, with `max` saturating at infinity for `*` and `+`.
Beyond powering the diagnostics in 4.5 and 4.6, Phase 2 consumes
these: the minimum bounds how far a scan window can advance on a
failed match, and a finite maximum tells the matcher how much
lookahead to buffer over a FASTA stream.

## 6. Out of scope

No back-references (they would make the language non-regular and
break the automaton construction), no capture groups, no strand-aware
reverse-complement matching, and no position-weight-matrix scoring.
Reverse-complement search is the most natural future extension: it
requires no grammar change, only a second automaton built from the
complemented pattern.
