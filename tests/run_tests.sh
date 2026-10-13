#!/usr/bin/env bash
# MotifQL test suite: front end (Phase 1), automaton and matcher,
# edit distance, compile-time specificity analysis (Phase 2), and
# per-part budgets, sensitivity and budget synthesis (Phase 3).
#
# Each case asserts only the ACCEPT/REJECT verdict, not exact message
# text, so error wording can be improved without breaking the suite.
# The messages themselves are checked by eye via `make demo`.

BIN=./build/motifql
PASS=0
FAIL=0

check() {
    local file="$1"
    local expected="$2"     # "accept" or "reject"
    local label="$3"

    $BIN --symbols "$file" > /dev/null 2>&1
    local code=$?

    local actual="accept"
    [ $code -ne 0 ] && actual="reject"

    if [ "$actual" = "$expected" ]; then
        printf "  PASS  %-22s %s\n" "$expected" "$label"
        PASS=$((PASS + 1))
    else
        printf "  FAIL  expected %-13s got %s -- %s\n" "$expected" "$actual" "$label"
        FAIL=$((FAIL + 1))
    fi
}

echo "== MotifQL Phase 1 front-end tests =="

check tests/promoter.motif        accept  "DNA promoter with composition and IUPAC codes"
check tests/protein.motif         accept  "protein alphabet, N is asparagine not wildcard"
check tests/rna.motif             accept  "RNA alphabet accepts uracil"

check tests/semantic_errors.motif reject  "invalid symbol, bounds, dangling ref, empty match"
check tests/cycle.motif           reject  "indirect reference cycle a -> b -> c -> a"
check tests/self_cycle.motif      reject  "direct self-reference"
check tests/syntax_errors.motif   reject  "syntax and lexical errors with recovery"
check tests/cross_alphabet.motif  reject  "DNA rejects a protein-only residue"

echo
echo "== MotifQL Phase 2 automaton / matcher tests =="

check_output_contains() {
    local file="$1"
    local needle="$2"
    local label="$3"
    if (cd tests && ../build/motifql --symbols "$(basename "$file")" 2>/dev/null | grep -qF "$needle"); then
        printf "  PASS  %s\n" "$label"
        PASS=$((PASS + 1))
    else
        printf "  FAIL  expected output containing: %s -- %s\n" "$needle" "$label"
        FAIL=$((FAIL + 1))
    fi
}

# Exact literal, hand-verified position (GAATTC planted at 1-based offset 13).
check_output_contains tests/restriction_sites.motif \
    "13-18 (len 6)  GAATTC" \
    "EcoRI site found at its known planted position"

# IUPAC class in the middle of a literal (GC-N-GC), hand-verified position.
check_output_contains tests/restriction_sites.motif \
    "27-31 (len 5)  GCAGC" \
    "degenerate GC-N-GC site resolved and found at its planted position"

# Composed motif via @references and a bounded-length spacer,
# hand-verified position and exact matched length.
check_output_contains tests/promoter.motif \
    "17-56 (len 40)  TATAAAAGGGGCCCCGGGGCCCCAAAAAAAAAAGCCATGG" \
    "composite promoter (TATA + spacer + Kozak) found at its planted position"

# Mismatch tolerance: the exact site should still be found...
check_output_contains tests/promoter.motif \
    "78-83 (len 6)  AATAAA" \
    "polyA signal found as an exact match under 'allow 1 mismatch'"

# ...and a real Hamming-distance-1 near match should ALSO be found,
# confirming the mismatch budget is doing real work, not being ignored.
check_output_contains tests/promoter.motif \
    "40-45 (len 6)  AAAAAA" \
    "polyA signal's 1-mismatch tolerance also matches a Hamming-1 near-miss"

# Literal-prefix extraction: a fully-literal motif's prefix should be
# the whole motif; a motif starting with an ambiguity code should
# extract no prefix at all.
check_output_contains tests/restriction_sites.motif \
    'prefix "GAATTC"' \
    "literal-prefix extraction recovers the whole motif when it is fully literal"

echo
echo "== Edit distance (alphabet-typed Levenshtein) =="

check tests/edit_errors.motif     reject  "3 edits on a 3-residue motif could match the empty sequence"
check tests/edit_syntax.motif     reject  "unknown tolerance word after 'allow k' is a syntax error"

check_stderr_contains() {
    local file="$1"; local needle="$2"; local label="$3"; shift 3
    if (cd tests && ../build/motifql --no-search "$@" "$(basename "$file")" 2>&1 >/dev/null | grep -qF "$needle"); then
        printf "  PASS  %s\n" "$label"; PASS=$((PASS + 1))
    else
        printf "  FAIL  expected diagnostic containing: %s -- %s\n" "$needle" "$label"; FAIL=$((FAIL + 1))
    fi
}

check_stderr_contains tests/edit_errors.motif "only 2 residue(s) would remain" \
    "loose-but-legal edit budget produces a warning, not an error" --symbols

# Planted variants in edits.fa (1-based): exact site at 11, a one-base
# DELETION variant at 29, a one-base INSERTION variant at 46, a
# substitution variant at 65.
check_output_contains tests/edits.motif "29-35 (len 7)  TATAAAG" \
    "deletion variant found only under 'allow 1 edits' (match is 7 long)"
check_output_contains tests/edits.motif "46-54 (len 9)  TATACAAAG" \
    "insertion variant found under 'allow 1 edits' (match is 9 long)"
check_output_contains tests/edits.motif "65-72 (len 8)  TATGAAAG" \
    "substitution variant found by both mismatch and edit tolerance"

if (cd tests && ../build/motifql --symbols edits.motif | sed -n '/tata_mm in/,/hit(s)/p' | grep -q "29-35"); then
    printf "  FAIL  mismatch-only tolerance wrongly matched the deletion variant\n"; FAIL=$((FAIL + 1))
else
    printf "  PASS  mismatch-only tolerance does not match the deletion variant\n"; PASS=$((PASS + 1))
fi

if command -v python3 >/dev/null 2>&1; then
    if python3 tools/oracle_check.py 200 42 > /tmp/motifql_oracle.txt 2>&1; then
        printf "  PASS  %s\n" "$(tail -1 /tmp/motifql_oracle.txt) (brute-force oracle)"; PASS=$((PASS + 1))
    else
        printf "  FAIL  brute-force oracle disagreed:\n"; cat /tmp/motifql_oracle.txt; FAIL=$((FAIL + 1))
    fi
fi

echo
echo "== Compile-time specificity analysis =="

check_spec_contains() {
    local file="$1"; local needle="$2"; local label="$3"
    if (cd tests && ../build/motifql --specificity --no-search "$(basename "$file")" 2>/dev/null | grep -qF "$needle"); then
        printf "  PASS  %s\n" "$label"; PASS=$((PASS + 1))
    else
        printf "  FAIL  expected specificity output containing: %s -- %s\n" "$needle" "$label"; FAIL=$((FAIL + 1))
    fi
}

# Closed-form checks under a uniform background:
#   TATAWAWR exact: (1/4)^5 * (1/2)^3 = 1/8192 -> 122.1 per Mb
#   AATAAA, 1 mismatch: (1 + 6*3) / 4^6 = 19/4096 -> 4638.7 per Mb
#   [^AG]{3,}: first three residues must be C/T -> 1/8 -> 125000 per Mb
check_spec_contains tests/promoter.motif "tata             exact        9      1.221e-04      122.1" \
    "exact IUPAC motif: estimate equals the closed form 1/8192"
check_spec_contains tests/promoter.motif "polyA            1 mismatch   12     4.639e-03      4638.7" \
    "1-mismatch motif: estimate equals the closed form 19/4096"
check_spec_contains tests/promoter.motif "pyrimidine_run   exact        4      1.250e-01      125000.0" \
    "unbounded repeat: estimate equals the closed form 1/8"
check_stderr_contains tests/promoter.motif "expects about 4639 chance hits per Mb" \
    "over-permissive motif triggers a specificity warning at compile time" --specificity

echo
echo "== Phase 3: per-part budgets, sensitivity, budget synthesis =="

check_output_contains tests/perpart.motif "70-109 (len 40)  TATAAAGGGGG" \
    "per-part budget: TATA deletion inside a composite is found"
if (cd tests && ../build/motifql --symbols perpart.motif | sed -n '/^strict in/,/hit(s)/p' | grep -q "70-109"); then
    printf "  FAIL  composite without a part budget wrongly matched the TATA deletion\n"; FAIL=$((FAIL + 1))
else
    printf "  PASS  composite without a part budget does not match the TATA deletion\n"; PASS=$((PASS + 1))
fi
if (cd tests && ../build/motifql --symbols perpart.motif | sed -n '/^promoter in/,/hit(s)/p' | grep -q "129-"); then
    printf "  FAIL  Kozak deletion matched although Kozak has no budget\n"; FAIL=$((FAIL + 1))
else
    printf "  PASS  a part's budget stays local: Kozak deletion is not matched\n"; PASS=$((PASS + 1))
fi

# Sensitivity, closed form with substitutions only: TATAWAWR keeps its
# 5 fixed positions with p = 0.9 each and its 3 two-base positions with
# p = 0.9 + 0.1/3 each: 0.9^5 * (0.9333...)^3 = 0.4801.
if (cd tests && ../build/motifql --sensitivity --no-search --indel 0 --samples 20000 edits.motif 2>/dev/null \
      | awk '/^tata_exact/ { gsub(/[\[\],]/, " "); if ($5 <= 0.4801 && $6 >= 0.4801) ok=1 } END { exit !ok }'); then
    printf "  PASS  sensitivity 95%% interval contains the closed form 0.4801\n"; PASS=$((PASS + 1))
else
    printf "  FAIL  sensitivity interval misses the closed form 0.4801\n"; FAIL=$((FAIL + 1))
fi

check_synth() {
    local file="$1"; local needle="$2"; local label="$3"
    if (cd tests && ../build/motifql --symbols "$(basename "$file")" 2>&1 | grep -qF "$needle"); then
        printf "  PASS  %s\n" "$label"; PASS=$((PASS + 1))
    else
        printf "  FAIL  expected output containing: %s -- %s\n" "$needle" "$label"; FAIL=$((FAIL + 1))
    fi
}
check_synth tests/synthesis.motif "note: chose tata = 1 edit, kozak = 0 edits" \
    "synthesis picks the most sensitive budgets that meet 2 per 100kb"
check_synth tests/synthesis.motif "10 skipped by monotone pruning" \
    "monotone pruning skips assignments above an infeasible one"
check_synth tests/synthesis.motif "70-109 (len 40)" \
    "synthesised budgets are used by the search"
check tests/synth_errors.motif   reject  "orphan, shared and zero-rate contracts rejected"
check tests/contract.motif       reject  "fixed budget that breaks its 'require' clause rejected"
check tests/infeasible.motif     reject  "no budget can meet the target: rejected with a reason"


echo
echo "== Paper extensions: position-restricted errors, learning from known sites =="
check_output_contains tests/positions.motif "4-9 (len 6)  TATAAA" "'at 1,2,5': variant at position 1 found"
check_output_contains tests/positions.motif "31-36 (len 6)  AATATA" "'at 1,2,5': variant at position 5 found"
if (cd tests && ../build/motifql --symbols positions.motif | grep -q "AAGAAA"); then
    printf "  FAIL  variant at position 3 matched although 3 is not listed\n"; FAIL=$((FAIL + 1))
else
    printf "  PASS  'at 1,2,5': variant at position 3 is not matched\n"; PASS=$((PASS + 1))
fi
check tests/positions_errors.motif reject "'at' on a non-linear motif, out-of-range position, 'at ?' without '?' budget"
if (cd tests && ../build/motifql --symbols --sites learn_sites.fa --empirical learn_bg.fa learn.motif 2>/dev/null \
      | grep -qF "note: chose pas = 1 mismatch at 2"); then
    printf "  PASS  learning from known sites: compiler finds the one position where they vary\n"; PASS=$((PASS + 1))
else
    printf "  FAIL  learning from known sites did not choose '1 mismatch at 2'\n"; FAIL=$((FAIL + 1))
fi

echo
echo "passed: $PASS   failed: $FAIL"
[ $FAIL -eq 0 ] || exit 1
