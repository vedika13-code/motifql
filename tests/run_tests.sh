#!/usr/bin/env bash
# Phase 1 front-end test suite.
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
echo "passed: $PASS   failed: $FAIL"
[ $FAIL -eq 0 ] || exit 1
