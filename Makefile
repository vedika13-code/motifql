CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -g -O2
LDLIBS = -lm
SRC = src/alphabet.c src/lexer.c src/ast.c src/symtab.c src/diagnostics.c src/parser.c src/sema.c \
      src/nfa.c src/dfa.c src/codegen.c src/fasta.c src/matcher.c src/specificity.c src/sensitivity.c src/synthesis.c src/main.c
BIN = build/motifql

all: $(BIN)

$(BIN): $(SRC) src/*.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRC) -o $(BIN) $(LDLIBS)

test: $(BIN)
	@bash tests/run_tests.sh

clean:
	rm -rf build

.PHONY: all test demo clean

demo: $(BIN)
	@echo "--- tokens ---"; ./$(BIN) --tokens tests/promoter.motif | head -20
	@echo "--- ast ---"; ./$(BIN) --ast tests/promoter.motif
	@echo "--- symbols ---"; ./$(BIN) --symbols tests/promoter.motif
	@echo "--- nfa (polyA, shows the mismatch-tolerance layering) ---"; ./$(BIN) --nfa tests/promoter.motif | sed -n '/motif polyA/,/^$$/p'
	@echo "--- dfa ---"; cd tests && ../$(BIN) --dfa promoter.motif
	@echo "--- bytecode + search ---"; cd tests && ../$(BIN) --bytecode --symbols promoter.motif
