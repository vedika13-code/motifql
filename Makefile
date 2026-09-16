CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -g -O0
SRC = src/alphabet.c src/lexer.c src/ast.c src/symtab.c src/diagnostics.c src/parser.c src/sema.c src/main.c
BIN = build/motifql

all: $(BIN)

$(BIN): $(SRC) src/*.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRC) -o $(BIN)

test: $(BIN)
	@bash tests/run_tests.sh

demo: $(BIN)
	@echo "--- tokens ---"; ./$(BIN) --tokens tests/promoter.motif | head -20
	@echo "--- ast ---"; ./$(BIN) --ast tests/promoter.motif
	@echo "--- symbols ---"; ./$(BIN) --symbols tests/promoter.motif

clean:
	rm -rf build

.PHONY: all test demo clean
