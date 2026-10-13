#!/bin/bash
# Compile the MotifQL compiler to WebAssembly (WASI) for the browser playground.
# Needs wasi-sdk: https://github.com/WebAssembly/wasi-sdk/releases (tested with wasi-sdk 24).
set -e
WASI_SDK=${WASI_SDK:-/opt/wasi-sdk}
cd "$(dirname "$0")/.."
"$WASI_SDK/bin/clang" --target=wasm32-wasip1 --sysroot="$WASI_SDK/share/wasi-sysroot" \
    -std=c11 -O2 -Wall -Wextra -o web/motifql.wasm src/*.c -lm
echo "wrote web/motifql.wasm"
