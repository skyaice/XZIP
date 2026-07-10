#!/bin/sh
set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
binary="$root_dir/build/xzip"

test -x "$binary"
output=$($binary --help 2>&1)
printf '%s\n' "$output" | grep -q "xzip compress"
printf '%s\n' "$output" | grep -q "xzip decompress"

echo "CLI smoke test passed."
