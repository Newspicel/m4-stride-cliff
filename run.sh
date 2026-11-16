#!/bin/sh
# Builds the reproducer and runs every combination twice.
set -e
cd "$(dirname "$0")"
clang -O2 -o cliff cliff.c
echo "$(sysctl -n machdep.cpu.brand_string), macOS $(sw_vers -productVersion)"
echo "columns: pass 0 (first pass) to pass 4 (last pass)"
for round in 1 2; do
  for skew in 0 128; do
    for v in A B; do
      ./cliff $v --skew $skew
      ./cliff $v --skew $skew --dit
    done
  done
done
