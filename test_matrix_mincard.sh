#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MINCARD_DIR="$SCRIPT_DIR/../dataset/eds2mincard"

shopt -s nullglob
files=("$MINCARD_DIR"/*.eds)

if [[ ${#files[@]} -eq 0 ]]; then
    echo "Errore: nessun file .eds trovato in $MINCARD_DIR" >&2
    exit 1
fi

echo "Dataset mincard: $MINCARD_DIR"
echo "File EDS: ${#files[@]}"

OUTPUT_PREFIX="mincard_similarity" \
    "$SCRIPT_DIR/test_matrix.sh" "${files[@]}"
