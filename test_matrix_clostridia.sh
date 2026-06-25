#!/bin/bash

set -e

BASE="$HOME/Desktop/Programma/bacteria/clostridia/matrices"
SIM="$HOME/Desktop/Programma/EDS-BWT/EDSBWTsimilarity"
OUT="$BASE/results"

mkdir -p "$OUT"

make_matrix () {
    local DIR="$1"
    local NAME="$2"
    local CSV="$OUT/${NAME}_matrix.csv"

    mapfile -t FILES < <(find "$DIR" -maxdepth 1 -name "*.eds" | sort)

    echo -n "eds" > "$CSV"
    for f in "${FILES[@]}"; do
        b=$(basename "$f" .eds)
        echo -n ",$b" >> "$CSV"
    done
    echo >> "$CSV"

    for f1 in "${FILES[@]}"; do
        b1=$(basename "$f1" .eds)
        echo -n "$b1" >> "$CSV"

        for f2 in "${FILES[@]}"; do
            sim=$("$SIM" "$f1" "$f2" | tail -n 1 | grep -Eo '[0-9]+(\.[0-9]+)?' | tail -n 1)
            echo -n ",$sim" >> "$CSV"
        done

        echo >> "$CSV"
    done

    echo "Creato: $CSV"
}

make_matrix "$BASE/junctions" "junctions"
make_matrix "$BASE/mincard" "mincard"
