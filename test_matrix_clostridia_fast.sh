#!/bin/bash
set -e

BASE="$HOME/Desktop/Programma/bacteria/clostridia/matrices"
SIM="$HOME/Desktop/Programma/EDS-BWT/EDSBWTsimilarity"
OUT="$BASE/results"

mkdir -p "$OUT"

run_one () {
    local DIR="$1"
    local NAME="$2"
    local CSV="$OUT/${NAME}_matrix.csv"

    local A="$DIR/C_botulinum.eds"
    local B="$DIR/C_difficile.eds"
    local C="$DIR/C_perfringens.eds"

    echo "Calcolo $NAME..."

    echo "  C_botulinum vs C_difficile"
    AB=$("$SIM" "$A" "$B" | tail -n 1 | grep -Eo '[0-9]+(\.[0-9]+)?' | tail -n 1)

    echo "  C_botulinum vs C_perfringens"
    AC=$("$SIM" "$A" "$C" | tail -n 1 | grep -Eo '[0-9]+(\.[0-9]+)?' | tail -n 1)

    echo "  C_difficile vs C_perfringens"
    BC=$("$SIM" "$B" "$C" | tail -n 1 | grep -Eo '[0-9]+(\.[0-9]+)?' | tail -n 1)

    cat > "$CSV" <<EOF
eds,C_botulinum,C_difficile,C_perfringens
C_botulinum,1,$AB,$AC
C_difficile,$AB,1,$BC
C_perfringens,$AC,$BC,1
EOF

    echo "Creato: $CSV"
}

run_one "$BASE/junctions" "junctions"
run_one "$BASE/mincard" "mincard"
