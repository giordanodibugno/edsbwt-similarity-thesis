#!/bin/bash

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATASET_DIR="$SCRIPT_DIR/../dataset"
INDEX_DIR="$SCRIPT_DIR/edsbwt_form"
OUT_DIR="$SCRIPT_DIR/results"

SEARCH_BIN="$SCRIPT_DIR/EDSBWTsearch"
EDS2FASTA="$SCRIPT_DIR/eds_to_fasta"
DA_TOOL="$SCRIPT_DIR/da_to_everything"
GSUF_BIN="$SCRIPT_DIR/gsufsort/gsufsort"

mkdir -p "$OUT_DIR" "$INDEX_DIR"

check_executable() {
    local tool="$1"
    if [[ ! -x "$tool" ]]; then
        echo "Errore: tool richiesto non trovato o non eseguibile: $tool" >&2
        echo "Ricompila con: make -B similarity RECOVERBW=0 SDSL_INC=... SDSL_LIB=..." >&2
        exit 1
    fi
}

cleanup_partial_index() {
    local base="$1"
    rm -f "$base.fasta" "$base.len" "$base.info" "$base.empty.info"
    rm -f "$base.bwt" "$base.4.da"
    rm -f "$base.ebwt" "$base.lcp" "$base.da" "$base.posSA" "$base.SAP" "$base.bitvector"
    rm -f "${base}_info.aux" "${base}_alpha.txt" "${base}_tableOcc.txt"
    rm -f "${base}"_bwt_*.aux "${base}"_bv_*.aux
}

index_is_ready() {
    local eds="$1"
    local base="$2"

    [[ -f "$base.ebwt" && -f "${base}_info.aux" && "$base.ebwt" -nt "$eds" && "${base}_info.aux" -nt "$eds" ]]
}

build_index() {
    local eds="$1"
    local name="$2"
    local base="$INDEX_DIR/$name"

    if index_is_ready "$eds" "$base"; then
        echo "Indice gia presente: $name" >&2
        return 0
    fi

    echo "Creo indice: $name" >&2
    cleanup_partial_index "$base"
    "$EDS2FASTA" "$eds" "$base" > /dev/null 2>&1

    if [[ ! -f "$base.fasta" ]]; then
        echo "Errore: FASTA non creato per $eds" >&2
        cleanup_partial_index "$base"
        return 1
    fi

    if ! "$GSUF_BIN" "$base.fasta" --da --bwt --output "$base" > /dev/null 2>&1; then
        echo "Errore durante gsufsort per $eds" >&2
        cleanup_partial_index "$base"
        return 1
    fi

    if [[ ! -f "$base.bwt" || ! -f "$base.4.da" ]]; then
        echo "Errore: output gsufsort incompleto per $eds" >&2
        cleanup_partial_index "$base"
        return 1
    fi

    rm -f "$base.fasta" "$base.len" "$base.info"

    "$DA_TOOL" "$base" > /dev/null 2>&1

    rm -f "$base.bwt" "$base.4.da" "$base.empty.info"

    if [[ ! -f "${base}_info.aux" ]]; then
        echo "Errore: file ${base}_info.aux non creato" >&2
        cleanup_partial_index "$base"
        return 1
    fi
}

directional_similarity() {
    local index_base="$1"
    local query_eds="$2"

    "$SEARCH_BIN" "$index_base" "$query_eds" 2>&1 | awk -F: "/^Similarita totale EDS/ {gsub(/[[:space:]]/, \"\", \$2); value=\$2} END {if (value != \"\") print value}"
}

check_executable "$SEARCH_BIN"
check_executable "$EDS2FASTA"
check_executable "$DA_TOOL"
check_executable "$GSUF_BIN"

if [[ $# -gt 0 ]]; then
    files=("$@")
else
    files=("$DATASET_DIR"/*.eds)
fi

if [[ ${#files[@]} -eq 0 || ! -e "${files[0]}" ]]; then
    echo "Errore: nessun file .eds trovato in $DATASET_DIR" >&2
    exit 1
fi

names=()
indices=()
for f in "${files[@]}"; do
    name="$(basename "$f" .eds)"
    names+=("$name")
    indices+=("$INDEX_DIR/$name")
done

timestamp="$(date +%Y%m%d_%H%M%S)"
output_prefix="${OUTPUT_PREFIX:-similarity}"
out_csv="$OUT_DIR/${output_prefix}_matrix_$timestamp.csv"
out_long="$OUT_DIR/${output_prefix}_pairs_$timestamp.tsv"

echo "Matrice di similarita ottimizzata"
echo "File EDS: ${#files[@]}"
echo "Cartella indici: $INDEX_DIR"
echo "Output matrice: $out_csv"
echo "Output coppie:  $out_long"
echo ""

echo "=== Costruzione/riuso indici EDS-BWT ===" >&2
for i in "${!files[@]}"; do
    build_index "${files[$i]}" "${names[$i]}" || exit 1
done

printf "P1\tP2\tA(P1,P2)\tA(P2,P1)\tSimilarity\n" > "$out_long"
declare -A matrix

n=${#files[@]}

echo "" >&2
echo "=== Calcolo matrice ===" >&2
for ((i = 0; i < n; i++)); do
    for ((j = i; j < n; j++)); do
        b1="$(basename "${files[$i]}")"
        b2="$(basename "${files[$j]}")"
        echo "Confronto: $b1 vs $b2" >&2

        if [[ $i -eq $j ]]; then
            a12="1"
            a21="1"
            sim="1"
        else
            a12="$(directional_similarity "${indices[$j]}" "${files[$i]}")"
            a21="$(directional_similarity "${indices[$i]}" "${files[$j]}")"

            if [[ -z "$a12" || -z "$a21" ]]; then
                a12="ERROR"
                a21="ERROR"
                sim="ERROR"
                echo "  Errore nel confronto $b1 vs $b2" >&2
            else
                sim="$(awk -v x="$a12" -v y="$a21" "BEGIN {printf \"%.6g\", (x + y) / 2}")"
                echo "  Similarity: $sim" >&2
            fi
        fi

        matrix["$i,$j"]="$sim"
        matrix["$j,$i"]="$sim"
        printf "%s\t%s\t%s\t%s\t%s\n" "$b1" "$b2" "$a12" "$a21" "$sim" >> "$out_long"
    done
done

printf "eds" > "$out_csv"
for name in "${names[@]}"; do
    printf ",%s" "$name" >> "$out_csv"
done
printf "\n" >> "$out_csv"

for ((i = 0; i < n; i++)); do
    printf "%s" "${names[$i]}" >> "$out_csv"
    for ((j = 0; j < n; j++)); do
        printf ",%s" "${matrix[$i,$j]}" >> "$out_csv"
    done
    printf "\n" >> "$out_csv"
done

echo ""
echo "Finito."
echo "Matrice salvata in: $out_csv"
echo "Dettaglio coppie salvato in: $out_long"
