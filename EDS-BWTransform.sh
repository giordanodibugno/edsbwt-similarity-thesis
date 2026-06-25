#!/bin/bash

GSUFPATH="gsufsort"

NAMEFILE=$1
OUTPUT=$2

./eds_to_fasta $NAMEFILE.eds $OUTPUT

$GSUFPATH/"gsufsort" $OUTPUT.fasta --da --bwt --output $OUTPUT
rm $OUTPUT.fasta
./da_to_everything $OUTPUT

echo "File "$NAMEFILE" done."
