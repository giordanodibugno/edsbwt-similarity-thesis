#!/bin/bash

BCR=1

BCRPATH="BCR_LCP_GSA"
GSUFPATH="gsufsort"

NAMEFILE=$1
OUTPUT=$2

./eds_to_fasta $NAMEFILE.eds $OUTPUT

if [ $BCR -eq 1 ]
then
    $BCRPATH/"BCR_LCP_GSA" $OUTPUT.fasta $OUTPUT
	rm $OUTPUT.fasta
	rm $OUTPUT.len
	rm $OUTPUT.info
	./EOFpos_to_everything $OUTPUT
else
    $GSUFPATH/"gsufsort" $OUTPUT.fasta --da --bwt --output $OUTPUT
	rm $OUTPUT.fasta
	./da_to_everything $OUTPUT
fi 

echo "File "$NAMEFILE" done."
