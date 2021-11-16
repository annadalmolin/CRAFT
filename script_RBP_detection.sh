#!/bin/bash

## Script detecting RBP binding sites in input sequences

## input:
# circRNA sequences file: "../backsplice_sequence_1.fa"
# beRBP parameters: $PWM, $RBP
#	PWM: position weight matrix
#	RBP: RBP name
# prefix for genome file: "$PREFIX_UCSC"

## command: /scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC



PWM=$1
RBP=$2
PREFIX_UCSC=$3


mkdir /data/functional_predictions/RBP_detection
cd /data/functional_predictions/RBP_detection

cat ../backsplice_sequence_1.fa | grep ">" > headers.txt

cat ../backsplice_sequence_1.fa | grep -v ">" | cut -c1-20 > first_20_characters.txt
cat ../backsplice_sequence_1.fa | grep -v ">" > sequences.txt

paste sequences.txt first_20_characters.txt | sed -e 's/\t//' > sequences_plus_first_20_characters.txt

paste headers.txt sequences_plus_first_20_characters.txt | sed -e 's/\t/\n/' > backsplice_sequence_per_RBP.fa



# beRBP

cp backsplice_sequence_per_RBP.fa /tools/beRBP/work/temp/analysis_RBP.fasta

mkdir beRBP
cd beRBP

cd /tools/beRBP/work
../code/general_sPWM.sh analysis_RBP $PWM $RBP $PREFIX_UCSC &>temp/analysis_RBP.log

mv analysis_RBP/ /data/functional_predictions/RBP_detection/beRBP
mv temp/analysis_RBP.fasta /data/functional_predictions/RBP_detection/beRBP
mv temp/analysis_RBP.log /data/functional_predictions/RBP_detection/beRBP

cd /data/functional_predictions/RBP_detection/beRBP
cd analysis_RBP

cat resultMatrix.tsv | head -n1 > resultMatrix_a.tsv

for CIRC in $( cat ../../../circRNA_length.txt | cut -f1 )
do
	L=$( cat ../../../circRNA_length.txt | grep $CIRC | cut -f2 )
	cat resultMatrix.tsv | grep $CIRC | awk '$5 <= '$L' {print}' >> resultMatrix_a.tsv
done


cat resultMatrix_a.tsv | sed -e 's/seqID/circ_id/' > resultMatrix_b.txt

cat resultMatrix_a.tsv | cut -f1-6 > c1-6.txt
cat resultMatrix_a.tsv | cut -f7 | sed -e 's/\./,/' > c7.txt
paste c1-6.txt c7.txt | sed -e 's/ /\t/g' | sed -e 's/seqID/circ_id/' > resultMatrix_b.tsv


#select suggested binding threshold from beRBP ("voteFrac" > 0.35)
#head -n1 resultMatrix_b.tsv > resultMatrix_bind.tsv
#cat resultMatrix_b.tsv | awk '$4==1 {print}' | sort -k7nr >> resultMatrix_bind.tsv

rm resultMatrix_a.tsv
rm c1-6.txt
rm c7.txt

cd ..
cd ..


rm headers.txt
rm first_20_characters.txt
rm sequences.txt
rm sequences_plus_first_20_characters.txt

cd ..

