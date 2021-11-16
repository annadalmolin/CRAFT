#!/bin/bash

## Script detecting ORFs in input sequences

## input:
# circRNA sequences file: "/data/functional_predictions/backsplice_sequence_1.fa"
# ORFfinder parameters: $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
#	GEN_CODE: genetic code to use (1-31), see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details; default: 1
#	START_CODON: start codon to use (0 = "ATG" only, 1 = "ATG" and alternative initiation codons, 2 = any sense codon); default: 0
#	MIN_LENGTH: minimal length of the ORF (nt); allowed values are: 30, 75, 150. Default: 30
#	NESTED_ORFS: ignore nested ORFs (ORFs completely placed within another); allowed values: "TRUE", "FALSE". Default: "FALSE"
#	STRAND: output ORFs on specified strand only; allowed values: "both", "plus", "minus". Default: "plus"

## command: /scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND



GEN_CODE=$1
START_CODON=$2
MIN_LENGTH=$3
NESTED_ORFS=$4
STRAND=$5


mkdir /data/functional_predictions/ORF_detection
cd /data/functional_predictions/ORF_detection

cat ../backsplice_sequence_1.fa | sed -r 's/(^[ATCG]+)/\1\1/; s/:/_/' > backsplice_sequence_per_ORF_$MIN_LENGTH.fa

cp ../backsplice_sequence_1.txt backsplice_sequence_$MIN_LENGTH.txt


# create file with backsplice position
cat ../backsplice_circRNA_length_1.txt | grep -v "circ_id" > circ_id_length.txt
cat circ_id_length.txt | awk '{print $1,$2-1,$2+1}' | sed -e 's/ /\t/g' | sort -k1,1n -k2,2n > backsplice_position.bed


# ORFfinder

mkdir ORFfinder
cd ORFfinder


N=$( cat ../circ_id_length.txt | wc -l )


if [ $N -gt 0 ]
then

	/tools/ORFfinder -in ../backsplice_sequence_per_ORF_$MIN_LENGTH.fa -g $GEN_CODE -s $START_CODON -ml $MIN_LENGTH -n $NESTED_ORFS -strand $STRAND -out result_list_ORF_$MIN_LENGTH.txt -outfmt 0 -logfile ORF0_$MIN_LENGTH.log
	/tools/ORFfinder -in ../backsplice_sequence_per_ORF_$MIN_LENGTH.fa -g $GEN_CODE -s $START_CODON -ml $MIN_LENGTH -n $NESTED_ORFS -strand $STRAND -out result_list_CDS_$MIN_LENGTH.txt -outfmt 1 -logfile ORF1_$MIN_LENGTH.log
	/tools/ORFfinder -in ../backsplice_sequence_per_ORF_$MIN_LENGTH.fa -g $GEN_CODE -s $START_CODON -ml $MIN_LENGTH -n $NESTED_ORFS -strand $STRAND -out result_text_ORF_$MIN_LENGTH.txt -outfmt 2 -logfile ORF2_$MIN_LENGTH.log
	/tools/ORFfinder -in ../backsplice_sequence_per_ORF_$MIN_LENGTH.fa -g $GEN_CODE -s $START_CODON -ml $MIN_LENGTH -n $NESTED_ORFS -strand $STRAND -out result_table_ORF_$MIN_LENGTH.txt -outfmt 3 -logfile ORF3_$MIN_LENGTH.log

	cat result_list_ORF_$MIN_LENGTH.txt | grep ">" | sed -e 's/>lcl|ORF[0-9]*_//; s/:/\t/g; s/ /\t/; s/ /_/g; s/_/:/' > result_$MIN_LENGTH.txt			#start-1, end-1
	cat result_$MIN_LENGTH.txt | grep -v "partial" | awk '{print $1, $2, $3+1}' | sed -e 's/ /\t/g' | sort -k1,1n -k2,2n > result_$MIN_LENGTH.bed			#start-1, end


	# detect ORFs crossing the backsplice junction
	bedtools intersect -a ../backsplice_position.bed -b result_$MIN_LENGTH.bed -wa -wb > ORF_backsplice_$MIN_LENGTH.bed

	echo "circ_id ORF start end length_nt length_aa" | sed -e 's/ /\t/g' > ORF_backsplice2_$MIN_LENGTH.bed
	cat ORF_backsplice_$MIN_LENGTH.bed | awk '{print $4, '$MIN_LENGTH', $5, $6, $6-$5, ($6-$5)/3}' | sed -e 's/ /\t/g' >> ORF_backsplice2_$MIN_LENGTH.bed		#start-1, end

	cat result_$MIN_LENGTH.txt > result.txt																#start-1, end-1

	cat ORF_backsplice2_$MIN_LENGTH.bed > ORF_backsplice2.bed													#start-1, end

	rm result_$MIN_LENGTH.txt
	rm result_$MIN_LENGTH.bed
	rm ORF_backsplice_$MIN_LENGTH.bed
	rm ORF_backsplice2_$MIN_LENGTH.bed

fi


cat ORF_backsplice2.bed | sort -k1,1n | uniq > ORF_backsplice0.txt
cat ORF_backsplice0.txt | grep "circ_id" > ORF_backsplice.txt
cat ORF_backsplice0.txt | grep -v "circ_id" | awk '{print $1, $2, $3+1, $4, $5, $6}' | sed -e 's/ /\t/g' >> ORF_backsplice.txt						#start, end

rm ORF_backsplice2.bed



# rolling ORFs
echo "circ_id start end open" | sed -e 's/ /\t/g' > result_open.txt
cat result.txt | grep "partial" | sed -e 's/unnamed_protein_product,_partial/open/' | awk '{print $1, $2+1, $3+1, $4}' | sed -e 's/ /\t/g' >> result_open.txt

cat result.txt | grep "partial" | awk '{print $1, $2, $3+1}' | sed -e 's/ /\t/g' | sort -k1,1n -k2,2n > result_open.bed
# detect ORFs crossing the backsplice junction
bedtools intersect -a ../backsplice_position.bed -b result_open.bed -wa -wb | awk '{print $4, $5, $6}' | sed -e 's/ /\t/g'> ORF_backsplice_open0.txt

echo "circ_id start end open" | sed -e 's/ /\t/g' > ORF_backsplice_open.txt
cat ORF_backsplice_open0.txt | awk '{print $1, $2+1, $3, "open"}' | sed -e 's/ /\t/g' >> ORF_backsplice_open.txt

rm result.txt
rm result_open.txt
rm result_open.bed



# retrieve ORF sequence (nt)
cat result_list_CDS_$MIN_LENGTH.txt | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | sed '/^$/d' > result_list_CDS.fa
cat result_list_CDS.fa | grep ">" | sed -e 's/>lcl|//; s/ ORF.*$//; s/:/ /; s/-/_/; s/-/\t/; s/_/:/; s/_/-/' | awk '{print $1, $2, $3}' | sed -e 's/ /\t/g' | sed -e 's/\t/_/g' > head.txt
cat result_list_CDS.fa | grep -v ">" > seq.txt
paste head.txt seq.txt > result_list_CDS.txt

rm head.txt
rm seq.txt



# retrieve ORF sequence (aa)
cat result_list_ORF_$MIN_LENGTH.txt | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | sed '/^$/d' > result_list_ORF.fa
cat result_list_ORF.fa | grep ">" | sed -e 's/>lcl|ORF[0-9]*_//; s/ unnamed.*$//; s/:/\t/g; s/_/:/' | awk '{print $1, $2, $3}' | sed -e 's/ /\t/g' | sed -e 's/\t/_/g' > head.txt
cat result_list_ORF.fa | grep -v ">" > seq.txt
paste head.txt seq.txt > result_list_ORF.txt

rm head.txt
rm seq.txt



cd ..


rm circ_id_length.txt


cd ..

