#!/bin/bash

## Script detecting miRNA binding sites in input sequences

## input:
# circRNA sequences file: "/data/backsplice_sequence_1.fa"
# miRNA database: "/input/mature_$SPECIES.fa"
# miRanda parameters: $SCORE, $ENERGY
#	SCORE: best predictions are obtained with higher score
#	ENERGY: best predictions are obtained with lower energy
# AGO2 binding file (optional): "/data/input/AGO2_binding_sites.bed"
# file of region positions for each circRNA (optional): "/data/sequence_extraction/region_to_extract_1.bed"

## command: /scripts/script_miRNA_detection.sh /input/mature_$SPECIES.fa $SCORE $ENERGY



MIRNA=$1
AGO2="/data/input/AGO2_binding_sites.bed"
SEQUENCE="/data/sequence_extraction/region_to_extract_1.bed"


mkdir /data/functional_predictions/miRNA_detection
cd /data/functional_predictions/miRNA_detection

cat /data/sequence_extraction/backsplice_sequence_1.fa | grep ">" > headers.txt

cat /data/sequence_extraction/backsplice_sequence_1.fa | grep -v ">" | cut -c1-20 > first_20_characters.txt
cat /data/sequence_extraction/backsplice_sequence_1.fa | grep -v ">" > sequences.txt

paste sequences.txt first_20_characters.txt | sed -e 's/\t//' > sequences_plus_first_20_characters.txt

paste headers.txt sequences_plus_first_20_characters.txt | sed -e 's/\t/\n/' > backsplice_sequence_per_miRNA.fa



# miRanda

mkdir miRanda
cd miRanda

SCORE=$2
ENERGY=$3

/tools/miRanda/src/miranda $MIRNA ../backsplice_sequence_per_miRNA.fa -sc $SCORE -en $ENERGY -quiet -out output_miRanda.txt

echo "#circ_id miRNA_id score energy circ_start circ_end miRNA_start miRNA_end align_length circ_align_perc miRNA_align_perc" | sed -e 's/ /\t/g' > output_miRanda_b.txt
cat output_miRanda.txt | grep ">" | grep -v ">>" | sed -e 's/>//' |  awk '{print $2,$1,$3,$4,$7,$8,$5,$6,$9,$11,$10}' | sed -e 's/ /\t/g' >> output_miRanda_b.txt
cat output_miRanda_b.txt | sed -e 's/ +/\t/g; s/#//' > output_miRanda_b_per_R.txt


cat output_miRanda_b_per_R.txt | head -n1 > output_miRanda_c_per_R.txt

for CIRC in $( cat ../../circRNA_length.txt | cut -f1 )
do
	L=$( cat ../../circRNA_length.txt | grep $CIRC | cut -f2 )
	cat output_miRanda_b_per_R.txt | grep $CIRC | awk '$5 <= '$L' {print}' >> output_miRanda_c_per_R.txt
done

#cat ../../backsplice_sequence_1.fa | grep -v ">" | awk '{print length}' > circRNA_length.txt
#paste ../headers.txt circRNA_length.txt | sed -e 's/>//' > backsplice_circRNA_length.txt

rm output_miRanda_b.txt
rm output_miRanda_b_per_R.txt

cd ..



# PITA

mkdir PITA
cd PITA

/tools/pita_prediction/pita_prediction.pl -utr ../backsplice_sequence_per_miRNA.fa -mir $MIRNA -prefix pred -gxp 2>>pita.err | tee -a pita.log

echo "circ_id miRNA_id circ_start circ_end seed loop dGduplex dG5 dG3 dG0 dG1 dGopen ddG" | sed -e 's/ /\t/g' > pred_pita_results_b.txt
cat pred_pita_results.tab | grep -v "UTR" | awk '{print $1,$2,$4,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | sed -e 's/ /\t/g' >> pred_pita_results_b.txt

echo "circ_id miRNA num_sites score" | sed -e 's/ /\t/g' > pred_pita_results_targets_b.txt
cat pred_pita_results_targets.tab | grep -v "RefSeq" | awk '{print $1,$2,$3,$4}' | sed -e 's/ /\t/g' >> pred_pita_results_targets_b.txt

cat pred_pita_results_b.txt | head -n1 > pred_pita_results_c.txt

for CIRC in $( cat ../../circRNA_length.txt | cut -f1 )
do
	L=$( cat ../../circRNA_length.txt | grep $CIRC | cut -f2 )
	cat pred_pita_results_b.txt | grep $CIRC | awk '$3 <= '$L' {print}' >> pred_pita_results_c.txt
done

rm tmp_seqfile1
rm tmp_seqfile2
rm pred_pita_results_b.txt

cd ..


rm headers.txt
rm first_20_characters.txt
rm sequences.txt
rm sequences_plus_first_20_characters.txt



# overlap with AGO2 binging sites

if [ -e $AGO2 ]
then

	if [ -e genomic_sequence_position_conversion.txt ]
	then
		rm genomic_sequence_position_conversion.txt
	fi

	#1# mapping genomic position-sequence position

	for CIRC in $( cat /data/list_backsplice.txt | cut -f1 )
	do
		CIRC_BED=$( echo $CIRC | sed -e 's/:/\t/; s/-/\t/' | awk '{print $1, $2-1, $3}' | sed -e 's/ /:/; s/ /-/' )

		grep $CIRC_BED $SEQUENCE > region.bed
		N=$( cat region.bed | wc -l )
		D=0

		if [ -e genomic_position.txt ]
		then
			rm genomic_position.txt
		fi

		if [ -e sequence_position.txt ]
		then
			rm sequence_position.txt
		fi

		for ((i=1; i<=$N; i++))
		do
		        A=$( cat region.bed | sed -n $i' p' | cut -f 2 )
		        A=$(( $A+1 ))
		        B=$( cat region.bed | sed -n $i' p' | cut -f 3 )
		        seq $A 1 $B >> genomic_position.txt

		        C=$(( $D+1 ))
		        D=$(( $D + $B-$A+1 ))
		        seq $C 1 $D >> sequence_position.txt
		done

		M=$( cat genomic_position.txt | wc -l )

		if [ -e circ.txt ]
		then
			rm circ.txt
		fi

		for ((i=1; i<=$M; i++))
		do
			echo $CIRC >> circ.txt
		done

		paste circ.txt genomic_position.txt sequence_position.txt >> genomic_sequence_position_conversion.txt
	done

	rm region.bed
	rm circ.txt
	rm genomic_position.txt
	rm sequence_position.txt



	# miRanda

	#2# intersection of the two files

	cd miRanda

	FILE_PRED="output_miRanda_c_per_R.txt"
	N=$( cat $FILE_PRED | wc -l )
	echo "gen_start" "gen_end" "id" | sed -e 's/ /\t/' > genomic.txt

	for ((i=2; i<=$N; i++))
	do
		CIRC=$( cat $FILE_PRED | sed -n $i' p' | cut -f 1 )
		SEQ_START=$( cat $FILE_PRED | sed -n $i' p' | cut -f 5 )
		SEQ_END=$( cat $FILE_PRED | sed -n $i' p' | cut -f 6 )

		cat ../genomic_sequence_position_conversion.txt | grep $CIRC > circ_conversion.txt
		GEN_START=$( cat circ_conversion.txt | awk '$3=='$SEQ_START' {print $2}')
		GEN_END=$( cat circ_conversion.txt | awk '$3=='$SEQ_END' {print $2}')
		echo $GEN_START $GEN_END "id"$(($i-1)) | sed -e 's/ /\t/g' >> genomic.txt

		CHR=$( echo $CIRC | sed -e 's/:.*//g' )
		STRAND=$( grep $CIRC /data/list_backsplice.txt | cut -f 2 )
		echo $CHR $(($GEN_START-1)) $GEN_END "id"$(($i-1)) "." $STRAND | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n >> output_miRanda_b_gen.bed
	done

	cat genomic.txt | awk '$3 != "" {print}' > genomic_1.txt
	cat genomic.txt | awk '$3 == "" {print $1, $1+22, $2}' | sed -e 's/ /\t/g' > genomic_2.txt
	cat genomic_1.txt genomic_2.txt | sed -e 's/id//' | sort -k3,3n | awk '{print $1, $2, "id"$3}' | sed -e 's/ /\t/g' > genomic.txt

	cat output_miRanda_b_gen.bed | awk '($2 != -1) && ($4 != ".") {print}' > output_miRanda_b_gen_1.bed
	cat output_miRanda_b_gen.bed | awk '($2 != -1) && ($4 == ".") {print $1, $2, $2+22, $3, $4, $5}' | sed -e 's/ /\t/g' > output_miRanda_b_gen_2.bed
	cat output_miRanda_b_gen_1.bed output_miRanda_b_gen_2.bed | sed -e 's/id//' | sort -k4,4n | awk '{print $1, $2, $3, "id"$4, $5, $6}' | sed -e 's/ /\t/g' > output_miRanda_c_gen.bed


	paste $FILE_PRED genomic.txt > output_miRanda_d_per_R.txt

	bedtools intersect -a $AGO2 -b output_miRanda_c_gen.bed -wo > output_miRanda_AG02.bed.txt

	cat output_miRanda_AG02.bed.txt | awk '$13>=5 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' | sed -e 's/ /\t/g' > output_miRanda_AG02.bed


	#3# select only AGO2 overlapping sites from the original file

	cat output_miRanda_AG02.bed | cut -f10 | sort | uniq > id_AGO2.txt

	cat output_miRanda_c_per_R.txt | head -n1 > output_miRanda_per_R.txt
	grep -f id_AGO2.txt output_miRanda_d_per_R.txt | cut -f 1-11 >> output_miRanda_per_R.txt

	echo "MiRanda miRNA binding sites overlapped by AGO2:" $(( $( cat output_miRanda_per_R.txt | wc -l ) -1))


	rm genomic.txt
	rm genomic_1.txt
	rm genomic_2.txt
	rm circ_conversion.txt
	rm output_miRanda_b_gen.bed
	rm output_miRanda_b_gen_1.bed
	rm output_miRanda_b_gen_2.bed
	rm output_miRanda_c_gen.bed
	rm output_miRanda_AG02.bed.txt
	rm output_miRanda_AG02.bed
	rm id_AGO2.txt
	rm output_miRanda_d_per_R.txt

	cd ..



	# PITA

	#2# intersection of the two files

	cd PITA

	FILE_PRED="pred_pita_results_c.txt"
	N=$( cat $FILE_PRED | wc -l )
	echo "gen_start" "gen_end" "id" | sed -e 's/ /\t/' > genomic.txt

	for ((i=2; i<=$N; i++))
	do
		CIRC=$( cat $FILE_PRED | sed -n $i' p' | cut -f 1 )
		SEQ_START=$( cat $FILE_PRED | sed -n $i' p' | cut -f 3 )
		SEQ_END=$( cat $FILE_PRED | sed -n $i' p' | cut -f 4 )

		cat ../genomic_sequence_position_conversion.txt | grep $CIRC > circ_conversion.txt
		GEN_START=$( cat circ_conversion.txt | awk '$3=='$SEQ_START' {print $2}')
		GEN_END=$( cat circ_conversion.txt | awk '$3=='$SEQ_END' {print $2}')
		echo $GEN_START $GEN_END "id"$(($i-1)) | sed -e 's/ /\t/' >> genomic.txt

		CHR=$( echo $CIRC | sed -e 's/:.*//g' )
		STRAND=$( grep $CIRC /data/list_backsplice.txt | cut -f 2)
		echo $CHR $(($GEN_START-1)) $GEN_END "id"$(($i-1)) "." $STRAND | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n >> pred_pita_results_b_gen.bed
	done

	cat genomic.txt | awk '$3 != "" {print}' > genomic_1.txt
	cat genomic.txt | awk '$3 == "" {print $1, $1+22, $2}' | sed -e 's/ /\t/g' > genomic_2.txt
	cat genomic_1.txt genomic_2.txt | sed -e 's/id//' | sort -k3,3n | awk '{print $1, $2, "id"$3}' | sed -e 's/ /\t/g' > genomic.txt

	cat pred_pita_results_b_gen.bed | awk '($2 != -1) && ($4 != ".") {print}' > pred_pita_results_b_gen_1.bed
	cat pred_pita_results_b_gen.bed | awk '($2 != -1) && ($4 == ".") {print $1, $2, $2+22, $3, $4, $5}' | sed -e 's/ /\t/g' > pred_pita_results_b_gen_2.bed
	cat pred_pita_results_b_gen_1.bed pred_pita_results_b_gen_2.bed | sed -e 's/id//' | sort -k4,4n | awk '{print $1, $2, $3, "id"$4, $5, $6}' | sed -e 's/ /\t/g' > pred_pita_results_c_gen.bed


	paste $FILE_PRED genomic.txt > pred_pita_results_d_per_R.txt

	bedtools intersect -a $AGO2 -b pred_pita_results_c_gen.bed -wo > pred_pita_results_AG02.bed.txt

	cat pred_pita_results_AG02.bed.txt | awk '$13>=5 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' | sed -e 's/ /\t/g' > pred_pita_results_AG02.bed


	#3# select only AGO2 overlapping sites from the original file

	cat pred_pita_results_AG02.bed | cut -f10 | sort | uniq > id_AGO2.txt

	cat pred_pita_results_c.txt | head -n1 > pred_pita_results_per_R.txt
	grep -f id_AGO2.txt pred_pita_results_d_per_R.txt | cut -f 1-13 >> pred_pita_results_per_R.txt

	echo "PITA miRNA binding sites overlapped by AGO2:" $(( $( cat pred_pita_results_per_R.txt | wc -l ) -1))


	rm genomic.txt
	rm genomic_1.txt
	rm genomic_2.txt
	rm circ_conversion.txt
	rm pred_pita_results_b_gen.bed
	rm pred_pita_results_b_gen_1.bed
	rm pred_pita_results_b_gen_2.bed
	rm pred_pita_results_c_gen.bed
	rm pred_pita_results_AG02.bed.txt
	rm pred_pita_results_AG02.bed
	rm id_AGO2.txt
	rm pred_pita_results_d_per_R.txt

	cd ..


	rm genomic_sequence_position_conversion.txt

else

	mv miRanda/output_miRanda_c_per_R.txt miRanda/output_miRanda_per_R.txt
	mv PITA/pred_pita_results_c.txt PITA/pred_pita_results_per_R.txt

fi

cd ..

