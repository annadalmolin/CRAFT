#!/bin/bash

## Script extracting the putative circRNA sequence from the backsplice coordinates

## input:
# circRNA backsplice file ("/data/list_backsplice.txt"): "$CIRC_LIST"
# annotation file (e.g. "/data/input/Homo_sapiens.GRCh38.93.gtf"): "$ANNOTATION_ENS"
# genome file (e.g. "/data/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"): "$GENOME_ENS"

## command: /scripts/script_sequence_extraction.sh $CIRC_LIST $ANNOTATION_ENS $GENOME_ENS



CIRC_LIST=$1
ANNOTATION_ENS=$2
GENOME_ENS=$3



## 1) backsplice file

# circIDs list

cat $CIRC_LIST | sed -e 's/:/\t/; s/-/\t/' | awk '{print $1,$2-1,$3,"backsplice",".",$4}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > circRNAs.bed


# start and stop position list

#circRNAs on strand "+" or "-" (strand 1):
cat circRNAs.bed | awk '(($6 == "+") || ($6 == "-")) {print}' > circRNAs_1.bed
cat circRNAs_1.bed | sed -r 's/^([0-9A-Z.]*)\t([0-9]*)\t([0-9]*)\t([a-z]*)\t(.)\t([-+|.]*)/echo "\1\t\2\t$((\2+1))\tstart\t\5\t\6\n\1\t$((\3-1))\t\3\tstop\t\5\t\6"/e' | sort -k1,1 -k2,2n -k3,3n | uniq > circRNAs_start_stop_1.bed

# circRNAs on strand "-|+" (intergenics or with uncertain annotation) (strand 2):
cat circRNAs.bed | awk '($6 == "-|+") {print}' | sed -r 's/^([0-9A-Z.]*)\t([0-9]*)\t([0-9]*)\t([a-z]*)\t(.)\t([-+|.]*)/echo "\1\t\2\t\3\t\4\t\5\t\"-\"\n\1\t\2\t\3\t\4\t\5\t\"+\""/e' > circRNAs_2.bed
cat circRNAs_2.bed | sed -r 's/^([0-9A-Z.]*)\t([0-9]*)\t([0-9]*)\t([a-z]*)\t(.)\t([-+|.]*)/echo "\1\t\2\t$((\2+1))\tstart\t\5\t\6\n\1\t$((\3-1))\t\3\tstop\t\5\t\6"/e' | sort -k1,1 -k2,2n -k3,3n | uniq > circRNAs_start_stop_2.bed





## 2) annotation files

# "merged_exons.bed"

if [ ! -e "merged_exons.bed" ]
then
	grep -w "exon" $ANNOTATION_ENS | cut -f 1,4,5,7 | awk '{print $1,$2-1,$3,"exon",".",$4}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > annotation_exon.bed
	bedtools merge -s -i annotation_exon.bed | awk '{print $1,$2,$3,"exon",".",$4}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > merged_exons.bed
fi


# "introns.bed"

if [ ! -e "introns.bed" ]
then
	grep -w "gene" $ANNOTATION_ENS | cut -f 1,4,5,7 | awk '{print $1,$2-1,$3,"gene",".",$4}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > annotation_gene.bed
	bedtools subtract -s -a annotation_gene.bed -b merged_exons.bed | sed -e 's/gene/intron/' | sort -k1,1 -k2,2n -k3,3n | uniq > introns.bed
fi


# merge of "merged_exons.bed" and "introns.bed"

if [ ! -e "exon_intron_sorted.bed" ]
then
	cat merged_exons.bed introns.bed | sort -k1,1 -k2,2n -k3,3n | uniq > exon_intron_sorted.bed
fi


# "intergenics.bed"

if [ ! -e "intergenics_for.bed" ] | [ ! -e "intergenics_rev.bed" ]
then
	cat $GENOME_ENS | grep -e ">" | sed -e 's/>//; s/ /\t/g; s/:/\t/g' | cut -f1,8 | sort -k1,1 > annotation_chr.genome
	cat annotation_gene.bed | awk '$6=="+" {print}' | sed 's/ /\t/' > annotation_gene_for.bed
	cat annotation_gene.bed | awk '$6=="-" {print}' | sed 's/ /\t/' > annotation_gene_rev.bed
	bedtools complement -i annotation_gene_for.bed -g annotation_chr.genome | awk '{print $1,$2,$3,"intergenic",".","+"}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > intergenics_for.bed	# strand "+" per l'intergenico
	bedtools complement -i annotation_gene_rev.bed -g annotation_chr.genome | awk '{print $1,$2,$3,"intergenic",".","-"}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > intergenics_rev.bed	# strand "-" per l'intergenico
fi


# merge of "exon_intron_sorted.bed" and "intergenics.bed"

if [ ! -e "exon_intron_intergenics_sorted.bed" ]
then
	cat exon_intron_sorted.bed intergenics_for.bed intergenics_rev.bed | sort -k1,1 -k2,2n -k3,3n | uniq > exon_intron_intergenics_sorted.bed
fi





## 3) retrieve backsplice regions

# strand 1:

bedtools intersect -s -loj -sorted -a circRNAs_start_stop_1.bed -b exon_intron_intergenics_sorted.bed -wa -wb | sort -k1,1 -k2,2n -k3,3n | uniq > circRNAs_start_stop_exon_intron_intergenic_1.bed
bedtools intersect -s -loj -a circRNAs_1.bed -b circRNAs_start_stop_exon_intron_intergenic_1.bed -wa -wb | sort -k1,1 -k2,2n -k3,3n | uniq > combined_circRNAs_1.bed

cat combined_circRNAs_1.bed | awk '$16 == "intron" {print $13,$14,$15,$16,$17,$18 }' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > backspliced_introns_1.bed
cat combined_circRNAs_1.bed | awk '$16 == "intergenic" {print $13,$14,$15,$16,$17,$18 }' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > backspliced_intergenics_1.bed
cat merged_exons.bed backspliced_introns_1.bed backspliced_intergenics_1.bed | sort -k1,1 -k2,2n -k3,3n | uniq > merged_exons_backspliced_introns_intergenics_1.bed


# strand 2:

bedtools intersect -s -loj -sorted -a circRNAs_start_stop_2.bed -b exon_intron_intergenics_sorted.bed -wa -wb | sort -k1,1 -k2,2n -k3,3n | uniq > circRNAs_start_stop_exon_intron_intergenic_2.bed
bedtools intersect -s -loj -a circRNAs_2.bed -b circRNAs_start_stop_exon_intron_intergenic_2.bed -wa -wb | sort -k1,1 -k2,2n -k3,3n | uniq > combined_circRNAs_2.bed

cat combined_circRNAs_2.bed | awk '$16 == "intron" {print $13,$14,$15,$16,$17,$18 }' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > backspliced_introns_2.bed
cat combined_circRNAs_2.bed | awk '$16 == "intergenic" {print $13,$14,$15,$16,$17,$18 }' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > backspliced_intergenics_2.bed
cat merged_exons.bed backspliced_introns_2.bed backspliced_intergenics_2.bed | sort -k1,1 -k2,2n -k3,3n | uniq > merged_exons_backspliced_introns_intergenics_2.bed





## 4) retrieve the regions intersected by the backsplice, for each backsplice

## each intersected fragment put on a different row, knowing if it is "exon", "intron" or "intergenic"

# strand 1:
bedtools intersect -s -a circRNAs_1.bed -b merged_exons_backspliced_introns_intergenics_1.bed -wb | cut -f1,2,3,10,11,12 | sed -e 's/exon/exon_intersected/; s/intron/intron_intersected/; s/intergenic/intergenic_intersected/' | sort -k1,1 -k2,2n -k3,3n | uniq > exon_intron_intergenic_intersected_1.bed

# strand 2:
bedtools intersect -s -a circRNAs_2.bed -b merged_exons_backspliced_introns_intergenics_2.bed -wb | cut -f1,2,3,10,11,12 | sed -e 's/exon/exon_intersected/; s/intron/intron_intersected/; s/intergenic/intergenic_intersected/' | sort -k1,1 -k2,2n -k3,3n | uniq > exon_intron_intergenic_intersected_2.bed



## each intersected region put on a different row, merging adjacent or overlapping fragments

# strand 1:
bedtools intersect -s -a circRNAs_1.bed -b exon_intron_intergenic_intersected_1.bed | sort -k1,1 -k2,2n -k3,3n | uniq | bedtools merge -s -i stdin | sort -k1,1 -k2,2n -k3,3n | uniq | awk '{print $1,$2,$3,"region_intersected",".",$4}' | sed -e 's/ /\t/g' > region_intersected_1.bed

# strand 2:
bedtools intersect -s -a circRNAs_2.bed -b exon_intron_intergenic_intersected_2.bed | sort -k1,1 -k2,2n -k3,3n | uniq | bedtools merge -s -i stdin | sort -k1,1 -k2,2n -k3,3n | uniq | awk '{print $1,$2,$3,"region_intersected",".",$4}' | sed -e 's/ /\t/g' > region_intersected_2.bed



## each backsplice with intersected regions put on a different row

# strand 1:
bedtools intersect -s -a region_intersected_1.bed -b circRNAs_1.bed -wb | awk '{print $7,$8,$9,$10,$11,$12, $1,$2,$3,$4,$5,$6}' | sed 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > backsplice_region_intersected_1.bed
cat backsplice_region_intersected_1.bed | sed -e 's/\t/:/ ; s/\t/-/' | cut -f1,2,5,6,7,8,10 > backsplice_region_intersected_1.txt

# strand 2:
bedtools intersect -s -a region_intersected_2.bed -b circRNAs_2.bed -wb | awk '{print $7,$8,$9,$10,$11,$12, $1,$2,$3,$4,$5,$6}' | sed 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > backsplice_region_intersected_2.bed
cat backsplice_region_intersected_2.bed | sed -e 's/\t/:/ ; s/\t/-/' | cut -f1,2,5,6,7,8,10 > backsplice_region_intersected_2.txt





## 5) extract the genomic circRNA sequences

## extract the sequence of each region

# strand 1:

cat backsplice_region_intersected_1.bed | sed -e 's/\t/:/ ; s/\t/-/' | cut -f1,5,6,7,10 | awk '{print $2,$3,$4,$1,".",$5}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > region_to_extract_1.bed

N1_plus=$( cat region_to_extract_1.bed | awk '$6 == "+" {print}' | wc -l )

if [ $N1_plus -gt 0 ]
then
	cat region_to_extract_1.bed | awk '$6 == "+" {print}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > region_to_extract_for_1.bed
	bedtools getfasta -fi $GENOME_ENS -bed region_to_extract_for_1.bed -s -name -tab -fo out_region_tmp_for_1.fa
	cat out_region_tmp_for_1.fa | nl | sort -k2,2 -k1,1n | awk '{print $2, $3}' | sed -e 's/ /\t/g' > out_region_for_1.fa
fi

N1_minus=$( cat region_to_extract_1.bed | awk '$6 == "-" {print}' | wc -l )

if [ $N1_minus -gt 0 ]
then
	cat region_to_extract_1.bed | awk '$6 == "-" {print}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > region_to_extract_rev_1.bed
	bedtools getfasta -fi $GENOME_ENS -bed region_to_extract_rev_1.bed -s -name -tab -fo out_region_tmp_rev_1.fa
	cat out_region_tmp_rev_1.fa | nl | sort -k2,2 -k1,1n | awk '{print $2, $3}' | sed -e 's/ /\t/g' | tac > out_region_rev_1.fa
fi


# strand 2:

cat backsplice_region_intersected_2.bed | sed -e 's/\t/:/ ; s/\t/-/' | cut -f1,5,6,7,10 | awk '{print $2,$3,$4,$1,".",$5}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > region_to_extract_2.bed

N2_plus=$( cat region_to_extract_2.bed | awk '$6 == "+" {print}' | wc -l )

if [ $N2_plus -gt 0 ]
then
	cat region_to_extract_2.bed | awk '$6 == "+" {print}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > region_to_extract_for_2.bed
	bedtools getfasta -fi $GENOME_ENS -bed region_to_extract_for_2.bed -s -name -tab -fo out_region_tmp_for_2.fa
	cat out_region_tmp_for_2.fa | nl | sort -k2,2 -k1,1n | awk '{print $2, $3}' | sed -e 's/ /\t/g' > out_region_for_2.fa
fi

N2_minus=$( cat region_to_extract_2.bed | awk '$6 == "-" {print}' | wc -l )

if [ $N2_minus -gt 0 ]
then
	cat region_to_extract_2.bed | awk '$6 == "-" {print}' | sed -e 's/ /\t/g' | sort -k1,1 -k2,2n -k3,3n | uniq > region_to_extract_rev_2.bed
	bedtools getfasta -fi $GENOME_ENS -bed region_to_extract_rev_2.bed -s -name -tab -fo out_region_tmp_rev_2.fa
	cat out_region_tmp_rev_2.fa | nl | sort -k2,2 -k1,1n | awk '{print $2, $3}' | sed -e 's/ /\t/g' | tac > out_region_rev_2.fa
fi



## retrieve the sequence of each backsplice

# strand 1:

if [ $N1_plus -gt 0 ]
then
	bedtools groupby -i out_region_for_1.fa -g 1 -c 2 -o concat > backsplice_sequence_bed_for_1.txt
fi

if [ $N1_minus -gt 0 ]
then
	bedtools groupby -i out_region_rev_1.fa -g 1 -c 2 -o concat > backsplice_sequence_bed_rev_1.txt
fi

if [ $N1_plus -gt 0 ] && [ $N1_minus -gt 0 ]
then
	cat backsplice_sequence_bed_for_1.txt backsplice_sequence_bed_rev_1.txt | sort -k1,1n > backsplice_sequence_bed_1.txt
else
	if [ $N1_plus -gt 0 ]
	then
		cat backsplice_sequence_bed_for_1.txt | sort -k1,1n > backsplice_sequence_bed_1.txt
	else
		if [ $N1_minus -gt 0 ]
		then
			cat backsplice_sequence_bed_rev_1.txt | sort -k1,1n > backsplice_sequence_bed_1.txt
		fi
	fi
fi

if [ $N1_plus -gt 0 ] || [ $N1_minus -gt 0 ]
then
	cat backsplice_sequence_bed_1.txt | sed -r 's/^([0-9A-Z]*):([0-9]*)-([0-9]*)\t([A-Z]*)/echo "\1:$((\2+1))-\3\t\4"/e' > backsplice_sequence_1.txt	#doesn't work for sequences > 150000 nt
	cat backsplice_sequence_1.txt | sed -e 's/^/>/; s/\t/\n/' > backsplice_sequence_1.fa
fi


# strand 2:

if [ $N2_plus -gt 0 ]
then
	bedtools groupby -i out_region_for_2.fa -g 1 -c 2 -o concat > backsplice_sequence_bed_for_2.txt
fi

if [ $N2_minus -gt 0 ]
then
	bedtools groupby -i out_region_rev_2.fa -g 1 -c 2 -o concat > backsplice_sequence_bed_rev_2.txt
fi

if [ $N2_plus -gt 0 ] && [ $N2_minus -gt 0 ]
then
	cat backsplice_sequence_bed_for_2.txt backsplice_sequence_bed_rev_2.txt | sort -k1,1n > backsplice_sequence_bed_2.txt
else
	if [ $N2_plus -gt 0 ]
	then
		cat backsplice_sequence_bed_for_2.txt | sort -k1,1n > backsplice_sequence_bed_2.txt
	else
		if [ $N2_minus -gt 0 ]
		then
			cat backsplice_sequence_bed_rev_2.txt | sort -k1,1n > backsplice_sequence_bed_2.txt
		fi
	fi
fi

if [ $N2_plus -gt 0 ] || [ $N2_minus -gt 0 ]
then
	cat backsplice_sequence_bed_2.txt | sed -r 's/^([0-9A-Z]*):([0-9]*)-([0-9]*)\t([A-Z]*)/echo "\1:$((\2+1))-\3\t\4"/e' > backsplice_sequence_2.txt	#doesn't work for sequences > 150000 nt
	cat backsplice_sequence_2.txt | sed -e 's/^/>/; s/\t/\n/' > backsplice_sequence_2.fa
fi





## 6) retrieve circRNA length

# strand 1:

if [ $N1_plus -gt 0 ] || [ $N1_minus -gt 0 ]
then
	cat backsplice_sequence_1.fa | grep ">" | sed -e 's/>//' > circ_id.txt
	cat backsplice_sequence_1.fa | grep -v ">" | awk '{print length}' > circ_length.txt

	echo "circ_id length" | sed -e 's/ /\t/g' > backsplice_circRNA_length_1.txt
	paste circ_id.txt circ_length.txt >> backsplice_circRNA_length_1.txt
fi


# strand 2:

if [ $N2_plus -gt 0 ] || [ $N2_minus -gt 0 ]						
then
	cat backsplice_sequence_2.fa | grep ">" | sed -e 's/>//' > circ_id.txt
	cat backsplice_sequence_2.fa | grep -v ">" | awk '{print length}' > circ_length.txt

	echo "circ_id length" | sed -e 's/ /\t/g' > backsplice_circRNA_length_2.txt
	paste circ_id.txt circ_length.txt >> backsplice_circRNA_length_2.txt
fi





## 7) remove intermediary and empty files

rm annotation*
rm backsplice_int*
rm backsplice_region_*
rm backsplice_sequence_bed*
rm circRNAs.*
rm circ_id.txt
rm circ_length.txt
rm combined*
rm exon_intron_intergenic_intersected*
#rm exon_intron_intergenic_sorted.bed
rm exon_intron_sorted.bed
rm int*
rm merged*
rm out_region*
rm region_intersected*
rm region_to_extract_for_?.bed
rm region_to_extract_rev_?.bed

find . -type f -empty -delete

