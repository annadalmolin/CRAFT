#!/bin/bash

## Main script for functional predictions of a set of circRNAs

## input:
# circRNA backsplice file: "list_backsplice.txt"
# file with the path of input files: "path_files.txt"
# file with CRAFT parameters: "params.txt"
#
#"path_files.txt":
#1 row: path to Ensembl annotation file
#2 row: path to Ensembl genome file
#3 row: path to miRNA sequences file
#
#"params.txt"
#1 row: PRED		  # PRED can be "M" for miRNA prediction, "R" for RBP prediction, "O" for ORF prediction; "MR", "MO", "RO" or "MRO" for a combination of the previous
#2 row: SPECIES	  # one of the species in miRBase: hsa, mmu, etc.
#3 row: PARAMS_miRanda	  # in order, separated by tab: miRanda_score miRanda_energy; default: 80 -15
#4 row: PARAMS_beRBP	  # in order, separated by tab: "PWM" "RBP" "prefix"; multiple PWMs (separated by ", ") and associated RBP (separated by ", ") are also allowed. Default: "all" "all"
#5 row: PREFIX_UCSC      # prefix of the UCSC genome used by beRBP; ex. "hg38"
#6 row: PARAMS_ORFfinder # in order, separated by tab: GEN_CODE START_CODON MIN_LENGTH NESTED_ORFS STRAND
#			  # GEN_CODE: genetic code to use (1-31), see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details; default: 1
#			  # START_CODON: start codon to use (0 = "ATG" only, 1 = "ATG" and alternative initiation codons, 2 = any sense codon); default: 0
#			  # MIN_LENGTH: minimal length of the ORF (nt); allowed values are: 30, 75, 150. Default: 30
#			  # NESTED_ORFS: ignore nested ORFs (ORFs completely placed within another); allowed values: "TRUE", "FALSE". Default: "FALSE"
#			  # STRAND: output ORFs on specified strand only; allowed values: "both", "plus", "minus". Default: "plus"
#7 row: PARAMS_CIRC	  # i.e. "l=50000, score_miRNA=125, energy_miRNA=-25, dGduplex_miRNA=-22, dGopen_miRNA=-10, voteFrac_RBP=0.3"
#8 row: PARAMS_GENERAL	  # i.e. "l=50000, score_miRNA=125, energy_miRNA=-25, dGduplex_miRNA=-22, dGopen_miRNA=-10, voteFrac_RBP=0.3"

## command: ./scripts/pipeline_predictions.sh /data/list_backsplice.txt /data/path_files.txt /data/params.txt



CIRC_LIST=$1
PATH_FILE=$2
PARAMS=$3



# 1) sequence extraction

ANNOTATION_ENS=$( cat $PATH_FILE | sed -n '1 p' )
GENOME_ENS=$( cat $PATH_FILE | sed -n '2 p' )

if [ ! -d "/data/sequence_extraction" ]
then

	mkdir /data/sequence_extraction
	cd /data/sequence_extraction

	/scripts/script_sequence_extraction.sh $CIRC_LIST $ANNOTATION_ENS $GENOME_ENS

	cd ..

else
	echo "Putative sequence/s already extracted or provided by the user."
fi



# 2) functional predictions

if [ ! -d "/data/functional_predictions" ]
then
	mkdir /data/functional_predictions
fi

cd /data/functional_predictions

if [ ! -e "backsplice_sequence_1.fa" ]
then
	mv /data/sequence_extraction/backsplice_sequence_1.fa .
fi

if [ ! -e "backsplice_sequence_1.txt" ]
then
	mv /data/sequence_extraction/backsplice_sequence_1.txt .
fi

if [ ! -e "backsplice_circRNA_length_1.txt" ]
then
	mv /data/sequence_extraction/backsplice_circRNA_length_1.txt .
	cat backsplice_circRNA_length_1.txt | grep -v "circ_id" > circRNA_length.txt
fi

if [ ! -e "backsplice_gene_name.txt" ]
then
	mv /data/input/backsplice_gene_name.txt .
fi


PRED=$( cat $PARAMS | sed -n '1 p' )

if [ $PRED == "M" -o $PRED == "MR" -o $PRED == "MO" -o $PRED == "MRO" ]
then

	SPECIES=$( cat $PARAMS | sed -n '2 p' )

	if [ ! -e /input/"mature_"$SPECIES".fa" ]
	then
		MIRNA_FILE="/input/mature_miRNA.txt"
		cat $MIRNA_FILE | grep $SPECIES | sed -e 's/\t/\n/' | sed -e 's/ MIMAT.*$//' > /data/input/.mature_$SPECIES.fa
		MIRNA_FILE="/data/input/.mature_"$SPECIES".fa"
	fi


	#PARAMS_miRanda=$( cat $PARAMS | sed -n '3 p' )
	SCORE=$( cat $PARAMS | sed -n '3 p' | cut -f1 )
	ENERGY=$( cat $PARAMS | sed -n '3 p' | cut -f2 )

	if [ -z $SCORE ]
	then
		SCORE=80
	fi

	if [ -z $ENERGY ]
	then
		ENERGY=-15
	fi

fi


if [ $PRED == "R" -o $PRED == "MR" -o $PRED == "RO" -o $PRED == "MRO" ]
then

	if [ ! -e "/tools/beRBP/lib/hg38.phyloP100way.bw" ]
	then
		wget -P /tools/beRBP/lib http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw
	fi


	#PARAMS_beRBP=$( cat $PARAMS | sed -n '4 p' )
	PWM=$( cat $PARAMS | sed -n '4 p' | cut -f1 )
	RBP=$( cat $PARAMS | sed -n '4 p' | cut -f2 )

	if [ -z $PWM ]
	then
		PWM="all"
	fi

	if [ -z $RBP ]
	then
		RBP="all"
	fi


	PREFIX_UCSC=$( cat $PARAMS | sed -n '5 p' )

fi


if [ $PRED == "O" -o $PRED == "MO" -o $PRED == "RO" -o $PRED == "MRO" ]
then

	#PARAMS_ORFfinder=$( cat $PARAMS | sed -n '6 p' )
	GEN_CODE=$( cat $PARAMS | sed -n '6 p' | cut -f1 )
	START_CODON=$( cat $PARAMS | sed -n '6 p' | cut -f2 )
	MIN_LENGTH=$( cat $PARAMS | sed -n '6 p' | cut -f3 )
	NESTED_ORFS=$( cat $PARAMS | sed -n '6 p' | cut -f4 )
	STRAND=$( cat $PARAMS | sed -n '6 p' | cut -f5 )

	if [ -z $GEN_CODE ]
	then
		GEN_CODE=1
	fi

	if [ -z $START_CODON ]
	then
		START_CODON=0
	fi

	if [ -z $MIN_LENGTH ]
	then
		MIN_LENGTH=30
	fi

	if [ -z $NESTED_ORFS ]
	then
		NESTED_ORFS="FALSE"
	fi

	if [ -z $STRAND ]
	then
		STRAND="plus"
	fi

fi


case $PRED in

    "M")
	if [ ! -d "/data/functional_predictions/miRNA_detection" ]
	then
		/scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY
		echo "MiRNA binding site prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/miRNA_detection/miRanda/output_miRanda_per_R.txt" ] | [ ! -e "/data/functional_predictions/miRNA_detection/PITA/pred_pita_results_per_R.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/miRNA_detection
			/scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY
			echo "MiRNA binding site prediction analysis completed."
		else
			echo "MiRNA binding site prediction analysis already performed."
		fi
	fi
	;;

    "R")
	if [ ! -d "/data/functional_predictions/RBP_detection" ]
	then
		if [ $SPECIES != "hsa" ]
		then
			echo "Human PWMs will be used instead for RBP binding."
		fi

		/scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC
		echo "RBP binding site prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/RBP_detection/beRBP/analysis_RBP/resultMatrix_b.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/RBP_detection

			if [ $SPECIES != "hsa" ]
			then
				echo "Human PWMs will be used instead for RBP binding."
			fi

			/scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC
			echo "RBP binding site prediction analysis completed."
		else
			echo "RBP binding site prediction analysis already performed."
		fi
	fi
	;;

    "O")
	if [ ! -d "/data/functional_predictions/ORF_detection" ]
	then
		/scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
		echo "ORF prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/ORF_detection/ORFfinder/ORF_backsplice.txt" ] | [ ! -e "/data/functional_predictions/ORF_detection/ORFfinder/ORF_backsplice_open.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/ORF_detection
			/scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
			echo "ORF prediction analysis completed."
		else
			echo "ORF prediction analysis already performed."
		fi
	fi
	;;

    "MR")
	if [ ! -d "/data/functional_predictions/miRNA_detection" ]
	then
		/scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY
		echo "MiRNA binding site prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/miRNA_detection/miRanda/output_miRanda_per_R.txt" ] | [ ! -e "/data/functional_predictions/miRNA_detection/PITA/pred_pita_results_per_R.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/miRNA_detection
			/scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY
			echo "MiRNA binding site prediction analysis completed."
		else
			echo "MiRNA binding site prediction analysis already performed."
		fi
	fi

	if [ ! -d "/data/functional_predictions/RBP_detection" ]
	then
		if [ $SPECIES != "hsa" ]
		then
			echo "Human PWMs will be used instead for RBP binding."
		fi

		/scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC
		echo "RBP binding site prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/RBP_detection/beRBP/analysis_RBP/resultMatrix_b.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/RBP_detection

			if [ $SPECIES != "hsa" ]
			then
				echo "Human PWMs will be used instead for RBP binding."
			fi

			/scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC
			echo "RBP binding site prediction analysis completed."
		else
			echo "RBP binding site prediction analysis already performed."
		fi
	fi
	;;

    "MO")
	if [ ! -d "/data/functional_predictions/miRNA_detection" ]
	then
		/scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY
		echo "MiRNA binding site prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/miRNA_detection/miRanda/output_miRanda_per_R.txt" ] | [ ! -e "/data/functional_predictions/miRNA_detection/PITA/pred_pita_results_per_R.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/miRNA_detection
			/scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY
			echo "MiRNA binding site prediction analysis completed."
		else
			echo "MiRNA binding site prediction analysis already performed."
		fi
	fi

	if [ ! -d "/data/functional_predictions/ORF_detection" ]
	then
		/scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
		echo "ORF prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/ORF_detection/ORFfinder/ORF_backsplice.txt" ] | [ ! -e "/data/functional_predictions/ORF_detection/ORFfinder/ORF_backsplice_open.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/ORF_detection
			/scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
			echo "ORF prediction analysis completed."
		else
			echo "ORF prediction analysis already performed."
		fi
	fi
	;;

    "RO")
	if [ ! -d "/data/functional_predictions/RBP_detection" ]
	then
		if [ $SPECIES != "hsa" ]
		then
			echo "Human PWMs will be used instead for RBP binding."
		fi

		/scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC
		echo "RBP binding site prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/RBP_detection/beRBP/analysis_RBP/resultMatrix_b.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/RBP_detection

			if [ $SPECIES != "hsa" ]
			then
				echo "Human PWMs will be used instead for RBP binding."
			fi

			/scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC
			echo "RBP binding site prediction analysis completed."
		else
			echo "RBP binding site prediction analysis already performed."
		fi
	fi

	if [ ! -d "/data/functional_predictions/ORF_detection" ]
	then
		/scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
		echo "ORF prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/ORF_detection/ORFfinder/ORF_backsplice.txt" ] | [ ! -e "/data/functional_predictions/ORF_detection/ORFfinder/ORF_backsplice_open.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/ORF_detection
			/scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
			echo "ORF prediction analysis completed."
		else
			echo "ORF prediction analysis already performed."
		fi
	fi
	;;

    "MRO")
	if [ ! -d "/data/functional_predictions/miRNA_detection" ]
	then
		/scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY
		echo "MiRNA binding site prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/miRNA_detection/miRanda/output_miRanda_per_R.txt" ] | [ ! -e "/data/functional_predictions/miRNA_detection/PITA/pred_pita_results_per_R.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/miRNA_detection
			/scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY
			echo "MiRNA binding site prediction analysis completed."
		else
			echo "MiRNA binding site prediction analysis already performed."
		fi
	fi

	if [ ! -d "/data/functional_predictions/RBP_detection" ]
	then
		if [ $SPECIES != "hsa" ]
		then
			echo "Human PWMs will be used instead for RBP binding."
		fi

		/scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC
		echo "RBP binding site prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/RBP_detection/beRBP/analysis_RBP/resultMatrix_b.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/RBP_detection

			if [ $SPECIES != "hsa" ]
			then
				echo "Human PWMs will be used instead for RBP binding."
			fi

			/scripts/script_RBP_detection.sh $PWM $RBP $PREFIX_UCSC
			echo "RBP binding site prediction analysis completed."
		else
			echo "RBP binding site prediction analysis already performed."
		fi
	fi

	if [ ! -d "/data/functional_predictions/ORF_detection" ]
	then
		/scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
		echo "ORF prediction analysis completed."
	else
		if [ ! -e "/data/functional_predictions/ORF_detection/ORFfinder/ORF_backsplice.txt" ] | [ ! -e "/data/functional_predictions/ORF_detection/ORFfinder/ORF_backsplice_open.txt" ]	# the analysis was not successful
		then
			rm -rf /data/functional_predictions/ORF_detection
			/scripts/script_ORF_prediction.sh $GEN_CODE $START_CODON $MIN_LENGTH $NESTED_ORFS $STRAND
			echo "ORF prediction analysis completed."
		else
			echo "ORF prediction analysis already performed."
		fi
	fi
	;;

    *)
	echo $PRED "is not a valid option. Try with one of the followings: M, R, O, MR, MO, RO, MRO."
	;;

esac


cd ..



# 3) graphical output

# prepare 1
if [ ! -d "/data/.R" ]
then
	mkdir /data/.R
fi

if [ ! -d "/data/.R/lib" ]
then
	mkdir /data/.R/lib
fi


# prepare 2
if [ ! -d "/data/.scripts" ]
then
	mkdir /data/.scripts
fi

cp /scripts/functional_predictions_* /data/.scripts/



if [ ! -d "/data/graphical_output" ]
then
	mkdir /data/graphical_output
fi

cd /data/graphical_output

cat $PARAMS | sed -n '1,2 p' > ../params_R.txt


# predictions for each circRNA singolarly

PARAMS_CIRC=$( cat $PARAMS | sed -n '7 p' )

if [ -z "$PARAMS_CIRC" ]
then
	PARAMS_CIRC="l=50000, score_miRNA=80, energy_miRNA=-15, dGduplex_miRNA=-15, dGopen_miRNA=-15, voteFrac_RBP=0.15"		# default
fi

for CIRC in $( cat $CIRC_LIST | cut -f1 )
do
	if [ ! -d $CIRC ]
	then
		mkdir $CIRC
	fi
	cd $CIRC

	R -e "rmarkdown::render('/data/.scripts/functional_predictions_single_circRNA.Rmd', output_file='functional_predictions_single_circRNA.html', output_dir='.', params = list(circ='$CIRC', $PARAMS_CIRC))"
	mv functional_predictions_single_circRNA.html functional_predictions_$CIRC.html

	cd ..
done



# prediction for all circRNAs together

PARAMS_GENERAL=$( cat $PARAMS | sed -n '8 p' )

if [ -z "$PARAMS_GENERAL" ]
then
	PARAMS_GENERAL="l=50000, score_miRNA=80, energy_miRNA=-15, dGduplex_miRNA=-15, dGopen_miRNA=-15, voteFrac_RBP=0.15"		# default
fi

if [ ! -d "general" ]
then
	mkdir general
fi
cd general

R -e "rmarkdown::render('/data/.scripts/functional_predictions_all_circRNAs.Rmd', output_file='functional_predictions_all_circRNAs.html', output_dir='.', params = list($PARAMS_GENERAL))"

cd ..


rm -rf /data/.scripts
rm ../params_R.txt


cd ..

