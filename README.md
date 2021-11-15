# CRAFT

CRAFT is a computational pipeline that predicts circRNA sequence and molecular interactions with miRNAs and RBPs, along with their coding potential. CRAFT provides a comprehensive graphical visualization of the results, links to several knowledge databases, extensive functional enrichment analysis and combination of predictions for different circRNAs. CRAFT is a useful tool to help the user explore the potential regulatory networks involving the circRNAs of interest and generate hypotheses about the cooperation of circRNAs into the regulation of biological processes.

## Installation

### Installation from the Docker image

The Docker image saves you from the installation burden. A Docker image of CRAFT is available from DockerHub at https://hub.docker.com/r/annadalmolin/craft; just pull it with the command:

	docker pull annadalmolin/craft:v1.0

## Usage

### Input data

Prepare your project directory with the following files:

- _list_backsplice.txt_: file with circRNA coordinates. The file format is a tab-separated text file, with circRNA backsplice coordinates in the first column and circRNA strand in the second. An example of _list_backsplice.txt_ is:

		4:143543509-143543972	+
		11:33286413-33287511	+
		15:64499292-64500166	+

- _path_files.txt_: file with the relative paths for Ensembl annotation and genome files. The file format is a text file with a path written in each row, __in the following order__:

	1. path to annotation file
	2. path to genome file

    An example of _path_files.txt_ is:
	
		/data/input/Homo_sapiens.GRCh38.104.gtf
		/data/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa

	The gene annotation (in GTF format) and the genome sequence (in FASTA format) files must be downloaded by the user from Ensembl database and placed into the _input/_ directory contained in the project directory. Annotation and genome files for _Homo sapiens_ (GRCh38) can be downloaded from http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/ and http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/, respectively.

- _params.txt_: file with the parameters to be setted in CRAFT. The file format is a text file with a/more parameter/s written in each row, __in the following order__:

	1. kind of prediction; it can be "M" for miRNA prediction, "R" for RBP prediction, "O" for ORF prediction, "MR", "MO", "RO" or "MRO" for a combination of the previous.
	2. investigated species; it can be one of the species in miRBase database: _hsa_ for _Homo sapiens_, _mmu_ for _Mus musculus_, etc.
	3. parameters for miRanda tool (optional); in a single row, they must be the _miRanda_score_ and the _miRanda_energy_, __in order__, separated by tab. The user must set or both parameters or neither of the two; default values are 80 (score) and -15 (energy).
	4. parameters for beRBP tool (optional); in a single row, in order and separated by a tab, they must be the _PWM/s_ and the _RBP/s_ investigated. The syntax is: "PWM" "RBP"; multiple PWMs (separated by ", ") and associated RBP (separated by ", ") are also allowed. The default is _"all" "all"_, searching for all PWMs and RBPs included in beRBP database. The user must set both parameters or none of the two.
	5. prefix of the genome and indexes downloaded from UCSC website; f.i. _hg38_ for _Homo sapiens_. The human genome file can be downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/ . Genome and indexes must be included in the _input/_ directory.
	6. parameters for ORFfinder tool (optional); __in order__, separated by tab, the user must specify: the genetic code to use, the start codon to use, the minimal ORF length, whether to ignore nested ORFs and the strand in which putative ORFs are searched. The user must set all parameters or none of them. The allowed options for each parameter are:
		1. genetic code: 1-31, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details; default: 1
		2. start codon: 0 = "ATG" only, 1 = "ATG" and alternative initiation codons, 2 = any sense codon; default: 0
		3. minimal ORF length (nt): allowed values are 30, 75, or 150; default: 30
		4. ignore nested ORFs (ORF completely placed within another). allowed values are "TRUE" or "FALSE"; default: "FALSE"
		5. strand (output ORFs on specified strand only): allowed values are "both", "plus" or "minus"; default: "plus"
	
	7. parameters for the graphical output for a single circRNA investigated (optional, but advised); the default parameters are: _l=50000, QUANTILE1=”FALSE”, thr1=0.95, score_miRNA=120, energy_miRNA=-22, QUANTILE2=”FALSE”, thr2=0.95, dGduplex_miRNA=-20, dGopen_miRNA=-11, QUANTILE3=”FALSE”, thr3=0.9, voteFrac_RBP=0.15, orgdb="org.Hs.eg.db", meshdb="MeSH.Hsa.eg.db", symbol2eg="org.Hs.egSYMBOL2EG", eg2uniprot="org.Hs.egUNIPROT", org="hsapiens"_. The user must specify __only__ the parameters to be changed with respect to the default, in a comma-separated list format; the parameter order does not matter.
	Available parameters:
		1. _l_: maximum length of circRNAs analyzed
		2. QUANTILE: whether to filter predictions based on a quantile threshold (thr); _QUANTILE1_ and _thr1_ are set for miRanda predictions, _QUANTILE2_ and _thr2_ for PITA predictions, _QUANTILE3_ and _thr3_ for beRBP predictions
		3. _score_miRNA_ and _energy_miRNA_: respectively, _score_ and _energy_ values of miRanda tool. Best predictions are obtained with higher _score_ and lower _energy_
		4. _dGduplex_miRNA_ and _dGopen_miRNA_: respectively, _dGduplex_ and _dGopen_ values of PITA tool. Best predictions are obtained with lower _dGduplex_ and higher _dGopen_
		5. _voteFrac_RBP_: _voteFrac_ value of beRBP tool. Best predictions are obtained with higher _voteFrac_
		6. _orgdb_ and _meshdb_: databases for miRNA enrichment analysis; the default values are “org.Hs.eg.db” and “MeSH.Hsa.eg.db”, respectively (_Homo sapiens_)
		7. _symbol2eg_ and _eg2uniprot_: databases for RBP enrichment analysis; the default values are “org.Hs.egSYMBOL2EG” and “org.Hs.egUNIPROT”, respectively (_Homo sapiens_)
		8. _org_: organism, in the form: human - ’hsapiens’, mouse - ’mmusculus’; the default value is for _Homo sapiens_

	8. parameters for the summary graphical output for all circRNAs investigated (optional, but advised); the default parameters are the same as the previous point. The user must specify __only__ the parameters to be changed with respect to the default, in a comma-separated list format; the parameter order does not matter. Available parameters: the same as before, except for _meshdb_ and _org_. It is advised to set point 7 and point 8 parameters with the same values.

	An example of _params.txt_ file is:

		M
		hsa


		hg38
		
		score_miRNA=125, energy_miRNA=-25, dGduplex_miRNA=-22, dGopen_miRNA=-10
		score_miRNA=125, energy_miRNA=-25, dGduplex_miRNA=-22, dGopen_miRNA=-10, voteFrac_RBP=0.3

and directory:

- _input/_: directory containing the following files:

	- genome and annotation files from Ensembl database, and genome and indexes files from UCSC databases (see before)
	- _backsplice_gene_name.txt_: file with circRNA gene names. __It must be created by the user__. The file format is a tab-separated text file, with circRNA backsplice in the first column and circRNA host gene name in the second; the official gene name has to be used. The header line is needed.  An example of _backsplice_gene_name.txt_ is:

			circ_id	gene_names
			4:143543509-143543972	SMARCA5
			11:33286413-33287511	HIPK3
			15:64499292-64500166	ZNF609

	- _AGO2_binding_sites.bed_ (optional): file with validated AGO2 binding sites. The file, in BED6 format, must have the following fields: chromosome, start genomic position (0-based), end genomic position, the string “AGO2_binding_site”, a dot, the strand. Keep attention to use the same genome reference version as that included in the _input/_ directory. An example of _AGO2_binding_sites.bed_ is:

			4    143543521    143543542    AGO2_binding_site    .    +
			4    143543530    143543559    AGO2_binding_site    .    +
			4    143543562    143543607    AGO2_binding_site    .    +
			
   		The number of miRNA binding sites overlapped with AGO2 binding sites is written in the standard output. Check it in order to decide to keep AGO2 overlapping or re-running the analysis without this information (f.i. when very few sites are overlapping).
			

### Running the analysis

To run CRAFT from the Docker container use:

	sudo docker run -it -v $(pwd):/data annadalmolin/craft:v1.0
	
All paths in _path_files.txt_ must be relative to the directory in the container where the volumes were mounted (f.i. _/data/input/file_name_, as detailed above). 
If you want the container to give your user permissions, you need to set the owner id with "-u `id -u`":

	sudo docker run -u `id -u` -it -v $(pwd):/data annadalmolin/craft:v1.0


### Output data

After CRAFT successful run end, you will find the following new directories in your project directory:

1. _sequence_extraction/_: contains intermediary files for the sequence reconstruction step
2. _functional_predictions/_: contains final files of sequence reconstruction step and the three directories for miRNA, RBP and ORF predictions, respectively
3. _graphical_output/_: contains the directory _general_ with the summary predictions of all circRNA analyzed, and a directory for each single circRNA with the specific investigation

- __sequence_extraction/__

	The output files for the sequence reconstruction step are:
	
	- _backsplice_sequence_1.fa_: file with the retrieved genomic sequence for each circRNA in FASTA format
	- _backsplice_sequence_1.txt_: tab-separated file with the retrieved genomic sequence for each circRNA in TXT format; the file appear with the circRNA backsplice coordinates in the first column and the sequence in the second
	- _backsplice_circRNA_length_1.txt_: tab-separated file with circRNA sequence length, with circRNA backsplice in the first column and circRNA length in the second

	All these files are found in the _functional_predictions/_ directory.
	
- __functional_predictions/__

	The output files of functional prediction step are (the final output of each tool is highlighted in bold):
	
	- _miRNA_detection/_:
		- _backsplice_sequence_per_miRNA.fa_: the sequence used for miRNA prediction, obtained repeating the first 20 nt of the sequence at the end of each circRNA
		- _miRanda/_:
			- _output_miRanda.txt_: original output of miRanda
			- _output_miRanda_c_per_R.txt_: output of miRanda (list of miRNA binding sites), not overlapping with AGO2 binding sites, if _AGO2_binding_sites.bed_ is provided, otherwise this file is missing
			- __*output_miRanda_per_R.txt*__: final output of miRanda (list of miRNA binding sites), overlapping with AGO2 binding sites if _AGO2_binding_sites.bed_ is provided, otherwise it contains the list of miRNA binding sites not overlapping with AGO2 binding sites
		- _PITA/_:
			- _pred_pita_results.tab, pred_pita_results_targets.tab, pita.err, pita.log, pred_pita_results.gxp_: original output of PITA
			- _pred_pita_results_targets_b.txt_: output for multiple sites
			- _pred_pita_results_c.txt_: output of PITA (list of miRNA binding sites), not overlapping with AGO2 binding sites, if _AGO2_binding_sites.bed_ is provided, otherwise this file is missing
			- __*pred_pita_results_per_R.txt*__: final output of PITA (list of miRNA binding sites), overlapping with AGO2 binding sites if _AGO2_binding_sites.bed_ is provided, otherwise it contains the list of miRNA binding sites not overlapping with AGO2 binding sites

	- _RBP_detection/_:
		- _backsplice_sequence_per_RBP.fa_: the sequence used for RBP prediction, obtained repeating the first 20 nt of the sequence at the end of each circRNA
		- _beRBP/_:
			- _analysis_RBP/_:
				- _resultMatrix.tsv_: original output of beRBP
				- _resultMatrix_b.tsv_: final output of beRBP in TSV (list of RBP binding sites)
				- _resultMatrix_b.txt_: final output of beRBP in TXT (list of RBP binding sites)

	- _ORF_detection/_:
		- _backsplice_sequence_per_ORF_MIN_LENGTH.fa_: the sequence used for ORF prediction (with minimal length of the ORF = _MIN_LENGTH_), obtained doubling circRNA sequence twice
		- _ORFfinder/_:
			- _result_list_ORF_MIN_LENGTH.txt, result_list_CDS_MIN_LENGTH.txt, result_text_ORF_MIN_LENGTH.txt, result_table_ORF_MIN_LENGTH.txt, ORF0_MIN_LENGTH.log, ORF1_MIN_LENGTH.log, ORF2_MIN_LENGTH.log, ORF3_MIN_LENGTH.log, ORF0_MIN_LENGTH.perf, ORF1_MIN_LENGTH.perf, ORF2_MIN_LENGTH.perf, ORF3_MIN_LENGTH.perf_: original output of ORFfinder (with minimal length of the ORF = _MIN_LENGTH_)
			- __*ORF_backsplice.txt*__ and __*ORF_backsplice0.txt*__: final output of ORFfinder (list of ORF detected), respectively with ORF start position in 1-based and in 0-based format
			- __*ORF_backsplice_open.txt*__ and __*ORF_backsplice_open0.txt*__: final output of ORFfinder (list of rolling ORF detected), respectively with ORF start position in 1-based and in 0-based format
			- __*result_list_CDS.fa*__ and _result_list_CDS.txt_: nucleotidic ORF sequence, respectively in FASTA and TXT format
			- __*result_list_ORF.fa*__ and _result_list_ORF.txt_: amino acid ORF sequence, respectively in FASTA and TXT format

- __graphical_output/__

	The output files for the graphical output step are:

	- _general/_: directory with the summary predictions of all circRNA analyzed:
		- _functional_predictions_all_circRNAs.html_: output HTML file summarizing all predictions of all circRNA tested (see CRAFT paper for more details)
		- single figures pulled out from the HTML file
		- _All_validated_TGs.csv_: table pulled out from the HTML file; it can be loaded into Cytoscape for network analysis
	- a directory for each single circRNA with it own predictions:
		- _functional_predictions_CIRC_ID.html_: output HTML file with the predictions related to _CIRC_ID_ (see CRAFT paper for more details)
		- single figures and tables pulled out from the HTML file


### CircRNA sequence provided by the user

If circRNA sequences are available to the user, CRAFT doesn’t perform the sequence reconstruction step. So, to let CRAFT use the provided circRNA sequences, the user must follow these steps:

1. create the _sequence_extraction/_ directory into the project directory
2. add the _backsplice_sequence_1.fa_, _backsplice_sequence_1.txt_ and _backsplice_circRNA_length_1.txt_ files, in the format described above, to _sequence_extraction/_
3. add the _backsplice_gene_name.txt_ file, in the format described above, to _sequence_extraction/_
4. if the user wants to filter for miRNA binding sites overlapped with AGO2 binding sites, he/she must also add the file _region_to_extract_1.bed_ to _sequence_extraction/_. The file in BED6 format must have six columns tab-separated: circRNA chromosome, 0-based start position, 1-based end position, backsplice coordinates, score, strand. Each row represents a single separated region from which the circRNA is arranged (exon, intron, part of exon/intron or intergenic region). An example of _region_to_extract_1.bed_ is:

		11	33286412	33287511	11:33286412-33287511	.	+
		15	64499291	64500166	15:64499291-64500166	.	+
		4	143543508	143543657	4:143543508-143543972	.	+
		4	143543852	143543972	4:143543508-143543972	.	+
		

## Additional notes

- Functional enrichments on validated target genes of miRNAs with predicted binding sites in circRNA sequences can be performed only for _Homo sapiens_ (hsa), _Mus musculus_ (mmu) and _Rattus norvegicus_ (rno) species.
- The output clearness and intelligibility improve at the growing of filtering stringency; f.i., if a figure is not understandable or CRAFT crashes due to too many predictions, simply re-run the graphical part of the analysis increasing CRAFT stringency.
