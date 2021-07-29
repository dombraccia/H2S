##### ======== WORKFLOW FOR TRANSCRIPTOMICS ANALYSIS FOR H2S PAPER ======== ###

## checking for expression of methanogens genes in the david 2014 paper
function = ["H2S_producing",
			"CH4_producing"]

dataset = ['david-2014',
		   'hpfs']

# for downloading metatx data: 
#	1. navigate to: https://www.ncbi.nlm.nih.gov/Traces/study/?
#	2. search for accession number PRJNA354235
#	3. in the `Select: Total` section, click the link to download the Metadata file
#	   named SraRunTable.txt
#	now the `download_txp_data` rule can be run.

rule download_txp_data:
	input: "data/{dataset}/SraRunTable.txt"
	output: "data/{dataset}/dummy/download.txt"
	shell: "code/cluster/download_SRA.sh {input} {output}"

rule quality_control:
	input: 
		datafile = 'data/{dataset}/rnaseq',
		dummy = 'data/{dataset}/dummy/download.txt'
	output: 'data/{dataset}/dummy/qc.txt'
	shell: 'sbatch code/cluster/quality_control.sh {input.datafile}'

rule index:
	input: "results/from-GutFunFind/Cysteine_Degradation_hmm.gene_seqs.fasta"
	output: directory("results/salmon_out/Cysteine_Degradation_hmm.gene_seqs_index")
	shell: "code/cluster/index_ntseqs.sh {input} {output}"

rule quantify:
	input: "results/salmon_out/Cysteine_Degradation_hmm.gene_seqs_index"
	output: "results/from-GutFunFind/dummy/quantify.txt"
	shell: "sbatch code/cluster/quantify_hpfs.sh {input} {output}"

rule process_salmon_out:
	input: "results/from-GutFunFind/dummy/quantify.txt"
	output: 
		tpm = "results/salmon_processed/tpm.tsv",
		counts = "results/salmon_processed/counts.tsv"
	shell: "code/cluster/process_salmon_out.py"

##### =========================== BASEMENT ============================== #####
# old workflow when only considering david 2014 data for metatranscriptomic 
# analysis 
##### =================================================================== #####

# for downloading david et al. 2014 data: 
#	1. navigate to: https://www.ncbi.nlm.nih.gov/Traces/study/?
#	2. search for accession number PRJNA202303
#	3. in the `Select: Total` section, click the link to download the Metadata file
#	   named SraRunTable.txt
#	now the `download_txp_data` rule can be run.

# rule download_txp_data:
# 	input: "data/david-2014/SraRunTable.txt"
# 	output: "data/david-2014/dummy/download.txt"
# 	shell: "code/cluster/download_SRA.sh {input} {output}"

# rule makedb:
# 	input: "data/from-uniprot/{function}_enzymes.fa"
# 	output: "data/from-uniprot/{function}_enzymes_db.dmnd"
# 	shell: "code/cluster/make_diamond_db.sh {input} {output}"

# rule align:
# 	input: "data/from-uniprot/{function}_enzymes_db.dmnd"
# 	output: directory("results/diamond_out/{function}") 
# 	shell: "code/cluster/align_david2014.sh {input} {output}"

# rule process_diamond:
# 	input: directory("results/diamond_out/{function}")
# 	output: directory("results/diamond_processed/{function}")
# 	shell: "code/cluster/process_diamond_H2S_producing.sh {input} {output}"

# rule fastp:
# 	input: "data/david-2014"
# 	output: directory("results/fastp_out")
# 	shell: "code/cluster/run_fastp.sh {input} {output}"
