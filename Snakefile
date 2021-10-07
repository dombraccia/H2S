##### ======== WORKFLOW FOR TRANSCRIPTOMICS ANALYSIS FOR H2S PAPER ======== ###

## checking for expression of methanogens genes in the david 2014 paper
function = ["CysDeg",
			"DSR",
			"ch4"]

dataset = ["david2014",
		   "hpfs"]

# for downloading HPFS metatx data: 
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

# input genes from homologue search
rule index_cys_deg_genes: 
	input: "results/from-GutFunFind/Cysteine_Degradation_hmm.gene_seqs.fasta"
	output: directory("results/salmon_out/Cysteine_Degradation_hmm.gene_seqs_index")
	shell: "code/cluster/index_ntseqs.sh {input} {output}"

# mined from Methanobrevibacter smithii genomes in UHGG
rule mine_ch4_prod_genes:
	input: 
		m_smithii_1 = "/fs/cbcb-lab/hall/data/UHGG/uhgg_catalogue/MGYG-HGUT-005/MGYG-HGUT-00522/pan-genome/pan-genome.fna",
		m_smithii_2 = "/fs/cbcb-lab/hall/data/UHGG/uhgg_catalogue/MGYG-HGUT-021/MGYG-HGUT-02163/pan-genome/pan-genome.fna",
		m_smithii_3 = "/fs/cbcb-lab/hall/data/UHGG/uhgg_catalogue/MGYG-HGUT-024/MGYG-HGUT-02446/pan-genome/pan-genome.fna",
		m_smithii_4 = "/fs/cbcb-lab/hall/data/UHGG/uhgg_catalogue/MGYG-HGUT-035/MGYG-HGUT-03516/pan-genome/pan-genome.fna"
	output: "data/ch4_production_gene_seqs.fasta"
	shell: 
		'''
		rm -f data/ch4_producing_gene_seqs.fasta
		python code/cluster/mine_ch4_prod_genes.py \
			{input.m_smithii_1} {input.m_smithii_2} {input.m_smithii_3} {input.m_smithii_4} {output}
		'''

# rule get_ntseqs:
# 	input: 
# 		bedfile = "data/dsrAB_gene_positions_uhgg.bed",
# 		taxa_hits = "data/dsrAB_taxa_hits.txt",
# 		uhgg_seqs = "/fs/cbcb-lab/hall/data/UHGG/uhgg_catalogue/",
# 		uhgg_accession_lineage = "/fs/cbcb-lab/hall/data/UHGG/uhgg_catalogue/MGnify_accession_Lineage.txt"
# 	output: "data/dsrAB_gene_seqs.fasta"
# 	shell:
# 		"""
# 		rm -f {output}
# 		python code/cluster/get_ntseqs_dsrAB.py {input.bedfile} {input.taxa_hits} {input.uhgg_seqs} {input.uhgg_acc_lin}
# 		"""

# for indexing dsrAB and genes involved in ch4 production 
rule index_curated_gene_catalogue: 
	input: "data/{function}.gene_seqs.fasta"
	output: directory("results/salmon_out/{function}_gene_seqs_index")
	shell: "sbatch code/cluster/index_ntseqs.sh {input} {output}"

# TODO: rule does not work, need to run the line of code next to "shell" in a terminal window
rule quantify:
	input: 
		index = "results/salmon_out/{function}_gene_seqs_index",
		dataset = "data/{dataset}/",
		script = "code/cluster/quantify_{dataset}.sh"
	output: "results/salmon_out/{dataset}/dummy/{function}.quantify.txt"
	shell: "sbatch {input.script} {input.index} {input.dataset} {output}"

rule process_salmon_out:
	input: 
		dummy = "results/salmon_out/{dataset}/dummy/{function}.quantify.txt",
		dataset = "data/{dataset}",
		function = directory("results/salmon_out/{function}_gene_seqs_index")
	output: 
		tpm = "results/salmon_processed/{dataset}/{function}_tpm.tsv",
		counts = "results/salmon_processed/{dataset}/{function}_counts.tsv"
	shell: "python code/cluster/process_salmon_out.py {input.dataset} {input.function}"

rule get_nreads:
	input: "data/{dataset}/rnaseq"
	output: "data/{dataset}/nreads.txt"
	shell: "sbatch code/cluster/get_nreads.sh {input}"

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
