##### ======== WORKFLOW FOR TRANSCRIPTOMICS ANALYSIS FOR H2S PAPER ======== ###

rule download_txp_data:
	input: "data/david-2014/SraRunTable.txt"
	output: "data/david-2014/dummy/download.txt"
	shell: "code/cluster/download_SRA.sh {input} {output}"

rule makedb:
	input: "data/from-uniprot/H2S_producing_enzymes.fa"
	output: "data/from-uniprot/H2S_enzymes_db.dmnd"
	shell: "code/cluster/make_diamond_db.sh {input} {output}"

rule fastp:
	input: "data/david-2014"
	output: directory("results/fastp_out")
	shell: "code/cluster/run_fastp.sh {input} {output}"

rule align:
	input: "data/from-uniprot/H2S_enzymes_db.dmnd"
	output: directory("results/diamond_out") 
	shell: "code/cluster/align_david2014.sh {input} {output}"

rule process_diamond:
	input: directory("results/diamond_out")
	output: directory("results/diamond_processed")
	shell: "code/cluster/process_diamond.sh {input} {output}"
