##### ======== WORKFLOW FOR TRANSCRIPTOMICS ANALYSIS FOR H2S PAPER ======== ###

## checking for expression of methanogens genes in the david 2014 paper
function = ["H2S_producing",
			"CH4_producing"]

rule download_txp_data:
	input: "data/david-2014/SraRunTable.txt"
	output: "data/david-2014/dummy/download.txt"
	shell: "code/cluster/download_SRA.sh {input} {output}"

rule makedb:
	input: "data/from-uniprot/{function}_enzymes.fa"
	output: "data/from-uniprot/{function}_enzymes_db.dmnd"
	shell: "code/cluster/make_diamond_db.sh {input} {output}"

rule align:
	input: "data/from-uniprot/{function}_enzymes_db.dmnd"
	output: directory("results/diamond_out/{function}") 
	shell: "code/cluster/align_david2014.sh {input} {output}"

rule process_diamond:
	input: directory("results/diamond_out/{function}")
	output: directory("results/diamond_processed/{function}")
	shell: "code/cluster/process_diamond_CH4_producing.sh {input} {output}"

## BASEMENT: where I keep potentially unuseful rules
rule fastp:
	input: "data/david-2014"
	output: directory("results/fastp_out")
	shell: "code/cluster/run_fastp.sh {input} {output}"
