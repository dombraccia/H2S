##### ======== WORKFLOW FOR TRANSCRIPTOMICS ANALYSIS FOR H2S PAPER ======== ###

rule download_txp_data:
	input: "data/SRR848994.txt"
	output: "data/from-SRA/dummy/download.txt"
	shell: "code/cluster/download_SRA.sh {input} {output}"
		