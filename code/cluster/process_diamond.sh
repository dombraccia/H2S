#!/bin/bash

## first part: taking the top hit for each read and saving in processed dir
## second part: distilling the read mappings into a count table for plotting

echo "- initializing variables"
diamond_out=$1/*.m8
diamond_processed=$2
mkdir -p results/diamond_processed

echo "- writing to processed directory"
for file in $diamond_out; do
	current_file=$(basename $file)
	cat $file | sort -u -k1,1 > $diamond_processed/$current_file
done

echo "- distilling data to a count table"
samples=$(ls results/diamond_processed | cut -d '_' -f1 | uniq)
diamond_processed=$2/*.m8
declare -i loop_count; loop_count=1
for sample in $samples; do
	## initializing read files 
	read1=results/diamond_processed/${sample}_1.m8
	read2=results/diamond_processed/${sample}_2.m8

	## counting occurences of each gene in each read file and adding them together
	declare -i mst_count; mst_count=$(cat $read1 | grep -o "mst" | wc -l)
	mst_count+=$(cat $read2 | grep -o "mst" | wc -l)
	declare -i cse_count; cse_count=$(cat $read1 | grep -o "cse" | wc -l)
	cse_count+=$(cat $read2 | grep -o "cse" | wc -l)
	declare -i cbs_count; cbs_count=$(cat $read1 | grep -o "cbs" | wc -l)
	cbs_count+=$(cat $read2 | grep -o "cbs" | wc -l)
	declare -i cyd_count; cyd_count=$(cat $read1 | grep -o "cyd" | wc -l)
	cyd_count+=$(cat $read2 | grep -o "cyd" | wc -l)
	declare -i mgl_count; mgl_count=$(cat $read1 | grep -o "mgl" | wc -l)
	mgl_count+=$(cat $read2 | grep -o "mgl" | wc -l)
	declare -i dsrA_count; dsrA_count=$(cat $read1 | grep -o "dsrA" | wc -l)
	dsrA_count+=$(cat $read2 | grep -o "dsrA" | wc -l)
	declare -i dsrB_count; dsrB_count=$(cat $read1 | grep -o "dsrB" | wc -l)
	dsrB_count+=$(cat $read2 | grep -o "dsrB" | wc -l)
	
	if [ $loop_count == 1 ]; then
		printf "%b" "Run\t" "MST\t" "CSE\t" "CBS\t" "CYD\t" "MGL\t" "dsrA\t" "dsrB\n" >> results/diamond_processed/counts.tsv
	fi

	## writing counts to tsv
	printf "%b" "$sample\t" "$mst_count\t" "$cse_count\t" "$cbs_count\t" \
		"$cyd_count\t" "$mgl_count\t" "$dsrA_count\t" "$dsrB_count\n" \
		>> results/diamond_processed/counts.tsv
	
	(( loop_count++ ))
done
