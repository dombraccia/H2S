#!/bin/bash

## first part: taking the top hit for each read and saving in processed dir
## second part: distilling the read mappings into a count table for plotting

echo "- initializing variables"
diamond_out=$1/*.m8
diamond_processed=$2
mkdir -p $2

echo "- writing to processed directory"
for file in $diamond_out; do
	current_file=$(basename $file)
	cat $file | sort -u -k1,1 > $diamond_processed/$current_file
done

echo "- distilling data to a count table"
samples=$(ls $diamond_processed | cut -d '_' -f1 | uniq)
diamond_processed_m8=$2/*.m8
declare -i loop_count; loop_count=1
for sample in $samples; do
	## initializing read files 
	echo "-- on sample: $sample"
	read1=$diamond_processed/${sample}_1.m8
	read2=$diamond_processed/${sample}_2.m8

	## counting occurences of each gene in each read file and adding them together
	declare -i mst_count; mst_count=$(cat $read1 | grep -o "mst" | wc -l)
	mst_count+=$(cat $read2 | grep -o "mst" | wc -l)
	declare -i cse_count; cse_count=$(cat $read1 | grep -o "cse" | wc -l)
	cse_count+=$(cat $read2 | grep -o "cse" | wc -l)
	declare -i cbs_count; cbs_count=$(cat $read1 | grep -o "cbs" | wc -l)
	cbs_count+=$(cat $read2 | grep -o "cbs" | wc -l)
	declare -i cyd_count; cyd_count=$(cat $read1 | grep -o "cyd" | wc -l)
	cyd_count+=$(cat $read2 | grep -o "cyd" | wc -l)
	declare -i metC_count; metC_count=$(cat $read1 | grep -o "metC" | wc -l)
	metC_count+=$(cat $read2 | grep -o "metC" | wc -l)
	declare -i cysM_count; cysM_count=$(cat $read1 | grep -o "cysM" | wc -l)
	cysM_count+=$(cat $read2 | grep -o "cysM" | wc -l)
	declare -i cysK_count; cysK_count=$(cat $read1 | grep -o "cysK" | wc -l)
	cysK_count+=$(cat $read2 | grep -o "cysK" | wc -l)
	declare -i malY_count; malY_count=$(cat $read1 | grep -o "malY" | wc -l)
	malY_count+=$(cat $read2 | grep -o "malY" | wc -l)
	declare -i yhaO_count; yhaO_count=$(cat $read1 | grep -o "yhaO" | wc -l)
	yhaO_count+=$(cat $read2 | grep -o "yhaO" | wc -l)
	declare -i yhaM_count; yhaM_count=$(cat $read1 | grep -o "yhaM" | wc -l)
	yhaM_count+=$(cat $read2 | grep -o "yhaM" | wc -l)
	declare -i decR_count; decR_count=$(cat $read1 | grep -o "decR" | wc -l)
	decR_count+=$(cat $read2 | grep -o "decR" | wc -l)
	declare -i iscS_count; iscS_count=$(cat $read1 | grep -o "iscS" | wc -l)
	iscS_count+=$(cat $read2 | grep -o "iscS" | wc -l)
	declare -i mgl_count; mgl_count=$(cat $read1 | grep -o "mgl" | wc -l)
	mgl_count+=$(cat $read2 | grep -o "mgl" | wc -l)
	declare -i tnaA_count; tnaA_count=$(cat $read1 | grep -o "tnaA" | wc -l)
	tnaA_count+=$(cat $read2 | grep -o "tnaA" | wc -l)
	declare -i dsrA_count; dsrA_count=$(cat $read1 | grep -o "dsrA" | wc -l)
	dsrA_count+=$(cat $read2 | grep -o "dsrA" | wc -l)
	declare -i dsrB_count; dsrB_count=$(cat $read1 | grep -o "dsrB" | wc -l)
	dsrB_count+=$(cat $read2 | grep -o "dsrB" | wc -l)
	
	if [ $loop_count == 1 ]; then
		printf "%b" "Run\t" "mst\t" "cse\t" "cbs\t" "cyd\t" "metC\t" "cysM\t" \
		"cysK\t" "malY\t" "yhaO\t" "yhaM\t" "decR\t" "iscS\t" "mgl\t" "tnaA\t" \
		"dsrA\t" "dsrB\n" >> $diamond_processed/counts.tsv
	fi

	## writing counts to tsv
	printf "%b" "$sample\t" "$mst_count\t" "$cse_count\t" "$cbs_count\t" \
		"$cyd_count\t" "$metC_count\t" "$cysM_count\t" "$cysK_count\t" "$malY_count\t" \
		"$yhaO_count\t" "$yhaM_count\t" "$decR_count\t" "$iscS_count\t" "$mgl_count\t" \
		"$tnaA_count\t" "$dsrA_count\t" "$dsrB_count\n" \
		>> $diamond_processed/counts.tsv
	
	(( loop_count++ ))
done
