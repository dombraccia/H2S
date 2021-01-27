#!/bin/bash

## first part: taking the top hit for each read and saving in processed dir
## second part: distilling the read mappings into a count table for plotting

echo "- initializing variables"
diamond_out=$1/*.m8
diamond_processed=$2
function=$(basename $1)
mkdir -p $2

echo "- writing to processed directory"
for file in $diamond_out; do
	current_file=$(basename $file)
	cat $file | sort -u -k1,1 > $diamond_processed/$current_file
done

echo "- distilling data to a count table"
samples=$(ls $diamond_processed | cut -d '_' -f1 | uniq)
# diamond_processed=$2/*.m8
declare -i loop_count; loop_count=1
for sample in $samples; do
	## initializing read files 
	read1=$diamond_processed/${sample}_1.m8
	read2=$diamond_processed/${sample}_2.m8

	## counting occurences of each gene in each read file and adding them together
	declare -i fwdB_count; fwdB_count=$(cat $read1 | grep -o "fwdB" | wc -l)
	fwdB_count+=$(cat $read2 | grep -o "fwdB" | wc -l)
	declare -i fwdF_count; fwdF_count=$(cat $read1 | grep -o "fwdF" | wc -l)
	fwdF_count+=$(cat $read2 | grep -o "fdwF" | wc -l)
	declare -i fwdA_count; fwdA_count=$(cat $read1 | grep -o "fwdA" | wc -l)
	fwdA_count+=$(cat $read2 | grep -o "fwdA" | wc -l)
	declare -i fwdE_count; fwdE_count=$(cat $read1 | grep -o "fwdE" | wc -l)
	fwdE_count+=$(cat $read2 | grep -o "fwdE" | wc -l)
	declare -i fwdD_count; fwdD_count=$(cat $read1 | grep -o "fwdD" | wc -l)
	fwdD_count+=$(cat $read2 | grep -o "fwdD" | wc -l)
	declare -i fwdG_count; fwdG_count=$(cat $read1 | grep -o "fwdG" | wc -l)
	fwdG_count+=$(cat $read2 | grep -o "fwdG" | wc -l)
	declare -i Ftr_count; Ftr_count=$(cat $read1 | grep -o "Ftr" | wc -l)
	Ftr_count+=$(cat $read2 | grep -o "Ftr" | wc -l)
	declare -i Hmd_count; Hmd_count=$(cat $read1 | grep -o "Hmd" | wc -l)
	Hmd_count+=$(cat $read2 | grep -o "Hmd" | wc -l)
	declare -i Mch_count; Mch_count=$(cat $read1 | grep -o "Mch" | wc -l)
	Mch_count+=$(cat $read2 | grep -o "Mch" | wc -l)
	declare -i mcrB_count; mcrB_count=$(cat $read1 | grep -o "mcrB" | wc -l)
	mcrB_count+=$(cat $read2 | grep -o "mcrB" | wc -l)
	declare -i mcrG_count; mcrG_count=$(cat $read1 | grep -o "mcrG" | wc -l)
	mcrG_count+=$(cat $read2 | grep -o "mcrG" | wc -l)
	declare -i mcrA_count; mcrA_count=$(cat $read1 | grep -o "mcrA" | wc -l)
	mcrA_count+=$(cat $read2 | grep -o "mcrA" | wc -l)
	declare -i Mer_count; Mer_count=$(cat $read1 | grep -o "Mer" | wc -l)
	Mer_count+=$(cat $read2 | grep -o "Mer" | wc -l)
	declare -i mtrA_count; mtrA_count=$(cat $read1 | grep -o "mtrA" | wc -l)
	mtrA_count+=$(cat $read2 | grep -o "mtrA" | wc -l)
	declare -i mtrD_count; mtrD_count=$(cat $read1 | grep -o "mtrD" | wc -l)
	mtrD_count+=$(cat $read2 | grep -o "mtrD" | wc -l)
	declare -i mtrE_count; mtrE_count=$(cat $read1 | grep -o "mtrE" | wc -l)
	mtrE_count+=$(cat $read2 | grep -o "mtrE" | wc -l)
	declare -i mtrB_count; mtrB_count=$(cat $read1 | grep -o "mtrB" | wc -l)
	mtrB_count+=$(cat $read2 | grep -o "mtrB" | wc -l)
	declare -i mtrC_count; mtrC_count=$(cat $read1 | grep -o "mtrC" | wc -l)
	mtrC_count+=$(cat $read2 | grep -o "mtrC" | wc -l)
	declare -i mtrF_count; mtrF_count=$(cat $read1 | grep -o "mtrF" | wc -l)
	mtrF_count+=$(cat $read2 | grep -o "mtrF" | wc -l)
	declare -i mtrG_count; mtrG_count=$(cat $read1 | grep -o "mtrG" | wc -l)
	mtrG_count+=$(cat $read2 | grep -o "mtrG" | wc -l)

	if [ $loop_count == 1 ]; then
		printf "%b" "Run\t" "fwdA\t" "fwdB\t" "fwdD\t" "fwdE\t" "fwdF\t" "fwdG\t" \
		"Ftr\t" "Hmd\t" "Mch\t" "mcrA\t" "mcrB\t" "mcrG\t" "Mer\t" "mtrA\t" "mtrB\t" \
		"mtrC\t" "mtrD\t" "mtrE\t" "mtrF\t" "mtrG\n" >> $diamond_processed/counts.tsv
	fi

	## writing counts to tsv
	printf "%b" "$sample\t" "$fwdA_count\t" "$fwdB_count\t" "$fwdD_count\t" "$fwdE_count\t" \
		"$fwdF_count\t" "$fwdG_count\t" "$Ftr_count\t" "$Hmd_count\t" "$Mch_count\t" "$mcrA_count\t" \
		"$mcrB_count\t" "$mcrG_count\t" "$Mer_count\t" "$mtrA_count\t" "$mtrB_count\t" "$mtrC_count\t" \
		"$mtrD_count\t" "$mtrE_count\t" "$mtrF_count\t" "$mtrG_count\n" \
		>> results/diamond_processed/$function/counts.tsv
	
	(( loop_count++ ))
done
