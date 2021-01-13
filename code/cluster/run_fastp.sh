#!/bin/bash

## description

echo "- initializing variables"
all_samples=$(ls $1/*.fastq | cut -d '/' -f3 | cut -d '_' -f 1 | uniq)
r1="_1"
r2="_2"
outdir=$2
mkdir -p results/fastp_out

echo "- trimming and merging read files with fastp"
for sample in $all_samples; do
	rf1="data/david-2014/$sample$r1.fastq"
	of1="results/fastp_out/$sample$r1.fastq"
	rf2="data/david-2014/$sample$r2.fastq"
	of2="results/fastp_out/$sample$r2.fastq"
	fastp -i $rf1 -o $of1 -I $rf2 -O $of2 \
		-m --merged_out $outdir/$sample.fastq --overlap_len_require 4 \
		--thread 16 -V
done

# echo "- running NGmerge"
# for sample in $all_samples; do
# 	readfile1="data/david-2014/$sample$r1.fastq"
# 	readfile2="data/david-2014/$sample$r2.fastq"
# 	external/NGmerge/NGmerge -1 $readfile1 -2 $readfile2 -o $outdir/$sample.fastq -q 39 -y
# done
