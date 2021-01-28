#!/bin/bash

files=*.fastq
for file in $files; do
	echo "- on file: $file"
	current_length=$(grep '^@' $file | wc -l)
	printf "%b" "$file\t" "$current_length\n" >> tmp_read_counts.tsv
done

rm sample_read_counts.tsv
touch sample_read_counts.tsv
samples=$(cut -f1 tmp_read_counts.tsv | cut -d '_' -f1 | uniq)
for sample in $samples; do
	echo "- on sample $sample"
	sample_count=$(grep "$sample" tmp_read_counts.tsv | cut -f2 | paste -sd+ - | bc)
	echo "- current sample count: $sample_count"
	printf "%b" "$sample\t" "$sample_count\n" >> sample_read_counts.tsv
done