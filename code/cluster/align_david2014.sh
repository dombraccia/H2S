#!/bin/bash

## description

echo "- initializing variables"
db=$1
outdir=$2
echo "-- TEST: outdir = $outdir"
mkdir -p $outdir
readfiles=data/david-2014/*.fastq

echo "- running DIAMOND on all readfiles"
for file in $readfiles; do
	echo "-- on file: $file"
	outfile=$(echo $file | cut -d '/' -f 3 | cut -d '.' -f 1)
	echo "--- writing to output file: $outdir/$outfile.m8"
	external/diamond blastx --db $db -q $file -o $outdir/$outfile.m8 \
		-b42.0 -c1 # better performance (need node w 500gb of memory)
done
