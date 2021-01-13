#!/bin/bash

## description

echo "- initializing variables"
diamond_out=$1/*.m8
diamond_processed=$2
mkdir -p results/diamond_processed

echo "- writing to processed directory"
for file in $diamond_out; do
	current_file=$(basename $file)
	cat $file | sort -u -k1,1 > $diamond_processed/$current_file
done
