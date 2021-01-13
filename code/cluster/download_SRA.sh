#!/bin/bash

##### ========================= DESCRIPTION ============================= #####
# A script for downloading transcriptomic data from the Lawrence A. David 2014
# publication. 
# 
# Can be generalized in the future if need be.
##### =================================================================== #####

echo "- initializing variables"
sampleIDs=$(cut -d , -f 1 $1)
dummy_outfile=$2

echo "- running fasterq-dump"
for sample in $sampleIDs; do
	external/sratoolkit.2.10.9-ubuntu64/bin/fasterq-dump $sample -S \
		-x \
		--temp data/david-2014/tmp \
		-O data/david-2014/
done

echo "- writing to dummy file"
echo "finished downloading raw transcriptomics data." >> $dummy_outfile
