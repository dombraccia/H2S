#!/bin/bash

##### ========================= DESCRIPTION ============================= #####
# A script for downloading transcriptomic data from the Lawrence A. David 2014
# publication. 
# 
# Can be generalized in the future if need be.
##### =================================================================== #####

echo "- initializing variables"
sampleID=$(cat $1)
echo "-- TEST: sampleID = $sampleID"
dummy_outfile=$2

echo "- running fasterq-dump"
external/sratoolkit.2.10.9-ubuntu64/bin/fasterq-dump $sampleID \
	-x \
	--temp data/from-SRA/tmp \
	-O data/from-SRA/

echo "- writing to dummy file"
echo "finished downloading raw transcriptomics data." >> $dummy_outfile
