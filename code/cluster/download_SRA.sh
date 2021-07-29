#!/bin/bash
#SBATCH -J snake # Job name 
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH -o logs/snake_%j.o # output file
#SBATCH -e logs/snake_%j.e # error file
#SBATCH --mail-type=fail,end # Get email for begin, end, and fail
#SBATCH --qos=throughput
#SBATCH --time=18:00:00
#SBATCH --mem=36gb
#SBATCH --array=1-756

##### ========================= DESCRIPTION ============================= #####
# A script for downloading transcriptomic data from the HPFS publication. 
# 
# Can be generalized in the future if need be.
##### =================================================================== #####

echo "- initializing variables"
# sampleIDs=$(grep ',rna,' $1 | cut -d , -f 1)
grep ',rna,' $1 > data/hpfs/SraRunTable_RNA.txt
dummy_outfile=$2

echo "- running fasterq-dump"
run=$(head -${SLURM_ARRAY_TASK_ID} data/hpfs/SraRunTable_RNA.txt | tail -1 | cut -d , -f 1)
# cut -d , -f 1 $run 
external/sratoolkit.2.10.9-ubuntu64/bin/fasterq-dump $run -S \
	-x \
	--temp data/hpfs/tmp \
	-O data/hpfs/rnaseq

echo "- writing to dummy file"
echo "finished downloading raw transcriptomics data." > $dummy_outfile
