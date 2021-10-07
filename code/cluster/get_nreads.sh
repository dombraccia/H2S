#!/bin/bash
#SBATCH -J nreads # Job name 
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH -o logs/nreads_%j.o # output file
#SBATCH -e logs/nreads_%j.e # error file
#SBATCH --mail-type=fail,end # Get email for begin, end, and fail
#SBATCH --qos=large
#SBATCH --time=6:00:00
#SBATCH --mem=120gb
#SBATCH --ntasks=16

## getting nreads from each sample of metatranscriptomics data

datadir=$1
dataset=$(echo $1 | cut -d '/' -f 2)

rm -f data/$dataset/nreads.txt

runs=$(ls $datadir | grep -v '_2')

for run in $runs; do
	nreads=$(grep '^@' $datadir/$run | wc -l)
	echo "$run" "$nreads" >> data/$dataset/nreads.txt
done