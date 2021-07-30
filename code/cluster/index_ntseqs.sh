#!/bin/bash
#SBATCH -J salmon # Job name 
#SBATCH -o salmon_%j_%a.o # output file
#SBATCH -e salmon_%j_%a.e # error file
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH --mail-type=fail,end # Get email for end, and fail
#SBATCH --time=1:00:00
#SBATCH --qos=large
#SBATCH --mem=128gb

echo "- loading in modules"
module purge
module load salmon

echo "- creating index of selected nt seqs"
time salmon index -t $1 -i $2
