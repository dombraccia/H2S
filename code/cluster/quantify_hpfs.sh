#!/bin/bash
#SBATCH -J salmon # Job name 
#SBATCH -o logs/salmon/salmon_%j_%a.o # output file
#SBATCH -e logs/salmon/salmon_%j_%a.e # error file
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH --mail-type=fail,end # Get email for end, and fail
#SBATCH --time=1:00:00
#SBATCH --qos=long
#SBATCH --mem=12gb
#SBATCH --array=1-1044

# [description]

echo '- loading modules'
module purge
module load salmon

echo '- loading in stdin from snakefile'
gene_seqs_index=$1
echo "-- the index is: $gene_seqs_index"
dataset=$(echo $2 | cut -d '/' -f 2)
echo "-- the dataset being maped is: $dataset"
dummy_out=$3
echo "-- the specified dummy outfile is: $dummy_out"
function=$(echo $gene_seqs_index | cut -d '/' -f3 | cut -d '_' -f 1)
echo "-- the inputted function is: $function"

# remaking dummy_out
rm -f $dummy_out

# ensuring that particular directory for this function has been created
mkdir -p results/salmon_out/$dataset/$function/
echo "-- IMPORTANT, WRITING TO OUTFILES TO: results/salmon_out/$dataset/$function/"

# iterating through all of the rnaseq runs over many nodes in slurm cluster
run=$(ls data/hpfs/rnaseq | head -${SLURM_ARRAY_TASK_ID} | tail -1 | cut -d '_' -f 1 | uniq)

if [ -f "data/hpfs/rnaseq/$run" ]; then
	echo "run $run contains single-end reads"
    run_name=$(echo $run | cut -d '.' -f 1)
    if [ ! -d "results/salmon_out/$dataset/$function/$run_name" ]; then
	    time salmon quant -i $gene_seqs_index -l A \
             -r data/hpfs/rnaseq/$run \
             -p 16 --validateMappings -o results/salmon_out/$dataset/$function/$run_name
    fi
else
	echo "run $run contains paired-end reads"
    if [ ! -d "results/salmon_out/$dataset/$function/$run" ]; then
	    time salmon quant -i $gene_seqs_index -l A \
             -1 data/hpfs/rnaseq/"$run"_1.fastq \
             -2 data/hpfs/rnaseq/"$run"_2.fastq \
             -p 16 --validateMappings -o results/salmon_out/$dataset/$function/$run
    fi

fi

echo "done with quantification step" > $dummy_out
