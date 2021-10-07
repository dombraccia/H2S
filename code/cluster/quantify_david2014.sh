#!/bin/bash
#SBATCH -J salmon # Job name 
#SBATCH -o logs/salmon_%j_%a.o # output file
#SBATCH -e logs/salmon_%j_%a.e # error file
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH --mail-type=fail,end # Get email for end, and fail
#SBATCH --time=1:00:00
#SBATCH --qos=long
#SBATCH --mem=12gb
#SBATCH --array=1-59

# DESCRIPTION: script to run salmon quant on data from diet intervention study (david 2014)

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

## making output directory for function if does not already exist
mkdir -p results/salmon_out/$dataset/$function/

runs=$(ls data/david2014/rnaseq | cut -d '_' -f 1 | uniq)
run=$(echo $runs | cut -d ' ' -f ${SLURM_ARRAY_TASK_ID})
# run=$(echo $runs | head -${SLURM_ARRAY_TASK_ID} | tail -1)
# run=$(ls data/$dataset/rnaseq | head -${SLURM_ARRAY_TASK_ID} | tail -1 | cut -d '_' -f 1 | uniq)


if [ -f "data/hpfs/rnaseq/$run" ]; then
	echo "run $run contains single-end reads"
    run_name=$(echo $run | cut -d '.' -f 1)
	time salmon quant -i $gene_seqs_index -l A \
         -r data/$dataset/rnaseq/$run \
         -p 16 --validateMappings -o results/salmon_out/$dataset/$function/$run_name
else
	echo "run $run contains paired-end reads"
	time salmon quant -i $gene_seqs_index -l A \
         -1 data/$dataset/rnaseq/"$run"_1.fastq \
         -2 data/$dataset/rnaseq/"$run"_2.fastq \
         -p 16 --validateMappings -o results/salmon_out/$dataset/$function/$run
fi

echo "done with quantification step" > $dummy_out
