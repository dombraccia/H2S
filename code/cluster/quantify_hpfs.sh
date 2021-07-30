#!/bin/bash
#SBATCH -J salmon # Job name 
#SBATCH -o logs/salmon/salmon_%j_%a.o # output file
#SBATCH -e logs/salmon/salmon_%j_%a.e # error file
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH --mail-type=fail,end # Get email for end, and fail
#SBATCH --time=1:00:00
#SBATCH --qos=long
#SBATCH --mem=12gb
#SBATCH --array=439-677

# [description]

echo '- loading modules'
module purge
module load salmon

echo '- loading in stdin from snakefile'
gene_seqs_index=$1
dataset=$2
dummy_out=$3
function=$(cut -d '/' -f 3 $1 | cut -d '_' -f 1,2)

run=$(ls data/hpfs/rnaseq | head -${SLURM_ARRAY_TASK_ID} | tail -1 | cut -d '_' -f 1 | uniq)

if [ -f "data/hpfs/rnaseq/$run" ]; then
	echo "run $run contains single-end reads"
    run_name=$(echo $run | cut -d '.' -f 1)
	time salmon quant -i $gene_seqs_index -l A \
         -r data/hpfs/rnaseq/$run \
         -p 16 --validateMappings -o results/salmon_out/$dataset/$function/$run_name
else
	echo "run $run contains paired-end reads"
	time salmon quant -i $gene_seqs_index -l A \
         -1 data/hpfs/rnaseq/"$run"_1.fastq \
         -2 data/hpfs/rnaseq/"$run"_2.fastq \
         -p 16 --validateMappings -o results/salmon_out/$dataset/$function/$run
fi

echo "done with quantification step" > $dummy_out
