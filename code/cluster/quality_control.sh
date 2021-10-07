#!/bin/bash
#SBATCH -J qc # Job name 
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH --mail-type=fail,end # Get email for begin, end, and fail
#SBATCH --qos=throughput
#SBATCH --time=18:00:00
#SBATCH --mem=36gb
#SBATCH --array=1-2088

# SBATCH options

# a script to perform quality control of fastq read files using fastp
# (downloaded to conda environment via `conda install -c bioconda fastp`)

rnaseq_dir=$1
echo "- the rnaseq dir is: $rnaseq_dir"
# single=$(ls data/hpfs/rnaseq/ | grep -v '_')
# paired=$(ls data/hpfs/rnaseq/ | grep '_' | cut -d '_' -f 1)

run=$(ls data/hpfs/rnaseq | head -${SLURM_ARRAY_TASK_ID} | tail -1 | cut -d '.' -f 1)

if [ `echo $run | grep -v '_' | wc -c` -gt 0 ]; then
	echo "- run $run contains single-end reads"
	fastp -i $rnaseq_dir/$run.fastq -o $rnaseq_dir/$run.trimmed.fastq
elif [ `echo $run | grep '_' | wc -c` -gt 0 ]; then
	run2=$(echo $run | cut -d '_' -f 1)
	echo "- run $run2 contains paired-end reads"
	fastp -i $rnaseq_dir/"$run2"_1.fastq -I $rnaseq_dir/"$run2"_2.fastq \
		  -o $rnaseq_dir/"$run2"_1.trimmed.fastq -O $rnaseq_dir/"$run2"_2.trimmed.fastq
fi

echo "done with quality control step" > data/hpfs/dummy/qc.txt

##### non-parallelized way of using fastp

# for file in $single
# do
# 	run=$(echo $file | cut -d '.' -f 1)
# 	echo "on single run: $run"
# 	fastp -i $rnaseq_dir/$file -o $rnaseq_dir/$run.trimmed.fastq
# done

# for run in $paired
# do
# 	echo "on paired run: $run"
# 	fastp -i $rnaseq_dir/"$run"_1.fastq -I $rnaseq_dir/"$run"_2.fastq \
# 		  -o $rnaseq_dir/"$run"_1.trimmed.fastq -O $rnaseq_dir/"$run"_2.trimmed.fastq \
# done
