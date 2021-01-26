#!/bin/bash
#SBATCH -J snake # Job name 
#SBATCH --mail-user=dbraccia@umd.edu # Email for job info
#SBATCH -o logs/snake_%j.o # output file
#SBATCH -e logs/snake_%j.e # error file
#SBATCH --mail-type=fail,end # Get email for begin, end, and fail
#SBATCH --qos=xlarge
#SBATCH --time=24:00:00
#SBATCH --mem=500gb

# SCRIPT FOR RUNNING SNAKEMAKE SCRIPTS ON CLUSTER

# $1: give desired output of snakefile
time snakemake -p --cores all $1
