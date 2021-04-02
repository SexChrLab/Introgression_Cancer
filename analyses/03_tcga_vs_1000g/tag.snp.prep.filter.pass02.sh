#!/bin/bash
#SBATCH --job-name=tag.snp.prep.filter.pass02 # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -q wildfire 
#SBATCH -p wildfire 
#SBATCH -t 48:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

conda activate introgression
#source activate introgression

cd /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/tcga_vs_1000g/scripts/

snakemake -s tag.snp.prep.filter.pass02.snakefile --rerun-incomplete -j 22 --cluster "sbatch -n 2 -t 24:00:00 -q wildfire -p wildfire --mem-per-cpu=8000 --mail-type=END,FAIL --mail-user=amtarave@asu.edu"

mv slurm* ../logs/

source deactivate introgression