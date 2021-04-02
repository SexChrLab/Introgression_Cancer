#!/bin/bash
#SBATCH --job-name=SAS.genes.perc.and.avg.introgression.analysis # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

conda activate introgression
#source activate introgression

cd /scratch/amtarave/introgression_pilot/elife/testing_scripts/

snakemake -s calc.perc.and.avg.introgression.analysis.genes.SAS.snakefile -j 23 --cluster "sbatch -n 1 -t 96:00:00 --mem-per-cpu=20000 --mail-type=END,FAIL --mail-user=amtarave@asu.edu"

mv slurm* ../logs/

source deactivate introgression
