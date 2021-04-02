#!/bin/bash
#SBATCH --job-name=call.tcga.exome # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 2
##SBATCH --mem-per-cpu=16000
#SBATCH -t 96:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

source activate varCalling

cd /scratch/amtarave/introgression_pilot/elife/TCGA_exome_processing/

#snakemake -s var.calling.snakefile --rerun-incomplete -j 40 --cluster "sbatch -n 2 --mem-per-cpu=20000 -t 96:00:00 --mail-type=END,FAIL --mail-user=amtarave@asu.edu"
snakemake -s var.calling.snakefile --rerun-incomplete -j 2 --cluster "sbatch -n 1 -t 96:00:00 --mail-type=END,FAIL --mail-user=amtarave@asu.edu"

mv slurm* logs/

source deactivate
