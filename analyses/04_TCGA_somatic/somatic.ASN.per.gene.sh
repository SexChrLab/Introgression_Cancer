#!/bin/bash
#SBATCH --job-name=somatic.ASN.per.gene # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -q wildfire
#SBATCH -p wildfire
#SBATCH -t 48:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

conda activate introgression

cd /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/scripts

cat /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/new/file.ind.workflow.names.ASN.mutect2.txt | cut -f2 | while read line;
do
    echo "${line}"
    # Get somatic VCF file name that corresponds to the sample we are analyzing
    # line is the sample name from TCGA exome germline vcf
    # filename is the name of the somatic vcf that corresponds to the TCGA sample
    filename=`grep -w ${line} /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/new/file.ind.workflow.names.ASN.mutect2.txt | cut -f1`

    # Lift over somatic VCF for the sample from GRCh38 to hg19
    java -Xmx12G -XX:-UseGCOverheadLimit -jar /home/amtarave/packages/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar LiftoverVcf --MAX_RECORDS_IN_RAM 100000 -I /data/CEM/wilsonlab/temp/somatic_mutations/${filename}/${filename}.vep.vcf -O /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/00_new_tcga_qc/liftoverSomVCFs/${filename}.liftover.hg38tohg19.vcf -C /home/amtarave/projects/TCGA_exome_popInf/chainfile/hg38ToHg19.over.chain.gz --REJECT /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/00_new_tcga_qc/liftoverSomVCFs/${filename}.rejected.variants.vcf -R /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/hg19.fa

    # Previously, we identified for each sample, whether they had evidence of Neanderthal
    # introgression for a given gene. This information was generated using: calc.prop.N.genes.per.sample.py
    # Here we are extracting the per gene results for each sample so that we can
    # determine if for a given gene if a sample has somatic mutations and introgression
    # at that gene.
    grep ${line} /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/tcga_vs_1000g/results/prop_Nean_tcga/ASN_tcga/All.chr_ASN_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_sample_per_gene.txt > /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/results/ASN_tcga/${line}.per.gene.results.txt
    awk '{ print "chr"$0 }' /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/results/ASN_tcga/${line}.per.gene.results.txt > /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/results/ASN_tcga/${line}.per.gene.results.chr.txt

    # intersect the per gene per sample results file with the somatic vcf. This
    # will give us an output of where the somatic mutations are relative to the
    # genes we have data for.
    bedtools intersect -a /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/results/ASN_tcga/${line}.per.gene.results.chr.txt -b /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/new/liftoverSomVCFs/${filename}.liftover.hg38tohg19.vcf > /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/results/ASN_tcga/${line}.per.gene.results.overlap.all.somatic.vcfs.txt
done
# There was an issue with one sample (TCGA-DD-AACA-10A-01D-A40U-10 and ). They
# had multiple somatic vcfs, and we only want to analyze the somatic vcfs from
# primary tumor, so we fix that here outside of the loop.
# Remove this sample from output
grep -v "TCGA-DD-AACA-10A-01D-A40U-10" ASN_gene_counts_results.txt > ASN_gene_counts_results.edit2.txt

# Process primary tumor only
# 3204b5e2-62b2-422e-b35a-88af89bbb83c Recurrent Tumor
# cf1b22ac-9ec0-4da6-bbe4-8ab1b585e2b2  Recurrent Tumor
# e8c0aafb-8a5e-4de3-a19b-76fb767fd064 Primary Tumor
# Lift over
java -Xmx12G -XX:-UseGCOverheadLimit -jar /home/amtarave/packages/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar LiftoverVcf --MAX_RECORDS_IN_RAM 100000 -I /data/CEM/wilsonlab/temp/somatic_mutations/e8c0aafb-8a5e-4de3-a19b-76fb767fd064/e8c0aafb-8a5e-4de3-a19b-76fb767fd064.vep.vcf -O /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/new/liftoverSomVCFs/e8c0aafb-8a5e-4de3-a19b-76fb767fd064.liftover.hg38tohg19.vcf -C /home/amtarave/projects/TCGA_exome_popInf/chainfile/hg38ToHg19.over.chain.gz --REJECT /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/new/liftoverSomVCFs/e8c0aafb-8a5e-4de3-a19b-76fb767fd064.rejected.variants.vcf -R /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/hg19.fa
# Intersect per gene introgression results with somatic vcf
bedtools intersect -a /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/results/ASN_tcga/TCGA-DD-AACA-10A-01D-A40U-10.per.gene.results.chr.txt -b /scratch/amtarave/introgression_pilot/elife/somatic_sites_analysis/new/liftoverSomVCFs/e8c0aafb-8a5e-4de3-a19b-76fb767fd064.liftover.hg38tohg19.vcf > /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/results/ASN_tcga/TCGA-DD-AACA-10A-01D-A40U-10.per.gene.results.overlap.all.somatic.vcfs.txt

# We now have for each sample a file with information on which genes have somatic
# mutations and if that gene has evidence or no evidence for Neanderthal
# introgression.

# We next want to iterate through each gene and count the total number of
# samples that have:
# 1) Neanderthal introgression and somatic mutations
# 2) Neanderthal introgression and no somatic mutations
# 3) No Neanderthal introgression and somatic mutations
# 4) No Neanderthal introgression and no somatic mutations

# Get list of genes present in output (genes to analyze)
cd /scratch/amtarave/introgression_pilot/elife/introg_analyses_TCGA_exome/somatic/results/ASN_tcga
cut -f4 *.per.gene.results.chr.txt | sort | uniq > ASN.genes.txt

echo -e "gene\tsamples_non_intro_no_somatic\tsamples_intro_no_somatic\tsamples_non_intro_somatic\tsamples_intro_somatic" > ASN_gene_sample_somatic_counts.txt

# For each gene:
cat ASN.genes.txt | while read line;
do
    echo "${line}"

    # 1. get the total number of samples with introgression and no Neanderthl
    # introgression for this gene
    totother=`grep -E $'\t'${line}$'\t' *.per.gene.results.chr.txt | cut -f4,6,10 | sort | uniq -c | grep other | wc -l`
    totneand=`grep -E $'\t'${line}$'\t' *.per.gene.results.chr.txt | cut -f4,6,10 | sort | uniq -c | grep neand | wc -l`

    # 2. get the number of samples with no introgression and somatic mutations and
    # samples with introgression and somatic mutations for this gene
    othersomatic=`grep -E $'\t'${line}$'\t' *.per.gene.results.overlap.all.somatic.vcfs.txt | cut -f4,6,10 | sort | uniq -c | grep other | wc -l`
    neandsomatic=`grep -E $'\t'${line}$'\t' *.per.gene.results.overlap.all.somatic.vcfs.txt | cut -f4,6,10 | sort | uniq -c | grep neand | wc -l`

    # 3. get the number of samples with no introgression and no somatic mutations
    # and samples with introgression and no somatic mutations for this gene
    othernosomatic=`expr ${totother} - ${othersomatic}`
    neandnosomatic=`expr ${totneand} - ${neandsomatic}`

    # output counts for this gene
    echo -e "${line}\t${othernosomatic}\t${neandnosomatic}\t${othersomatic}\t${neandsomatic}" >> ASN_gene_sample_somatic_counts.txt
done

# This output file (ASN_gene_sample_somatic_counts.txt) will be used for the
# Fisher's exact tests. These tests will be performed in R.
