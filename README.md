# Introgression_Cancer
Here, we investigate the evidence for a relationship between Neanderthal introgression and liver cancer susceptibility.


## Contents
1. Quantify the average proportion of Neanderthal introgression per gene
2. Germline TCGA exome processing
3. Compare between individuals with liver cancer patients and non-affected individuals
4. Analysis of somatic variation in TCGA liver cancer patients


### 1. Quantifying the average proportion of Neanderthal introgression

Using introgression maps from (Vernot et al. 2016) we calculated the average proportion of Neanderthal introgression for each gene and then compared the distribution of the average proportion of Neanderthal introgression across all genes and gene sets involved with cancer (tumor suppressor genes, oncogenes, COSMIC gene list, and DNA repair genes).

Link to downloaded introgression maps from (Vernot et al. 2016) (`introgressed_haplotypes.tar.gz`): https://drive.google.com/drive/folders/0B9Pc7_zItMCVWUp6bWtXc2xJVkk

Scripts for this step can be found in `analyses/01_avg_prop_Nean/`.

Gene lists used in this analyses can be found here: `gene_lists/`.

### 2. Germline TCGA exome processing
Germline bams from 411 TCGA liver cancer patients were previously aligned using HISAT2 to a custom sex-specific version of GRCh38. See: (Natri, Wilson, and Buetow 2019) for more details.

Here, using these alignment files, we called and joint genotyped all samples across chromosomes 1-22.

Scripts for this step can be found in `analyses/02_TCGA_exome_processing/`.


### 3. Compare the proportion of samples with Neanderthal introgression between individuals with liver cancer and non-affected individuals

For this step, we used tag SNPs (Vernot et al. 2016) downloaded from: https://drive.google.com/drive/folders/0B9Pc7_zItMCVM05rUmhDc0hkWmc to identify whether samples had evidence for Neanderthal introgression at a given gene for both the TCGA liver cancer germline exome data (processed in 2) and non-affected individuals of matched ancestries (East Asian and European) from the 1000 Genomes Phase 3 resource. 1000 Genomes variants can be downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/.

Scripts for this step can be found in `analyses/03_tcga_vs_1000g/`.


### 4. Analysis of somatic variation in TCGA liver cancer patients
For this analysis, we used results generated from step 3. For each sample and for each gene, we identified whether there was evidence of Neanderthal introgression using tag SNPs (see step 3 and `calc.prop.N.genes.per.sample.py` for more details). We then intersect these results with somatic variant call files (https://portal.gdc.cancer.gov/) for the corresponding sample. We then can get counts, for each gene, of the number of samples with:

1. Neanderthal introgression and somatic mutations
2. Neanderthal introgression and no somatic mutations
3. No Neanderthal introgression and somatic mutations
4. No Neanderthal introgression and no somatic mutations

Then with this information, for each gene, we performed a Fisher’s exact test to test whether the proportion of samples with somatic mutations is different between liver cancer patients with Neanderthal introgression and liver cancer patients without Neanderthal introgression.

Scripts for this step can be found in `analyses/04_TCGA_somatic/`.
