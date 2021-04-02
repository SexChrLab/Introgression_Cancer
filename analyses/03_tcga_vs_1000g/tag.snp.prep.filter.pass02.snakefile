import os

# Env: introgression

configfile: "tag.snp.prep.filter.config.json"

rule all:
    input:
        expand(os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga.vcf.gz"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect.vcf"), pop=config["population"], chrm=config["chromosomes"]),
        #expand(os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_rm_list.txt"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "tag_snps_filtered/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_intersect.bed"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "1000g_intro_haps/{pop}/chr{chrm}_{pop}_tcga_tagSNPs_intersect_intersect_1000g_introg_haps_dat.txt"), pop=config["population"], chrm=config["chromosomes"]),
#        expand(os.path.join(config["scratchdir"], "1000g_vcfs/{pop}/chr{chrm}_{pop}_1000g.vcf.gz"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "1000g_vcfs/{pop}/chr{chrm}_{pop}_1000g_tcga_tag_sites.recode.vcf"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/chr{chrm}_{pop}_1000g_tcga_tag_sites_mod_intro_alt_frq.txt"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/All.chr_{pop}_1000g_tcga_tag_sites_mod_intro_alt_frq.txt"), pop=config["population"]),
        expand(os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/All.chr_{pop}_tcga_tagSNPs_intersect_alt_frq_keep_list.txt"), pop=config["population"]),
        expand(os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/All.chr_{pop}_1000g_tcga_tag_sites_mod_intro_alt_frq_plots.pdf"), pop=config["population"]),
        expand(os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites.recode.vcf"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "tag_snps_filtered/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_intersect.bed"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_sample_per_gene.txt"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_gene.txt"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/All.chr_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_gene.txt"), pop=config["population"]),
        expand(os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/chr{chrm}_{pop}_1000g_introgression_per_sample_per_gene.txt"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/chr{chrm}_{pop}_1000g_introgression_per_gene.txt"), pop=config["population"], chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/All.chr_{pop}_1000g_introgression_per_gene.txt"), pop=config["population"])


# rule: extract samples from tcga vcfs (separated by chromosome)
rule extractSamplesTCGA:
    input:
        tcgavcf = os.path.join(config["tcgaVCFdir"], config["tcgaVCFfn"]),
        keeplist = os.path.join(config["scratchdir"], "lists/{pop}.tcga.txt")
    params:
        chrm = "{chrm}"
    output:
        tcgavcfout = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga.vcf.gz")
    shell:
        """
        vcftools --gzvcf {input.tcgavcf} --chr {params.chrm} --keep {input.keeplist} --recode --stdout | gzip -c > {output.tcgavcfout}
        """

# rule: Intersect population specific tag snp file with population specific vcf
rule intersectTagsVCFs:
    input:
        vcf = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga.vcf.gz"),
        tagbed = os.path.join(config["tagdir"], "all_tag_snps.{pop}.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.noChr.bed")
    output:
        vcfout = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect.vcf")
    shell:
        """
        bedtools intersect -a {input.vcf} -b {input.tagbed} -header > {output.vcfout}
        """

# rule: Intersect 1. Intersect tag snp bed file with sites in TCGA vcf. This 
# will subset the tag snp bed file to only keep sites that overlap the TCGA vcf
rule intersectTagSNPs:
    input:
        tagbed = os.path.join(config["tagdir"], "all_tag_snps.{pop}.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.noChr.bed"),
        vcf = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect.vcf")
    output:
        os.path.join(config["scratchdir"], "tag_snps_filtered/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_intersect.bed")
    shell:
        """
        bedtools intersect -a {input.tagbed} -b {input.vcf} > {output}
        """

# rule: Intersect 2. Intersect tag snp bed from last step with introgressed 
# haplotype data from 1000 genomes. Also subset columns from output
# Note: For ease of running this pipeline, I changed file name: 
# LL.callsetEAS.mr_0.99.full_data.bed to: LL.callsetASN.mr_0.99.full_data.bed
rule intersect1000gIntroHaps:
    input:
        introhaps = (config["introHapsDir"] + "LL.callset{pop}.mr_0.99.full_data.bed"),
        tags = os.path.join(config["scratchdir"], "tag_snps_filtered/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_intersect.bed")
    output:
        os.path.join(config["scratchdir"], "1000g_intro_haps/{pop}/chr{chrm}_{pop}_tcga_tagSNPs_intersect_intersect_1000g_introg_haps_dat.txt")
    shell:
        """
        bedtools intersect -a {input.introhaps} -b {input.tags} -wo | cut -f13,14,15,16,17,28 | uniq > {output}
        """

# rule: Extract sites from rmFilterTagsTCGA from 1000 genomes pop vcf
# We want to extract the tag snp sites from TCGA in the 1000 genomes pop samples
# This will make parsing this VCF faster in the next step.
rule extractTCGAsitesFrom1000g:
    input:
        vcf1 = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect.vcf"),
        vcf2 = os.path.join(config["1000gVCFdir"], "{pop}/chr{chrm}_{pop}_1000g.vcf.gz")
    params:
        stem = os.path.join(config["scratchdir"], "1000g_vcfs/{pop}/chr{chrm}_{pop}_1000g_tcga_tag_sites")
    output:
        os.path.join(config["scratchdir"], "1000g_vcfs/{pop}/chr{chrm}_{pop}_1000g_tcga_tag_sites.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input.vcf2} --positions {input.vcf1} --recode --out {params.stem}
        """

# rule: Get alt allele frq on modern haps and on introgressed haps
# Use script get.1000g.sample.allele.info.py
rule getAltFrqs:
    input:
        vcf = os.path.join(config["scratchdir"], "1000g_vcfs/{pop}/chr{chrm}_{pop}_1000g_tcga_tag_sites.recode.vcf"),
        haps = os.path.join(config["scratchdir"], "1000g_intro_haps/{pop}/chr{chrm}_{pop}_tcga_tagSNPs_intersect_intersect_1000g_introg_haps_dat.txt"),
        tagbed = os.path.join(config["scratchdir"], "tag_snps_filtered/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_intersect.bed")
    output:
        os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/chr{chrm}_{pop}_1000g_tcga_tag_sites_mod_intro_alt_frq.txt")
    shell:
        """
        python3 get.1000g.sample.allele.info.py --vcf {input.vcf} --haps {input.haps} --tagbed {input.tagbed} --out {output}
        """

# rule merge results
rule mergeResults:
    input:
        resultsf = expand(os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/chr{chrm}_{pop}_1000g_tcga_tag_sites_mod_intro_alt_frq.txt"), pop=config["population"], chrm=config["chromosomes"])
    params:
        rdir = os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/")
    output:
        os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/All.chr_{pop}_1000g_tcga_tag_sites_mod_intro_alt_frq.txt")
    shell:
        """
        cat {params.rdir}* > {output} 
        """

# rule: plot results and make keep lists
rule plotMakeLists:
    input:
        os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/All.chr_{pop}_1000g_tcga_tag_sites_mod_intro_alt_frq.txt")
    output:
        out1 = os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/All.chr_{pop}_tcga_tagSNPs_intersect_alt_frq_keep_list.txt"),
        out2 = os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/All.chr_{pop}_1000g_tcga_tag_sites_mod_intro_alt_frq_plots.pdf")
    shell:
        """
        Rscript alt.allele.frqs.mod.intro.haps.1000g.R {input} {output.out1} {output.out2}
        """

# rule: extract those sites in the TCGA vcfs from last rule 
rule extractAltFrqSitesTCGA:
    input:
        vcf = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect.vcf"),
        keep = os.path.join(config["scratchdir"], "results/1000g_alt_frqs/{pop}/All.chr_{pop}_tcga_tagSNPs_intersect_alt_frq_keep_list.txt")
    params:
        vcfstem = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites")
    output:
        os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites.recode.vcf")
    shell:
        """
        vcftools --vcf {input.vcf} --positions {input.keep} --recode --out {params.vcfstem}
        """


# rule: intersect tag snp file with this vcf to get a tagsnp file with only 
# sites in tcga vcf
rule intersectTagSNPswithTCGAfiltered:
    input:
        tagbed = os.path.join(config["tagdir"], "all_tag_snps.{pop}.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.noChr.bed"),
        vcf = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites.recode.vcf")
    output:
        os.path.join(config["scratchdir"], "tag_snps_filtered/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_intersect.bed")
    shell:
        """
        bedtools intersect -a {input.tagbed} -b {input.vcf} > {output}
        """

# rule: run calc.prop.N.genes.per.sample.py to get the proportion of neand 
# introgression per gene and introgression status per sample per gene
rule calcPropNTCGA:
    input:
        vcf = os.path.join(config["scratchdir"], "tcga_vcfs/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites.recode.vcf"),
        tagbed = os.path.join(config["scratchdir"], "tag_snps_filtered/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_intersect.bed"),
        coord = (config["ncbiGeneCoorDif"] + "chr{chrm}" + config["ncbiGeneCoorSuffix"])
    params:
        outstem = os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites")
    output:
        out1 = os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_sample_per_gene.txt"),
        out2 = os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_gene.txt")
    shell:
        """
        python3 calc.prop.N.genes.per.sample.py --vcf {input.vcf} --coord {input.coord} --tagbed {input.tagbed} --out {params.outstem}
        """

# rule merge results from last rule (calcPropNTCGA)
rule mergeResultsPropNTCGA:
    input:
        results = expand(os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_gene.txt"), pop=config["population"], chrm=config["chromosomes"])
    params:
        rdir = os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/")
    output:
        os.path.join(config["scratchdir"], "results/prop_Nean_tcga/{pop}_tcga/All.chr_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_gene.txt")
    shell:
        """
        cat {params.rdir}*_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_gene.txt > {output} 
        """

# rule: run calc.prop.N.genes.per.sample.1000g.py to get the proportion of neand 
# introgression per gene and introgression status per sample per gene
rule calcPropN1000g:
    input:
        vcf = os.path.join(config["scratchdir"], "1000g_vcfs/{pop}/chr{chrm}_{pop}_1000g_tcga_tag_sites.recode.vcf"),
        tagbed = os.path.join(config["scratchdir"], "tag_snps_filtered/{pop}_tcga/chr{chrm}_{pop}_tcga_tagSNPs_intersect_alt_frq_sites_intersect.bed"),
        coord = (config["ncbiGeneCoorDif"] + "chr{chrm}" + config["ncbiGeneCoorSuffix"])
    params:
        outstem = os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/chr{chrm}_{pop}_1000g")
    output:
        out1 = os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/chr{chrm}_{pop}_1000g_introgression_per_sample_per_gene.txt"),
        out2 = os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/chr{chrm}_{pop}_1000g_introgression_per_gene.txt")
    shell:
        """
        python3 calc.prop.N.genes.per.sample.1000g.py --vcf {input.vcf} --coord {input.coord} --tagbed {input.tagbed} --out {params.outstem}
        """

# rule merge results from last rule (calcPropN1000g)
rule mergeResultsPropN1000g:
    input:
        results = expand(os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/chr{chrm}_{pop}_1000g_introgression_per_gene.txt"), pop=config["population"], chrm=config["chromosomes"])
    params:
        rdir = os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/")
    output:
        os.path.join(config["scratchdir"], "results/1000g_neand_sites/{pop}/All.chr_{pop}_1000g_introgression_per_gene.txt")
    shell:
        """
        cat {params.rdir}*_1000g_introgression_per_gene.txt > {output} 
        """
