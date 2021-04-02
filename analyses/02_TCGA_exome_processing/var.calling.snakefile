import os

# Environment: varCalling

configfile: "var.calling.config.json"


rule all:
    input:
        expand(os.path.join(config["output_path"], "gvcfs/{samples}.g.vcf.gz"), samples=config["all"]),
        expand(os.path.join(config["output_path"], "vcfs/{chrs}.tcga.exome.normal.GRCh38.gatk.called.raw.vcf.gz"), chrs=config["all_chrs"]),
        expand(os.path.join(config["output_path"], "vcfs/{chrs}.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.vcf.gz"), chrs=config["all_chrs"]),
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.vcf.gz"),
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.vcf.gz"),
        os.path.join(config["output_path"], "stats/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.bcftools.stats.PSC.txt"),
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.vcf.gz"),
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1.vcf.gz"),
        #os.path.join(config["output_path"], "stats/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.bcftools.stats.PSC.txt"),
        os.path.join(config["output_path"], "stats/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1.bcftools.stats.PSC.txt"),
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19.vcf"),
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19_rejected.vcf"),
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19.vcf.gz"),
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19.vcf.gz.tbi")


#-------------------------------------------------------------------------------
# Main step 1 #
# Generate gvcfs for all samples
#-------------------------------------------------------------------------------
rule make_gvcfs:
    input:
        bam = os.path.join(config["bam_path"], "{samples}_sorted_add_mark.bam")
    params:
        ref1 = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
        ref2 = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"]
        #chrms = "{chrms}"
    output:
        gvcf = os.path.join(config["output_path"], "gvcfs/{samples}.g.vcf.gz")
    run:
        if "-XX" in wildcards.samples or "_XX" in wildcards.samples:
            shell("gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref1} -I {input.bam} -O {output.gvcf} -ERC GVCF ")
        if "-XY" in wildcards.samples or "_XY" in wildcards.samples:
            shell("gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref2} -I {input.bam} -O {output.gvcf} -ERC GVCF ")


#-------------------------------------------------------------------------------
# Main step 2 #
# Combing gvcfs and joint genotype
#-------------------------------------------------------------------------------
rule combine_gvcfs_all_dbImport:
    input:
        gvcf =  expand(os.path.join(config["output_path"], "gvcfs/{samples}.g.vcf.gz"), samples=config["all"])
    params:
         ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"], # I am using the Ref_GRCh38_Y_PARsMasked ref becasue I do not have the default ref. This shouldn't matter because I am only calling variants on the autosomes.
         gvcfs = expand(("-V " + config["output_path"] + "gvcfs/{samples}.g.vcf.gz"), samples=config["all"]),
         chrms = "{chrs}",
         #db = os.path.join(config["output_path"], "my_database_all_{chrs}"),
         db = "my_database_all_{chrs}",
         ploidy = 2
    output:
        os.path.join(config["output_path"], "vcfs/{chrs}.tcga.exome.normal.GRCh38.gatk.called.raw.vcf.gz")
    shell:
        """
        gatk --java-options "-Xmx50g -Xms4g" GenomicsDBImport -L {params.chrms} {params.gvcfs} --genomicsdb-workspace-path {params.db};
        gatk --java-options "-Xmx50g" GenotypeGVCFs -R {params.ref} -V gendb://{params.db} -L {params.chrms} -ploidy {params.ploidy} -O {output} 
        """


#-------------------------------------------------------------------------------
# Main step 3 #
# Filter variants. Keep biallelic snps for downstream analysis
#-------------------------------------------------------------------------------
rule select_biallelic_all_diploid:
    input:
        os.path.join(config["output_path"], "vcfs/{chrs}.tcga.exome.normal.GRCh38.gatk.called.raw.vcf.gz")
    output:
        os.path.join(config["output_path"], "vcfs/{chrs}.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrs}"
    shell:
        """gatk --java-options "-Xmx10g" SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-O {output} """
        """--select-type-to-include SNP """
        """--restrict-alleles-to BIALLELIC """


#-------------------------------------------------------------------------------
# Main step 4 #
# Combine variants into one file
#-------------------------------------------------------------------------------
rule combineVariants:
    input:
        expand(os.path.join(config["output_path"], "vcfs/{chrs}.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.vcf.gz"), chrs=config["all_chrs"])
    params:
         vcfs = expand(("I=" + config["output_path"] + "vcfs/{chrs}.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.vcf.gz"), chrs=config["all_chrs"])
    output:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.vcf.gz")
    shell:
        """
        picard MergeVcfs {params.vcfs} O={output}
        """


#-------------------------------------------------------------------------------
# Main step 5 #
# Change sample names in vcf to match entity_submitter_id
#-------------------------------------------------------------------------------
rule reheader:
    input:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.vcf.gz")
    params:
        samples = config["newsamplelist"]
    output:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.vcf.gz")
    shell:
        """
        bcftools reheader -s {params.samples} {input} -o {output}
        """


#-------------------------------------------------------------------------------
# Main step 6 #
# Run bcftools stats on vcf to get per sample counts
#-------------------------------------------------------------------------------
rule statsBiallSNPs:
    input:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.vcf.gz")
    output:
        os.path.join(config["output_path"], "stats/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """


#-------------------------------------------------------------------------------
# Main step 7 #
# Perform some initial QC filtering and get per sample counts
#-------------------------------------------------------------------------------
rule filterVCFQUAL:
    input:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.vcf.gz")
    params:
        QUAL = 30
    output:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.vcf.gz")
    shell:
        """
        bcftools view  -i '%QUAL>={params.QUAL}' {input} -Oz -o {output}
        """

rule filterVCFFMISSING:
    input:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.vcf.gz")
    params:
        F_MISSING = 0.1
    output:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1.vcf.gz")
    shell:
        """
        bcftools view  -i 'F_MISSING<={params.F_MISSING}' {input} -Oz -o {output}
        """

'''
rule statsFilteredSNPs1:
	input:
		os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.vcf.gz")
	output:
		os.path.join(config["output_path"], "stats/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.bcftools.stats.PSC.txt")
	shell:
		"""
		bcftools stats -s - {input} | grep PSC > {output}
		"""
'''

rule statsFilteredSNPs2:
	input:
		os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1.vcf.gz")
	output:
		os.path.join(config["output_path"], "stats/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1.bcftools.stats.PSC.txt")
	shell:
		"""
		bcftools stats -s - {input} | grep PSC > {output}
		"""


#-------------------------------------------------------------------------------
# Main step 8 #
# Liftover cleaned vcf to hg19
#-------------------------------------------------------------------------------
rule Liftover:
    input:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1.vcf.gz")
    output:
        lvars = os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19.vcf"),
        rvars = os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19_rejected.vcf")
    params:
        chain = config["chain"],
        ref = config["liftover_ref"]
    shell:
        """
        java -Xmx12G -XX:-UseGCOverheadLimit -jar /home/amtarave/packages/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar LiftoverVcf --MAX_RECORDS_IN_RAM 100000 -I {input} -O {output.lvars} -C {params.chain} --REJECT {output.rvars} -R {params.ref}
        """

#gatk --java-options "-Xmx12G -XX:-UseGCOverheadLimit" LiftoverVcf --MAX_RECORDS_IN_RAM 100000 -I {input} -O {output.lvars} -C {params.chain} --REJECT {output.rvars} -R {params.ref}

rule zipandindex:
    input:
        os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19.vcf")
    output:
        vcf = os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19.vcf.gz"),
        tbi = os.path.join(config["output_path"], "vcfs/ALL.chrs.tcga.exome.normal.GRCh38.gatk.called.biallelic.snps.samples.renamed.QUAL30.FMISSING0.1_liftover_hg38ToHg19.vcf.gz.tbi")
    shell:
        """
        bgzip -c {input} > {output.vcf}; 
        tabix -p vcf {output.vcf}
        """