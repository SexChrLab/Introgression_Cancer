import os

# Env: introgression

configfile: "calc.perc.and.avg.introgression.analysis.genes.EUR.config.json"

rule all:
    input:
        expand(os.path.join(config["scratchdir"], "intersect_files_genes/EUR/chr{chrm}.genes.haplotype.bps.tsv"), chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results_new/genes/EUR/chr{chrm}.genes.haplotype.bps.introgression.per.hap.txt"), chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "results_new/genes/EUR/chr{chrm}.genes.haplotype.bps.avg.introgression.txt"), chrm=config["chromosomes"])


rule bedtoolsint:
    input:
        bedfn = os.path.join(config["coordsdir"], "chr{chrm}.human.GRCh37.GenesAndGenePredictions.NCBIRefSeq.genes.sorted.bed"),
        llfn = os.path.join(config["introdir"], config["introfn"])
    output:
        os.path.join(config["scratchdir"], "intersect_files_genes/EUR/chr{chrm}.genes.haplotype.bps.tsv")
    shell:
        "bedtools intersect -a {input.bedfn} -b {input.llfn} -wo > {output}"

rule runScript:
    input:
        coorsfn = os.path.join(config["coordsdir"], "chr{chrm}.human.GRCh37.GenesAndGenePredictions.NCBIRefSeq.genes.sorted.bed"),
        bedtfn = os.path.join(config["scratchdir"], "intersect_files_genes/EUR/chr{chrm}.genes.haplotype.bps.tsv")
    params:
        scrptdir = config["scrptdir"]
    output:
        out1 = os.path.join(config["scratchdir"], "results_new/genes/EUR/chr{chrm}.genes.haplotype.bps.introgression.per.hap.txt"),
        out2 = os.path.join(config["scratchdir"], "results_new/genes/EUR/chr{chrm}.genes.haplotype.bps.avg.introgression.txt")
    shell:
        "python3 {params.scrptdir}calc.perc.and.avg.introgression.genes.py --coords {input.coorsfn} --bedtfn {input.bedtfn} --out1 {output.out1} --out2 {output.out2} || echo 0"


