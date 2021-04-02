#!/usr/bin/env python3
import csv
from optparse import  OptionParser

"""
# BACKGROUND #
# I want to go through tag snp stes and determine, for each TCGA sample whether 
# there are sites and individuals with a lot of homozygous Neanderthal sites 


# Input 1 - TCGA VCF. This VCF has already been filtered to keep the sites that 
# overlap the tag snp file. 
# /scratch/amtarave/introgression_pilot/elife/nead_sites_tcga/tcga_vcfs/EUR_tcga/chr1_EUR_tcga_tagSNPs_intersect_Nsites_alt_frq_sites.recode.vcf

# Input 2 - tag snp file. This has information for the Neanderthal and Denisovan
# allele for each tag snp.
cd /scratch/amtarave/introgression_pilot/elife/nead_sites_tcga/tag_snps_filtered/EUR_tcga/
head chr1_EUR_tcga_tagSNPs_intersect_Nsites_alt_frq_sites_intersect.bed 
1       3166567 3166568 C       T       1       0.00955 0.0     0.07061 0.15675 0.00298 0.07407 0.07157 T       C       1_3066509_3209504
1       3383285 3383286 G       A       1       0.02548 0.0     0.14409 0.19643 0.0497  0.12963 0.10429 A       G       1_3352783_3450985
1       3427297 3427298 G       A       1       0.02548 0.0     0.13977 0.13095 0.08748 0.09259 0.06442 A       G       1_3352783_3450985
1       3431554 3431555 C       T       1       0.02548 0.0     0.14121 0.12897 0.08648 0.11111 0.06442 T       C       1_3352783_3450985
1       3436879 3436880 G       A       1       0.02229 0.0     0.14986 0.12996 0.08648 0.09259 0.17178 A       G       1_3352783_3450985
1       3437013 3437014 G       A       1       0.02229 0.0     0.14986 0.13095 0.08648 0.09259 0.17178 A       G       1_3352783_3450985
1       3439256 3439257 A       G       1       0.02229 0.0     0.1513  0.09722 0.08648 0.05556 0.1636  G       A       1_3352783_3450985
1       3440444 3440445 G       A       1       0.02229 0.0     0.14986 0.13492 0.08549 0.09259 0.15235 A       G       1_3352783_3450985
1       3445285 3445286 G       C       1       0.02229 0.0     0.15274 0.13492 0.08748 0.18519 0.16053 C       C       1_3352783_3450985
1       3446249 3446250 G       A       1       0.02229 0.0     0.15274 0.13492 0.08748 0.09259 0.16155 A       G       1_3352783_3450985
# Column information
chr
start
stop
Ancestral allele
Derived allele
Anc/Der code
AFA allele frequency
AFR allele frequency
AMR allele frequency
EAS allele frequency
EUR allele frequency
PNG allele frequency
SAS allele frequency
Neanderthal base
Denisova base
haplotype tag


# output: list of sites where no sample has a Neanderthal allele
"""


################################################################################
USAGE = """
python3 get.tcga.neand.homozygous.sites.py	--vcf <path and name to vcf >
								--tagbed <path and name to tag snp bed file >
								--out <path and name for output file>

vcf == path and name to fist pass filtered vcf
tagbed == path and name to tag snp bed file
out == path and name for  output file
"""

parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcf', help = 'path and name to vcf')
parser.add_option('--tagbed',dest='tagbed', help = 'path and name to tag snp bed file')
parser.add_option('--out',dest='out', help = 'path and stem name for  output file')

(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.vcf is None:
	parser.error('path and name to vcf not given')
if options.tagbed is None:
	parser.error('path and name to tag snp bed file not give')
if options.out is None:
	parser.error('path and stem name for output file not give')
################################################################################
#outfile = open('chr1_EUR_tcga_tagSNPs_rm_list.txt', 'w')
outfile = open(options.out, 'w')

header = "chr\tpos\tmod/mod\tnean/nean\thet\tmissing\n"
outfile.write(header)

#with open('/scratch/amtarave/introgression_pilot/elife/nead_sites_tcga/tag_snps/introgressed_tag_snp_frequencies/all_tag_snps.EUR.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.noChr.bed', 'r') as hapdatin:
with open(options.tagbed, 'r') as hapdatin:
	hapdatfull = [line.strip().split() for line in hapdatin]
	#hapdat = [line.strip().split() for line in hapdatin]


#with open('/scratch/amtarave/introgression_pilot/elife/nead_sites_tcga/tcga_vcfs/EUR_tcga/chr1_EUR_tcga_tagSNPs_intersect.vcf', 'r') as vcf:
with open(options.vcf, 'r') as vcffile:
    #vcf_dat = [line.strip().split() for line in vcf]
    for line in vcffile:
    	if '#'in line:
    		continue
    	else:
    		line = line.strip().split()
    		# get Neanderthal allele info from tagbed file
    		hapdatsub = [ln for ln in hapdatfull if ln[0] == line[0] and ln[2] == line[1]] # matching chr and pos in tagbed and vcf
    		#print(line)
    		#print(hapdatsub)
    		# for every site get number of homoz mod, het mod and N, and homoz N, and no data
    		nallele = hapdatsub[0][13] # neanderthal allele at this site
    		ref = line[3]
	    	alt = line[4]
    		# is the neanderthal allele the ref or alt allele in vcf and if so get homoz N genotype
    		if nallele == ref:
    			nhomgeno = "0/0"
    		elif nallele == alt:
    			nhomgeno = "1/1"
    		#print("nhomgeno", nhomgeno)
    		# Now that we know what the Neanderthal homozygous genotype is we can
    		# iterate through each sample and ask whether it matches the Neanderthal 
    		# homozygous genotype, if not is it heterozygous, eles its modern/modern
    		homozNeand = 0
    		homozMod = 0
    		het = 0
    		missing = 0
    		for n in line[9:]:
    			gt = n.split(":")
    			gt = gt[0]
    			#print(gt)
    			if gt == nhomgeno:
    				homozNeand += 1
    				#print(gt)
    			elif gt == "./.":
    				missing += 1
    			elif gt == "0/1" or gt == "1/0":
    				het += 1
    			else:
    				homozMod += 1
    		#print("chr", line[0], "pos:", line[1], "Mod/Mod:", homozMod, "Nean/Nean:", homozNeand, "Het:", het, "Missing:", missing)
    		#print(het)
    		outfile.write('%s\t%s\t%i\t%i\t%i\t%i\n' % (str(line[0]), str(line[1]), homozMod, homozNeand, het, missing))
    		"""if len(hapdatsub) == 0: # for some reason there are sites in vcf not in tag bed...
    			# continue
    			# want to remove this sites probably
    			print(line)
    		else:
	    		nallele = hapdatsub[0][13] # neanderthal allele at this site
	    		#print(nallele)
	    		ref = line[3]
	    		alt = line[4]
	    		# is the neanderthal allele the ref or alt allele in vcf
	    		if nallele == ref:
	    			nallele = "0"
	    		elif nallele == alt:
	    			nallele = "1"
	    		#print(nallele)
	    		# iterate through each genotype, mark whether a sample has a neand allele
	    		# if by the end of the iteration no counts have been added to nallelecount
	    		# then it gets marked as a site to remove
	    		nallelecount = 0
	    		for n in line[9:]:
	    			gt = n.split(":")
	    			if gt[0] == "./.": # no data for tag snp for this sample
	    				continue
	    			else: # this sample has data
	    				gt = gt[0].split("/")
	    				#print(gt)
	    				gt1 = gt[0]
	    				gt2 = gt[1]
	    				#print(gt1, gt2)
	    				if gt1 == nallele or gt2 == nallele:
	    					nallelecount += 1
	    		if nallelecount == 0:
	    			#print(nallele)
	    			#print("no neanderthal allele: ", line)
	    			outfile.write('%s\t%s\n' % (str(line[0]), str(line[1])))
	    			#print(line[0], line[1])"""
outfile.close()
