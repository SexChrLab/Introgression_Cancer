#!/usr/bin/env python3
import csv
from optparse import  OptionParser

"""
# BACKGROUND #
# I want to go through all tag snp stes and determine, for each TCGA sample which 
# alleles are neadnerthal. If there is a tag snp site where no samples have a 
# Neanderthal allele, I also want to make note of that so I can further filter
# the TCGA data


# Input 1 - TCGA VCF. This VCF has already been filtered to keep the sites that 
# overlap the tag snp file. 
# /scratch/amtarave/introgression_pilot/elife/nead_sites_tcga/tcga_vcfs/EUR_tcga/chr1_EUR_tcga_tagSNPs_intersect.vcf

# Input 2 - tag snp file. This has information for the Neanderthal and Denisovan
# allele for each tag snp.
head all_tag_snps.EUR.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.noChr.bed
1       2270126 2270127 C       T       1       0.0     0.00198 0.0     0.00595 0.0     0.0     0.00102 T       C       1_2270126_2300081
1       2273653 2273654 T       A       1       0.0     0.0     0.0     0.00595 0.0     0.0     0.00102 A       T       1_2270126_2300081
1       2285111 2285112 G       A       1       0.0     0.0     0.0     0.00694 0.0     0.0     0.01329 A       G       1_2270126_2300081
1       2292649 2292650 A       C       1       0.0     0.0     0.0     0.00694 0.0     0.0     0.00102 C       A       1_2270126_2300081
1       2299058 2299059 C       T       1       0.0     0.0     0.0     0.00694 0.0     0.0     0.00102 T       C       1_2270126_2300081
1       2300080 2300081 C       T       1       0.0     0.0     0.0     0.00694 0.0     0.0     0.00102 T       C       1_2270126_2300081
1       2855996 2855997 G       A       1       0.0     0.0     0.00144 0.06448 0.0     0.07407 0.07055 A       G       1_2855996_2879424
1       2858726 2858727 G       A       1       0.0     0.0     0.00144 0.06448 0.0     0.07407 0.07055 A       G       1_2855996_2879424
1       2861780 2861781 G       A       1       0.0     0.0     0.00144 0.06448 0.0     0.07407 0.07055 G/A     G       1_2855996_2879424
1       2863783 2863784 G       A       1       0.0     0.0     0.00144 0.06448 0.0     0.03704 0.06953 A       G       1_2855996_2879424
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
python get.tcga.neand.allele.info.py	--vcf <path and name to vcf >
								--tagbed <path and name to tag snp bed file >
								--out <path and name for two output files>

vcf == path and name to fist pass filtered vcf
tagbed == path and name to tag snp bed file
out == path and name for two output files
"""

parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcf', help = 'path and name to vcf')
parser.add_option('--tagbed',dest='tagbed', help = 'path and name to tag snp bed file')
parser.add_option('--out',dest='out', help = 'path and stem name for two output files')

(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.vcf is None:
	parser.error('path and name to vcf not given')
if options.tagbed is None:
	parser.error('path and name to tag snp bed file not give')
if options.out is None:
	parser.error('path and stem name for two output files not give')
################################################################################
#outfile = open('chr1_EUR_tcga_tagSNPs_rm_list.txt', 'w')
outfile = open(options.out, 'w')


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
    		if len(hapdatsub) == 0: # for some reason there are sites in vcf not in tag bed...
    			continue
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
	    			#print(line[0], line[1])
outfile.close()
