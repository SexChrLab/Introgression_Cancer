#!/usr/bin/env python3
import csv
from optparse import  OptionParser
################################################################################
USAGE = """
python3 get.1000g.sample.allele.info.py	--vcf <path and name to vcf >
								--haps <path and name to introgressed haplotype data >
								--tagbed <path and name to subsetted tag snp bed file, see step 4 in readme >
								--out <path and name for output file>

vcf == path and name to 1000 genome filtered vcf
haps == path and name to introgressed haplotype data
out == path and name for output file
"""

parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcf', help = 'path and name to vcf')
parser.add_option('--haps',dest='haps', help = 'path and name to introgressed haplotype data')
parser.add_option('--tagbed',dest='tagbed', help = 'path and name to subsetted tag snp bed file, see step 4 in readme')
parser.add_option('--out',dest='out', help = 'path and stem name for output file')

(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.vcf is None:
	parser.error('path and name to vcf not given')
if options.haps is None:
	parser.error('path and name to introgressed haplotype data not give')
if options.tagbed is None:
	parser.error('path and name to subsetted tag snp bed file not give')
if options.out is None:
	parser.error('path and stem name for output file not give')
################################################################################
# in 1: EAS.sample.haps.tag.snps.dat.txt
# in 2: chr1_ASN_1000g.tcga.tag.sites.vcf.recode.vcf
outfile = open(options.out, 'w')

header = "chr\tstart_pos\tend_pos\talt_frq_mod\talt_frq_intro\n"
print(header)
outfile.write(header)
with open(options.vcf, 'r') as vcf:
		vcf_dat = []
		for line in vcf:
			if "##" in line:
				continue
			elif "#" in line:
				sampleheader = line.strip().split()
			else:
				line = line.strip().split()
				vcf_dat.append(line)

with open(options.haps, 'r') as hapdatin:
	hapdatfull = [line.strip().split() for line in hapdatin]

with open(options.tagbed, 'r') as coords:
	for coord in coords:
		coord = coord.strip().split()
		# subset hapdatfull to only keep sites in coord
		# then proceed as before
		hapdat = [ln for ln in hapdatfull if ln[2] == coord[0] and ln[3] == coord[1] and ln[4] == coord[2]]
		hapdat_uniq = dict() 
		for i in range(0,len(hapdat)):
			if hapdat[i][1] not in hapdat_uniq:
				hapdat_uniq[hapdat[i][1]] = hapdat[i]
			else:
				#print(hapdat[i])
				if hapdat[i][0] == 'neand' and hapdat_uniq[hapdat[i][1]][0] != 'neand':
					hapdat_uniq[hapdat[i][1]][0] = hapdat[i][0]
		#print(len(hapdat))
		#print(len(hapdat_uniq))
		#print(hapdat_uniq)

		mod_ref_obs = 0.0
		mod_alt_obs = 0.0
		neand_ref_obs = 0.0
		neand_alt_obs = 0.0
		tot_mod = 0.0
		tot_neand = 0.0

		for key in hapdat_uniq:
			for line in vcf_dat:
				if '#'in line:
					continue
				else:
					#line = line.strip().split()
					if hapdat_uniq[key][2] == line[0] and hapdat_uniq[key][4] == line[1]: # find matching site in vcf
						#print(hapdat_uniq[key])
						hapid = hapdat_uniq[key][1].split('.')
						sampleid = hapid[1]
						samplehap = hapid[2]

						ref = line[3]
						alt = line[4]
						gt = line[sampleheader.index(sampleid)]
						gt = gt.split("|")

						if samplehap == '1':
							gt = gt[0]
							if gt == '0':
								allele = ref
							else: 
								allele = alt
						else:
							gt = gt[1]
							if gt == '0':
								allele = ref
							else: 
								allele = alt
						if hapdat_uniq[key][0] == 'null' and allele == ref:
							mod_ref_obs += 1
							tot_mod += 1
						elif hapdat_uniq[key][0] == 'null' and allele == alt:
							mod_alt_obs += 1
							tot_mod += 1
						elif hapdat_uniq[key][0] == 'neand' and allele == ref:
							neand_ref_obs += 1
							tot_neand += 1
						elif hapdat_uniq[key][0] == 'neand' and allele == alt:
							neand_alt_obs += 1
							tot_neand += 1
		if tot_mod > 0.0 and tot_neand > 0.0:
			outfile.write('%s\t%s\t%s\t%f\t%f\n' % (str(coord[0]), str(coord[1]), str(coord[2]), (mod_alt_obs/tot_mod), (neand_alt_obs/tot_neand)))
			#print(str(coord[0]), str(coord[1]), str(coord[2]), (mod_alt_obs/tot_mod), (neand_alt_obs/tot_neand))			
		else:
			continue
outfile.close()
