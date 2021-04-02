#!/usr/bin/env python3
# Angela Oill
# calc.perc.and.avg.introgression.py
# This script will take output from bedtools intersect (-wo) and ncbi RefSeq gene
# file and calculate:
# 1) the % of introgressed sequence per haplotype, and
# 2) the average Neanderthal introgression across exons across all haplotypes
import csv
from optparse import  OptionParser
import sys
################################################################################
USAGE = """
python calc.perc.and.avg.introgression.py	--coords <path and name to NCBIRefSeq gene coordinates >
                                --bedtfn <path and name to bedtools intersect results file >
								--out1 <path and name of results output file 1 >
                                --out2 <path and name of results output file 2 >

coords == path and name to NCBIRefSeq gene coordinates
bedtfn == path and name to bedtools intersect results file
out == path and stem name of results output file
"""

parser = OptionParser(USAGE)
parser.add_option('--coords',dest='coords', help = 'path and name to NCBIRefSeq gene coordinates')
parser.add_option('--bedtfn',dest='bedtfn', help = 'path and name to bedtools intersect results file')
parser.add_option('--out1',dest='out1', help = 'path and stem name of results output file 1')
parser.add_option('--out2',dest='out2', help = 'path and stem name of results output file 2')

(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.coords is None:
	parser.error('path and name to NCBIRefSeq gene coordinates not given')
if options.bedtfn is None:
	parser.error('path and name to bedtools intersect results file not given')
if options.out1 is None:
	parser.error('path and stem name of results output file 1 not given. .introgression.per.hap.txt')
if options.out2 is None:
	parser.error('path and stem name of results output file 2 not given. .avg.introgression.txt')
################################################################################
# Read in ncbi file
with open(options.coords, "r") as bed_file:
    windows = list(csv.reader(bed_file, delimiter='\t'))
    print("Opened and processed NCBIRefSeq file")
# Read in bedtools results file
with open(options.bedtfn, "r") as hap_fn: # dont have to open and process lines every time. i dont think...
    hap_fn_list = [line.strip().split() for line in hap_fn]
    print("Opened and initially processed bedtools intersect results file")

# Set up output files
# Output file 1
#outfilename1 = options.out + ".introgression.per.hap.txt"
outfilename1 = options.out1
outfile1 = open(outfilename1, "w")
# Add a header to this output file
header1 = "hap_name\tchromosome\tstart\tend\tgene_name\tneand\tden\tmod\tambig\tseq_neanderthal\n"
outfile1.write(header1)

# Output file 2
#outfilename2 = options.out + ".avg.introgression.txt"
outfilename2 = options.out2
outfile2 = open(outfilename2, "w")
# Add a header to this output file
header2 = "chromosome\tstart\tend\tgene_name\thaps\tseq_neanderthal\n"
outfile2.write(header2)

print("\nStarting main analysis")
# Iterate through each gene in bed file
for w in windows:
    chr = w[0]
    txStart = w[1]
    txEnd = w[2]
    name2 = w[3] # gene name 
    
    #exonStarts = w[9].split(',')
    #exonEnds = w[10].split(',')
    #exonCount = int(w[8])
    #chr = w[2]
    #txStart = w[4]
    #txEnd = w[5]
    #name1 = w[1]
    #name2 = w[12]
    #print("Starting analysis for gene: ", name1, name2, ". There are a total of ", exonCount, " exons to process.")
    toprint = ["Starting analysis for gene: " + str(name2) + "."]
    print("".join(toprint))
    # Initiate counters for number of total base pairs called Neanderthal,
    # Denisovan, Modern, Ambiguous for an indivdual haplotype.
    tot_Neand_bp_across_haps = 0.0
    tot_bp_across_haps = 0.0
    # make dictionary for each haplotype
    haps = dict()
    # Get bedtools entries that match the exon we are processing
    #hap_fn_list_subset = [ln for ln in hap_fn_list if str("chr" + chr) == ln[0] and (int(exonStarts[i])-1) == int(ln[1]) and int(exonEnds[i]) == int(ln[2])]
    hap_fn_list_subset = [ln for ln in hap_fn_list if str(chr) == ln[0] and int(txStart) == int(ln[1]) and int(txEnd) == int(ln[2])]

    #print(len(hap_fn_list_subset))
    # There are duplicates in hap_fn_list_subset because there are overlapping
    # windows. I need to figure out a way to account for this so an indivdual
    # haplotype does not get counted more than once.
    # Make a new list or dictionary with only one entry per haplotype. If the
    # indivdual haplotype has >1 entry and there is any evidence of neanderthal
    # introgression, count the haplotype for that exon as neanderthal.
    unique_hap_fn_list_subset = dict()
    for i in range(0, len(hap_fn_list_subset)):
        #print(hap_fn_list_subset[i][6], hap_fn_list_subset[i][7])
        if hap_fn_list_subset[i][8] not in unique_hap_fn_list_subset:
            unique_hap_fn_list_subset[hap_fn_list_subset[i][8]] = hap_fn_list_subset[i] # add indivdual hap as key
        else:
            if hap_fn_list_subset[i][7] == "neand" and unique_hap_fn_list_subset[hap_fn_list_subset[i][8]][7] != "neand":
                #print(hap_fn_list_subset[i][13], hap_fn_list_subset[i][12])
                #print(hap_fn_list_subset[i][13], unique_hap_fn_list_subset[hap_fn_list_subset[i][13]])
                unique_hap_fn_list_subset[hap_fn_list_subset[i][8]][7] = hap_fn_list_subset[i][7] # or "neand"
                #print("New", hap_fn_list_subset[i][13], unique_hap_fn_list_subset[hap_fn_list_subset[i][13]])
            else:
                continue
    #print(len(unique_hap_fn_list_subset))
    #print(unique_hap_fn_list_subset)

    # For each unique haplotype get number N, D, M, ambig
    # Add above info to dictionary
    for hkey in unique_hap_fn_list_subset:
        #print(hkey)
        if unique_hap_fn_list_subset[hkey][8] not in haps: # if the haplotype not in the dictionary
            # dictionary structure
            # indivdual.hap.id : #Neanbp, #Denbp, #Modbp, #Ambbp
            # if the haplotype is not currently in the dictionary
            # then we are just adding 1 to start and the other values
            # will be 0 because they are not observed yet.
            #print(hap_fn_list_subset_no_dups[i][13], hap_fn_list_subset_no_dups[i][12])
            if unique_hap_fn_list_subset[hkey][7] == "neand":
                haps[unique_hap_fn_list_subset[hkey][8]] = [int(unique_hap_fn_list_subset[hkey][9]), 0, 0, 0]
            elif unique_hap_fn_list_subset[hkey][7] == "den":
                haps[unique_hap_fn_list_subset[hkey][8]] = [0, int(unique_hap_fn_list_subset[hkey][9]), 0, 0]
            elif unique_hap_fn_list_subset[hkey][7] == "null":
                haps[unique_hap_fn_list_subset[hkey][8]] = [0, 0, int(unique_hap_fn_list_subset[hkey][9]), 0]
            else:
                haps[unique_hap_fn_list_subset[hkey][8]] = [0, 0, 0, int(unique_hap_fn_list_subset[hkey][9])]
        else: # If hap is in dictionary, add values for number N, D, M, A to dictionary for the existing haplotype
            #print(hap_fn_list_subset_no_dups[hkey][13], hap_fn_list_subset_no_dups[hkey][12])
            print("ERROR: duplicates")
            if unique_hap_fn_list_subset[hkey][7] == "neand":
                haps[unique_hap_fn_list_subset[hkey][8]][0] += int(unique_hap_fn_list_subset[hkey][9])
            elif unique_hap_fn_list_subset[hkey][7] == "den":
                haps[unique_hap_fn_list_subset[hkey][8]][1] += int(unique_hap_fn_list_subset[hkey][9])
            elif unique_hap_fn_list_subset[hkey][7] == "null":
                haps[unique_hap_fn_list_subset[hkey][8]][2] += int(unique_hap_fn_list_subset[hkey][9])
            else:
                haps[unique_hap_fn_list_subset[hkey][8]][3] += int(unique_hap_fn_list_subset[hkey][9])
    #print(haps) # this is a dictionary with the counts of Neand, Den, Mod, and
    # ambig across the gene

    # use dictionary holding haplotype as key and numbers of bps as values
    # iterate throug dictionary after collecting values for all bps across
    # all exons in isoform. When iterating, can calculate the average Neand
    # introgression
    #toprint = ["Adding results to ", outfilename1]
    #print("".join(toprint))
    #print(len(haps))
    if len(haps) > 0:
        tot_haps = 0.0
        for key in haps:
            #print(key)
            tot_haps += 1
            tot_Neand_bp_hap = float(haps[key][0])
            tot_bp_hap = float(haps[key][0]) + float(haps[key][1]) + float(haps[key][2]) + float(haps[key][3])
            tot_Neand_bp_across_haps += tot_Neand_bp_hap
            tot_bp_across_haps += tot_bp_hap
            per_intro_seq = tot_Neand_bp_hap/tot_bp_hap
            #print(tot_Neand_bp_across_haps, tot_bp_across_haps)
            #print("% seq neanderthal for hap: ", key, " = ", tot_Neand_bp_hap/tot_bp_hap)
            #print(key, chr, txStart, txEnd, name1, name2, exonCount, int(haps[key][0]), int(haps[key][1]), int(haps[key][2]), int(haps[key][3]), per_intro_seq)
            outfile1.write('%s\t%s\t%s\t%s\t%s\t' % (key, chr, txStart, txEnd, name2))
            outfile1.write('%i\t%i\t%i\t%i\t%f\n' % (int(haps[key][0]), int(haps[key][1]), int(haps[key][2]), int(haps[key][3]), per_intro_seq))
        #toprint = ["Adding results to ", outfilename2, "\n"]
        #print("".join(toprint))
        avg_N_intro = (tot_Neand_bp_across_haps/tot_haps)/(tot_bp_across_haps/tot_haps)
        outfile2.write('%s\t%s\t%s\t%s\t' % (chr, txStart, txEnd, name2))
        outfile2.write('%i\t%f\n' % (tot_haps, avg_N_intro))
    else:
        outfile2.write('%s\t%s\t%s\t%s\t' % (chr, txStart, txEnd, name2))
        outfile2.write('%i\t%s\n' % (0, "NA"))
    #print("Average Neanderthal introgression for isoform: ", w[1], w[12], " = ", avg_N_intro)

outfile1.close()
outfile2.close()

sys.exit("\nDONE")
#sys.exit(0)
#os._exit(os.EX_OK)
