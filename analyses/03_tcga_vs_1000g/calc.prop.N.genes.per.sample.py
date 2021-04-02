#!/usr/bin/env python3
import csv
from optparse import  OptionParser
"""
# Background:
# For this script, I want to take as input 1) TCGA vcf with tag snp sites that
# pass filtering, 2) bed file with gene coordinates, and 3) tag snp bed file
# with information on what the Neanderthal allele is at each site

# With this input, I want to output the following:
# 1) For each sample and each gene:
chr posStart posEnd geneName sampleID numTagSites numNeansites numModsites call
# 2) For each gene across all samples, get the proportion of Neanderthal
# introgression across all samples
chr posStart posEnd geneName totNumSamples NumNean NumMod PropNean
"""
################################################################################
USAGE = """
python3 calc.prop.N.genes.per.sample.py	--vcf <path and name to vcf >
                                --coord <path and name to NCBI gene coordinates file>
								--tagbed <path and name to tag snp bed file >
								--out <path and stem name for output files>

vcf == path and name to TCGA filtered vcf
coord == path and name to NCBI gene coordinates file
tagbed == path and name to tag snp bed file
out == path and stem name for output file
"""

parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcf', help = 'path and name to TCGA vcf')
parser.add_option('--coord',dest='coord', help = 'path and name to NCBI gene coordinates file')
parser.add_option('--tagbed',dest='tagbed', help = 'path and name to tag snp bed file')
parser.add_option('--out',dest='out', help = 'path and stem name to output file')

(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.vcf is None:
	parser.error('path and name to TCGA vcf not given')
if options.coord is None:
	parser.error('path and name to NCBI gene coordinates file not give')
if options.tagbed is None:
	parser.error('path and name to tag snp bed file not give')
if options.out is None:
	parser.error('path and stem name for output file not give')
################################################################################
# Set up output files
# 1. Introgression call per sample per gene
outfn1 = (options.out + "_introgression_per_sample_per_gene.txt")
outfile1 = open(outfn1, 'w')
header1 = "chr\tposStart\tposEnd\tgeneName\ttotTagSites\tSampleID\ttagSiteswithData\tNeanTagSites\tpropTagSiteNean\tCall\n"
#header1 = "chr\tposStart\tposEnd\tgeneName\tSampleID\ttotTagSites\tNeanTagSites\tModTagSite\tCall\n"
outfile1.write(header1)

# 2. Proportion Neanderthal ancestry per gene across samples
outfn2 = (options.out + "_introgression_per_gene.txt")
outfile2 = open(outfn2, 'w')
header2 = "chr\tposStart\tposEnd\tgeneName\ttotSamples\tNeanSamples\tPropNean\n"
#header2 = "chr\tposStart\tposEnd\tgeneName\ttotSamples\tNeanSamples\tModSamples\tPropNean\n"
outfile2.write(header2)

with open(options.tagbed, "r") as tagBed:
    tagBedData = [line.strip().split() for line in tagBed]

with open(options.coord, "r") as geneCords:
    for gene in geneCords:
        #print("gene", gene)
        # get chr, gene start, gene end
        info = gene.strip().split()
        #print(info)
        chrm = info[0]
        chrmNum = chrm[3:] # just want the number
        geneStart = info[1]
        geneEnd = info[2]
        geneID = info[3]
        #print(chrmNum, geneStart, geneEnd, geneID)
        #gene = gene.strip().split()
        #gene = gene.strip()
        #print(gene)
        # get all entries (tag snps) in tagBedData that are in gene
        #tagBedDataSubset = [ln for ln in tagBedData if ln[15] == gene]
        tagBedDataSubset = [ln for ln in tagBedData if int(ln[0]) == int(chrmNum) and int(ln[2]) >= int(geneStart) and int(ln[2]) <= int(geneEnd)]
        #print(len(tagBedDataSubset))
        #print(tagBedDataSubset)
        totTag = len(tagBedDataSubset)
        # make sample dictionary
        geneSampleInfo = dict()
        for tag in tagBedDataSubset:
            #print("tag", tag[15])
            # go through every sample in vcf for that specific tagsnp in genelotype
            with open(options.vcf, "r") as vcf:
                for line in vcf:
                    if "##" in line:
                        continue
                    elif "#" in line:
                        sampleheader = line.strip().split()
                    else:
                        line = line.strip().split() # process VCF line
                        #print("line", line[1])
                        # is this line the tag snp the program is currently processing?
                        if line[0] == tag[0] and line[1] == tag[2]:
                            #print("match")
                            #print("tag pos", tag[0], tag[2], "ancestral allele", tag[3], "derived allele", tag[4], "neand base", tag[13], "den base", tag[14])
                            neandBase = tag[13]
                            #denBase = tag[14]
                            #ancestralBase = tag[3]
                            #derivedBase = tag[4]
                            # get info for each sample in vcf
                            #print(line)
                            ref = line[3]
                            alt = line[4]
                            #samplenumber = 0
                            samplenumber = 8
                            for n in line[9:]:
                                samplenumber += 1
                                gt = n.split(":")
                                if gt[0] == "./.": # no data for tag snp for this sample
                                    continue
                                else:
                                    # this sample has data so start adding this sample to geneSampleInfo dictionary
                                    if "|" in gt[0]:
                                        gt = gt[0].split("|")
                                    else:
                                        gt = gt[0].split("/")
                                    #print(gt)
                                    #print(n)
                                    gt1 = gt[0]
                                    gt2 = gt[1]
                                    #print(gt1, gt2)
                                    if str(gt1) == "0":
                                        gt1 = ref
                                    else:
                                        gt1 = alt
                                    if str(gt2) == "0":
                                        gt2 = ref
                                    else:
                                        gt2 = alt
                                    #print(gt1, gt2)
                                    sampleName = sampleheader[samplenumber]
                                    #print(sampleName)

                                    # We probably want neand, den, modern, and ambig
                                    # but for now we will do a crude calculation
                                    if sampleName not in geneSampleInfo:
                                        if gt1 == neandBase or gt2 == neandBase:
                                            #print(gt1, gt2, "neandBase: ", neandBase)
                                            geneSampleInfo[sampleName] = [1, 1] # first entry is number of tag snps with data, second entry is number of bases where at least 1 allele is Neanderthal
                                        else: # niether base is neanderthal
                                            geneSampleInfo[sampleName] = [1, 0]
                                    else: # the sample is already in dictionary so now we just keep adding info
                                        # add to tagWithData, add to Ntags
                                        if gt1 == neandBase or gt2 == neandBase:
                                            #print(gt1, gt2, "neandBase: ", neandBase)
                                            geneSampleInfo[sampleName][0] += 1
                                            geneSampleInfo[sampleName][1] += 1
                                        else: # niether base is neanderthal
                                            geneSampleInfo[sampleName][0] += 1 # only add to number of tag snps with data
                                    """
                                    # In this part of the code I will make the dictionary have the following entries
                                    # Key = sample ID
                                    # Value 1 = number of tag snps with data
                                    # Value 2 = number of bases where at least 1 allele is Neanderthal and it is does not match the ancestral allele
                                    # Value 3 = number of bases where at least 1 allele is Denisovan and it is does not match the ancestral allele
                                    # Value 4 = If there is a Neanderthal allele or Den allele that matches the ancestral allele, this will be called ambig
                                    # Value 5 = number of bases that niether matches Neanderthal or Denisovan allele
                                    if sampleName not in geneSampleInfo:
                                        if (gt1 == neandBase or gt2 == neandBase) and (gt1 != ancestralBase or gt2 != ancestralBase):
                                            #print("Call = Neanderthal", gt1, gt2, "neandBase: ", neandBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName] = [1, 1, 0, 0, 0]
                                        elif (gt1 == denBase or gt2 == denBase) and (gt1 != ancestralBase or gt2 != ancestralBase):
                                            #print("Call = Denisovan", gt1, gt2, "denBase: ", denBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName] = [1, 0, 1, 0, 0]
                                        elif (gt1 == neandBase or gt2 == neandBase) and (gt1 == ancestralBase or gt2 == ancestralBase):
                                            #print("Call = Ambig", gt1, gt2, "neandBase: ", neandBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName] = [1, 0, 0, 1, 0]
                                        elif (gt1 == denBase or gt2 == denBase) and (gt1 == ancestralBase or gt2 == ancestralBase):
                                            #print("Call = Ambig", gt1, gt2, "denBase: ", denBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName] = [1, 0, 0, 1, 0]
                                        else: # niether base is neanderthal or ambiguous so is mod
                                            #print("Call = Mod", gt1, gt2, "neandBase: ", neandBase, "denBase: ", denBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName] = [1, 0, 0, 0, 1]
                                    else: # the sample is already in dictionary so now we just keep adding info
                                        if (gt1 == neandBase or gt2 == neandBase) and (gt1 != ancestralBase or gt2 != ancestralBase):
                                            #print("Call = Neanderthal", gt1, gt2, "neandBase: ", neandBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName][0] += 1
                                            geneSampleInfo[sampleName][1] += 1
                                        elif (gt1 == denBase or gt2 == denBase) and (gt1 != ancestralBase or gt2 != ancestralBase):
                                            #print("Call = Denisovan", gt1, gt2, "denBase: ", denBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName][0] += 1
                                            geneSampleInfo[sampleName][2] += 1
                                        elif (gt1 == neandBase or gt2 == neandBase) and (gt1 == ancestralBase or gt2 == ancestralBase):
                                            #print("Call = Ambig", gt1, gt2, "neandBase: ", neandBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName][0] += 1
                                            geneSampleInfo[sampleName][3] += 1
                                        elif (gt1 == denBase or gt2 == denBase) and (gt1 == ancestralBase or gt2 == ancestralBase):
                                            #print("Call = Ambig", gt1, gt2, "denBase: ", denBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName][0] += 1
                                            geneSampleInfo[sampleName][3] += 1
                                        else: # niether base is neanderthal or ambiguous so is mod
                                            #print("Call = Mod", gt1, gt2, "neandBase: ", neandBase, "denBase: ", denBase, "ancestral allele: ", ancestralBase)
                                            geneSampleInfo[sampleName][0] += 1
                                            geneSampleInfo[sampleName][4] += 1"""
                        else:
                            continue
        """
        #print(chrmNum, geneStart, geneEnd, geneID, totTag)
        #print(geneSampleInfo)
        if len(geneSampleInfo) > 0:
            neandTot = 0
            for key in geneSampleInfo:
                totTagPerSample = float(geneSampleInfo[key][0])
                totNeandTagPerSample = float(geneSampleInfo[key][1])
                totDenTagPerSample = float(geneSampleInfo[key][2])
                totAmbTagPerSample = float(geneSampleInfo[key][3])
                propTagNeand = totNeandTagPerSample/totTagPerSample
                if propTagNeand > 0: # we can change this to be more restrictive.
                #if totNeandTagPerSample > 0 and totDenTagPerSample == 0 and totAmbTagPerSample == 0:
                    call = "neand"
                    neandTot += 1
                    #print(geneSampleInfo[key])
                    #print(key, call)
                elif totNeandTagPerSample == 0 and totDenTagPerSample > 0 and totAmbTagPerSample == 0:
                    call = "den"
                    #print(geneSampleInfo[key])
                    #print(key, call)
                elif totNeandTagPerSample > 0 and totDenTagPerSample > 0 or totAmbTagPerSample > 0:
                    call = "ambig"
                    #print(geneSampleInfo[key])
                    #print(key, call)
                else:
                    call = "mod"
                    #print(geneSampleInfo[key])
                    #print(key, call)
                #outfile.write('%s\t%s\t%s\t%s\t%i\t%s\t%i\t%i\t%f\t%s\n' % (str(chrmNum), geneStart, geneEnd, geneID, int(totTag), key, int(geneSampleInfo[key][0]), int(geneSampleInfo[key][1]), float(propTagNeand), call))
                outfile.write('%s\t%s\t%s\t%s\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%f\t%s\n' % (str(chrmNum), geneStart, geneEnd, geneID, int(totTag), key, int(geneSampleInfo[key][0]), int(geneSampleInfo[key][1]), int(geneSampleInfo[key][2]), int(geneSampleInfo[key][3]), int(geneSampleInfo[key][4]), float(propTagNeand), call))
                #print(chrmNum, geneStart, geneEnd, geneID, totTag, key, geneSampleInfo[key][0], geneSampleInfo[key][1], propTagNeand, call)
                #print(chrmNum, geneStart, geneEnd, geneID, totTag, key, geneSampleInfo[key][0], geneSampleInfo[key][1], geneSampleInfo[key][2], geneSampleInfo[key][3], geneSampleInfo[key][4], propTagNeand, call)
            outfile2.write('%s\t%s\t%s\t%s\t%s\t%s\t%f\n' % (chrmNum, geneStart, geneEnd, geneID, str(len(geneSampleInfo)), neandTot, float(neandTot)/float(len(geneSampleInfo))))
            print(chrmNum, geneStart, geneEnd, geneID, len(geneSampleInfo), neandTot, float(neandTot)/float(len(geneSampleInfo)))
        else:
            continue
        """
        if len(geneSampleInfo) > 0:
            neandTot = 0
            for key in geneSampleInfo:
                totTagPerSample = float(geneSampleInfo[key][0])
                totNeandTagPerSample = float(geneSampleInfo[key][1])
                propTagNeand = totNeandTagPerSample/totTagPerSample
                if propTagNeand > 0: # we can change this to be more restrictive.
                    call = "neand"
                    neandTot += 1
                else:
                    call = "other"
                outfile1.write('%s\t%s\t%s\t%s\t%i\t%s\t%i\t%i\t%f\t%s\n' % (str(chrmNum), geneStart, geneEnd, geneID, int(totTag), key, int(geneSampleInfo[key][0]), int(geneSampleInfo[key][1]), float(propTagNeand), call))
                #print(chrmNum, geneStart, geneEnd, geneID, totTag, key, geneSampleInfo[key][0], geneSampleInfo[key][1], propTagNeand, call)
            outfile2.write('%s\t%s\t%s\t%s\t%s\t%s\t%f\n' % (chrmNum, geneStart, geneEnd, geneID, str(len(geneSampleInfo)), neandTot, float(neandTot)/float(len(geneSampleInfo))))
            #print(chrmNum, geneStart, geneEnd, geneID, len(geneSampleInfo), neandTot, float(neandTot)/float(len(geneSampleInfo)))
        else:
            continue


# Close output file when done
outfile1.close()
outfile2.close()
