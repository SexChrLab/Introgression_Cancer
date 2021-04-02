# Angela Oill 
# 2020-12-10
# Comparing the distribution of proportion of Neanderthal ancestry at the gene 
# level across all genes and cancer gene sets.

# This script takes as input: 1) input path to gene set files, 2) path to results 
# folder, 3) pop codes to analyze, 4) filter threshold
# And outputs 3 tables: Table 1: Medians, Table 2: Pvalues, Table 3: Pvalues, 
# corrected for multiple testing

###############################################################################
# Load Libraries #
library(dplyr)
###############################################################################

###############################################################################
# Input Variables #
# 1. Cancer gene sets to compare to all genes
cosmicFile <- "../../data/gene_sets/COSMIC.txt"
dnaRepFile <- "../../data/gene_sets/HumanDNARepairGenes.txt"
oncoFile <- "../../data/gene_sets/oncogenes.GeneSymbols.txt"
tsgFile <- "../../data/gene_sets/Human_TSGs.GeneSymbols.txt"
# 2. Results folder path and stem file name
resultPath <- "../../results/1.1_Gene_level/"
resultsStem <- "chrAll.genes.haplotype.bps.avg.introgression.sorted.fixed."
# 3. Populations 
#popCodes <- c("EAS", "EUR", "PNG", "SAS")
popCodes <- c("EAS", "EUR", "SAS")
# Keep results with a minimum of 50 haplotypes
minFilter <- 50
# 4. Output
outName1 <- paste(resultPath, "gene.level.medians.filter.",  minFilter, ".txt", 
                  sep = "")
outName2 <- paste(resultPath, "gene.level.p.val.filter.", minFilter, ".txt", 
                  sep = "")
outName3 <- paste(resultPath, "gene.level.p.val.corrected.filter.", minFilter, 
  ".txt", sep = "")

#outName1 <- paste(resultPath, "gene.level.medians.filter.",  minFilter, ".ONETAIL.txt", 
#                  sep = "")
#outName2 <- paste(resultPath, "gene.level.p.val.filter.", minFilter, ".ONETAIL.txt", 
#                  sep = "")
#outName3 <- paste(resultPath, "gene.level.p.val.corrected.filter.", minFilter, 
#                  ".ONETAIL.txt", sep = "")
###############################################################################

###############################################################################
# Start analysis #
#--------------------------
# Read in cancer gene sets
#--------------------------
cosmicDat <- read.table(cosmicFile, sep = "\t", header = T)
dnaRepDat <- read.table(dnaRepFile, sep = "\t", header = T)
oncoDat <- read.table(oncoFile, sep = "\t", header = T)
tsgDat <- read.table(tsgFile, sep = "\t", header = T)

#--------------------------------
# Set up data frames for results
#--------------------------------
medianResults <- c()
wilTestResults <- c()
multipleCorrectedResults <- c()


#---------------------------------------------------------------
# Get medians and perform statistical tests for each population
#---------------------------------------------------------------
for (pop in popCodes) {
  #print(pop)
  # read in data
  popDat <- read.table(paste(resultPath, resultsStem, pop, ".txt", 
                             sep = ""), sep = "\t", header = T)
  
  # filter data
  popDat.filtered <- popDat %>% filter(haps > minFilter)
  
  # filter cancer gene sets
  cosmicDatFiltered <- merge(popDat.filtered, cosmicDat, by.x = "gene_name", 
                             by.y = "GeneSymbol")
  dnaRepDatFiltered <- merge(popDat.filtered, dnaRepDat, by.x = "gene_name", 
                             by.y = "GeneSymbol")
  oncoDatFiltered <- merge(popDat.filtered, oncoDat, by.x = "gene_name", 
                           by.y = "GeneSymbol")
  tsgDatFiltered <- merge(popDat.filtered, tsgDat, by.x = "gene_name", 
                          by.y = "GeneSymbol")
  
  # Get medians for all genes and cancer gene sets
  allMedian <- median(popDat.filtered$seq_neanderthal)
  cosmicMedian <- median(cosmicDatFiltered$seq_neanderthal)
  dnaRepMedian <- median(dnaRepDatFiltered$seq_neanderthal)
  oncoMedian <- median(oncoDatFiltered$seq_neanderthal)
  tsgMedian <- median(tsgDatFiltered$seq_neanderthal)
  mediansMerged <- cbind(pop, allMedian, cosmicMedian, dnaRepMedian, 
                         oncoMedian, tsgMedian)
  medianResults <- rbind(medianResults, mediansMerged)
  
  # Perform Wilcox rank sum test to test difference between all gene and cancer gene sets
  # all-cosmic
  if (nrow(cosmicDatFiltered) == 0) {
    cosmicWiltest <- ""
    cosmicWiltest$p.value <- "NA"
  } else if (nrow(cosmicDatFiltered) > 0) {
    cosmicWiltest <- 
      wilcox.test(popDat.filtered$seq_neanderthal, cosmicDatFiltered$seq_neanderthal)
    #cosmicWiltest <- 
    #  wilcox.test(cosmicDatFiltered$seq_neanderthal, popDat.filtered$seq_neanderthal, 
    #              alternative = "greater")
  }
  # all-dna repair
  if (nrow(dnaRepDatFiltered) == 0) {
    dnaRepWiltest <- ""
    dnaRepWiltest$p.value <- "NA"
  } else if (nrow(dnaRepDatFiltered) > 0) {
    dnaRepWiltest <- 
      wilcox.test(popDat.filtered$seq_neanderthal, dnaRepDatFiltered$seq_neanderthal)
    #dnaRepWiltest <- 
    #  wilcox.test(dnaRepDatFiltered$seq_neanderthal, popDat.filtered$seq_neanderthal, 
    #              alternative = "greater")
  }
  # all-oncogenes
  if (nrow(oncoDatFiltered) == 0) {
    oncoWiltest <- ""
    oncoWiltest$p.value <- "NA"
  } else if (nrow(oncoDatFiltered) > 0) {
    oncoWiltest <- 
      wilcox.test(popDat.filtered$seq_neanderthal, oncoDatFiltered$seq_neanderthal)
    #oncoWiltest <- 
    #  wilcox.test(oncoDatFiltered$seq_neanderthal, popDat.filtered$seq_neanderthal, 
    #              alternative = "greater")
  }
  # all-tsgs
  if (nrow(tsgDatFiltered) == 0) {
    tsgWiltest <- ""
    tsgWiltest$p.value <- "NA"
  } else if (nrow(tsgDatFiltered) > 0) {
    tsgWiltest <- 
      wilcox.test(popDat.filtered$seq_neanderthal, tsgDatFiltered$seq_neanderthal)
    #tsgWiltest <- 
    #  wilcox.test(tsgDatFiltered$seq_neanderthal, popDat.filtered$seq_neanderthal, 
    #              alternative = "greater")
  }
  wilTestpop <- cbind(pop, cosmicWiltest$p.value, dnaRepWiltest$p.value, oncoWiltest$p.value,
                      tsgWiltest$p.value)
  wilTestResults <- rbind(wilTestResults, wilTestpop)
  
  # Perform multiple test correction
  adjustedP <- p.adjust(wilTestpop[1,2:5], "bonferroni")
  multipleCorrectedPop <- 
    cbind(pop, adjustedP[1], adjustedP[2], adjustedP[3], adjustedP[4])
  multipleCorrectedResults <- 
    rbind(multipleCorrectedResults, multipleCorrectedPop)
}

#-----------------------------
# Add column names to results 
#-----------------------------
colnames(medianResults) <- c("Population", "All", "COSMIC", "DNArepair", 
                             "Oncogenes", "TSGs")
colnames(wilTestResults) <- c("Population", "All-COSMIC", "ALL-DNArepair", 
                              "ALL-Oncogenes", "All-TSGs")
colnames(multipleCorrectedResults) <- c("Population", "All-COSMIC", "ALL-DNArepair", 
                                        "ALL-Oncogenes", "All-TSGs")

#----------------
# Output results
#----------------
write.table(medianResults, outName1, sep = "\t", quote = F, row.names = F)
write.table(wilTestResults, outName2, sep = "\t", quote = F, row.names = F)
write.table(multipleCorrectedResults, outName3, sep = "\t", quote = F, row.names = F)

