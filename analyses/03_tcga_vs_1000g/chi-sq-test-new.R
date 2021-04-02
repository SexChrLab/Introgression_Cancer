# Angela Oill
# 01-05-2020
# Chi square test to compare proportion of introgression between TCGA and
# 1000 genomes
library("dplyr")
library("ggplot2")
library(cowplot)
library(ggpubr)

pop = "EUR"
# Read in data
dat.1000g <- read.table(
  paste("../../results/3.1_TCGA_vs_1000_genomes/prop_Neand/All.chr_", pop, "_1000g_introgression_per_gene.txt", sep = ""),
  header = T)
dat.tcga <- read.table(
  paste("../../results/3.1_TCGA_vs_1000_genomes/prop_Neand/All.chr_", pop,
        "_tcga_tagSNPs_intersect_alt_frq_sites_introgression_per_gene.txt",
        sep = ""),
  header = T)
# remove extra headers
dat.1000g.filter <- dat.1000g[ which(dat.1000g$PropNean!="PropNean"),]
dat.tcga.filter <- dat.tcga[ which(dat.tcga$PropNean!="PropNean"),]

# Merge data
dat.merge <- merge(dat.tcga.filter, dat.1000g.filter, by = "geneName")

# only sites where higher in tcga
dat.merge.higher.TCGA <- dat.merge[ which(as.numeric(as.character(dat.merge$PropNean.x)) >
                                        as.numeric(as.character(dat.merge$PropNean.y))),]
# there are 780 genes where the proportion of Neanderthal introgression is higher in TCGA than 1000 genomes

#dat.merge <- merge(dat.tcga, dat.1000g, by = "geneName")
# there are some rows in the merged file that dont have data but just the
# contents of the header. remove those rows
dat.merge.filter <- dat.merge[ which(dat.merge$PropNean.y!="PropNean"),]


# for every gene that overlaps between 1000 genomes and tcga perform chi square
all.results <- c()
for (i in 1:nrow(dat.merge.filter)) {
  #print(dat.merge[i,1])
  # number nean and not neand in tcga
  non.tcga <- as.numeric(as.character(dat.merge.filter[i,5])) - as.numeric(as.character(dat.merge.filter[i,6]))
  non.1000g <- as.numeric(as.character(dat.merge.filter[i,11])) - as.numeric(as.character(dat.merge.filter[i,12]))
  n.tcga <- as.numeric(as.character(dat.merge.filter[i,6]))
  n.1000g <- as.numeric(as.character(dat.merge.filter[i,12]))
  dattable <- matrix(c(n.1000g, n.tcga, non.1000g, non.tcga), ncol=2)
  colnames(dattable) = c('neand', 'non')
  rownames(dattable) = c('1000g', 'tcga')
  #print(dattable)
  testresult <- chisq.test(dattable)
  #print(testresult$p.value)
  # add results to table
  rline <- cbind(as.character(dat.merge.filter[i,2]), as.character(dat.merge.filter[i,3]), as.character(dat.merge.filter[i,4]),
                 as.character(dat.merge.filter[i,1]),
                 n.1000g, n.tcga, non.1000g, non.tcga,
                 as.character(dat.merge.filter[i,13]), as.character(dat.merge.filter[i,7]),
                 testresult$p.value)
  all.results <- rbind(all.results, rline)
}

# correct for multiple testing
colnames(all.results) <- c("chr", "posStart", "posEnd", "gene",
                           "Nean1000g", "NeanTCGA", "NonNean1000g", "NonNeanTCGA",
                           "propNean1000g", "propNeanTCGA",
                           "pval")
#adj.pVal <- p.adjust(as.numeric(all.results[,2]), method = "bonferroni")
all.results.df <- as.data.frame(all.results)
adj.pVal <- p.adjust(as.numeric(as.character(all.results.df$pval)), method = "bonferroni")
all.results.df$pvalAdj <- adj.pVal

# sort results
df.sorted <- arrange(all.results.df, chr, posStart)

df.sorted.sig <- df.sorted[ which(df.sorted$pvalAdj <= 0.05), ]
df.sorted.sig2 <- df.sorted[ which(as.numeric(as.character(df.sorted$pval)) <= 0.05), ]

# grab genes that are higher in TCGA than 1000 genomes
dat.higher.TCGA <- df.sorted[ which(as.numeric(as.character(df.sorted$propNeanTCGA)) >
                                            as.numeric(as.character(df.sorted$propNean1000g))),]
dat.higher.TCGA <- df.sorted.sig2[ which(as.numeric(as.character(df.sorted.sig2$propNeanTCGA)) >
                                      as.numeric(as.character(df.sorted.sig2$propNean1000g))),]
dat.lower.TCGA <- df.sorted.sig2[ which(as.numeric(as.character(df.sorted.sig2$propNeanTCGA)) <
                                           as.numeric(as.character(df.sorted.sig2$propNean1000g))),]
dat.equal.TCGA <- df.sorted.sig2[ which(as.numeric(as.character(df.sorted.sig2$propNeanTCGA)) ==
                                          as.numeric(as.character(df.sorted.sig2$propNean1000g))),]


# output results
write.table(df.sorted, paste("../../results/3.1_TCGA_vs_1000_genomes/prop_Neand/", pop, "_chi_results.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
write.table(df.sorted.sig, paste("../../results/3.1_TCGA_vs_1000_genomes/prop_Neand/", pop, "_chi_results.sig.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
write.table(df.sorted.sig2, paste("../../results/3.1_TCGA_vs_1000_genomes/prop_Neand/", pop, "_chi_results.sig.before.correction.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)

write.table(dat.higher.TCGA, paste("../../results/3.1_TCGA_vs_1000_genomes/prop_Neand/", pop, "_propN_higher_in_TCGA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
write.table(dat.lower.TCGA, paste("../../results/3.1_TCGA_vs_1000_genomes/prop_Neand/", pop, "_propN_lower_in_TCGA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
write.table(dat.equal.TCGA, paste("../../results/3.1_TCGA_vs_1000_genomes/prop_Neand/", pop, "_propN_equal_in_TCGA.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
