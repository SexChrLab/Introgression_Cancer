# Angela Oill
# 2021-01-06
# Fishers exact test

library("dplyr")
library(ggplot2)
library(cowplot)
library(ggpubr)


pop = "EUR"
# Read in data
dat.gene.counts <- read.table(
  paste("../../results/4.2_TCGA_somatic_genes/", pop, "_gene_sample_somatic_counts.txt", sep = ""),
  header = T)

all.results <- c()
for (i in 1:nrow(dat.gene.counts)) {
  #print(dat.gene.counts[i,])
  dattable <- matrix(c(as.numeric(as.character(dat.gene.counts[i,5])),
                       as.numeric(as.character(dat.gene.counts[i,3])),
                       as.numeric(as.character(dat.gene.counts[i,4])),
                       as.numeric(as.character(dat.gene.counts[i,2]))), ncol=2)
  colnames(dattable) = c('neand', 'non')
  rownames(dattable) = c('hasSomatic', 'noSomatic')
  #print(dattable)
  testresult <- fisher.test(dattable)
  #print(testresult$p.value)
  # add results to table
  rline <- cbind(dat.gene.counts[i,], testresult$p.value)
  all.results <- rbind(all.results, rline)
}

colnames(all.results) <- c("gene", "samples_non_intro_no_somatic", "samples_intro_no_somatic",
                           "samples_non_intro_somatic", "samples_intro_somatic",
                           "pval")
# correct for multiple testing
all.results.df <- as.data.frame(all.results)

adj.pVal <- p.adjust(as.numeric(as.character(all.results.df$pval)), method = "bonferroni")
all.results.df$pvalAdj <- adj.pVal

all.results.df$propInSom <- all.results.df$samples_intro_somatic / (all.results.df$samples_intro_somatic +
                                                                      all.results.df$samples_intro_no_somatic)
all.results.df$propNonInSom <- all.results.df$samples_non_intro_somatic / (all.results.df$samples_non_intro_somatic +
                                                                             all.results.df$samples_non_intro_no_somatic)

all.results.df.sig <- all.results.df[ which(all.results.df$pvalAdj <= 0.05), ]
all.results.df.sig.pre.cor <- all.results.df[ which(all.results.df$pval <= 0.05), ]

# NEW: Get number of genes where there are more introgressed samples with somatic
# mutations. Then get the number of gene where there are less introgressed samples
# with somatic mutations
more.int.som <- all.results.df.sig.pre.cor[ which(all.results.df.sig.pre.cor$samples_intro_no_somatic < all.results.df.sig.pre.cor$samples_intro_somatic), ]
less.int.som <- all.results.df.sig.pre.cor[ which(all.results.df.sig.pre.cor$samples_intro_no_somatic > all.results.df.sig.pre.cor$samples_intro_somatic), ]
equal.int.som <- all.results.df.sig.pre.cor[ which(all.results.df.sig.pre.cor$samples_intro_no_somatic == all.results.df.sig.pre.cor$samples_intro_somatic), ]


# output results
write.table(all.results.df, paste("../../results/4.2_TCGA_somatic_genes/", pop, "_gene_sample_somatic_fishers_results.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
write.table(all.results.df.sig, paste("../../results/4.2_TCGA_somatic_genes/", pop, "_gene_sample_somatic_fishers_results.sig.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
write.table(all.results.df.sig.pre.cor, paste("../../results/4.2_TCGA_somatic_genes/", pop, "_gene_sample_somatic_fishers_results.sig.no.mult.testing.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
