#!/usr/bin/Rscript
# Author: Angela Taravella Oill
# ___usage___
# Rscript alt.allele.frqs.mod.intro.haps.1000g.R <modern_and_introgressed_alt_frq_results.txt> <keep_sites_out_file.txt> <output_plot_name.pdf>
# Description:
# This script will take as input the alternative allele frequency files from the 
# 1000 genomes introgressed haplotype data and output sites where the alternative 
# allele is >=.9 on modern haplotype and <=.1 on the neanderthal haplotype. and 
# sites where the alternative allele is <=.1 on modern haplotype and =>.9 on the
# neanderthal haplotype. This script will also output a heatmap type plot showing 
# the alt allele frequency on both haplotypes pre and post filtering.

### LOAD LIBRARIES ###
library(ggplot2)
library(tidyr)
library(ggpubr)

### SET UP COMMAND INPUTS ###
args = commandArgs(trailingOnly=TRUE)
frqfn = args[1] # alt frq results text file
out1 = args[2] # keep list output
out2 = args[3] # figure output


#dat.in <- read.table("All.chr_EUR_1000g_tcga_tag_sites_mod_intro_alt_frq.txt", header = T)
dat.in <- read.table(frqfn, header = T)
dat.in.new <- dat.in %>% drop_na()
dat.in.new$mod <-as.numeric(as.character(dat.in.new[,4]))
dat.in.new$intro <-as.numeric(as.character(dat.in.new[,5]))
dat.in.sub <- subset(dat.in.new, mod>=0.9 & intro<=0.1 | 
                       intro>=0.9 & mod<=0.1)

dat.keep.list.vars <- names(dat.in.sub) %in% c("chr", "end_pos")
dat.keep.list <- dat.in.sub[dat.keep.list.vars]
colnames(dat.keep.list) <- c("#chr", "pos")

# output file 1
#write.table(dat.keep.list, "test.txt", sep = "\t", quote = F, row.names = F)
write.table(dat.keep.list, out1, sep = "\t", quote = F, row.names = F)

# plot and output file 2
pdf(out2, width = 8, height = 5, onefile=F)
#pdf("out.pdf", width = 8, height = 5, onefile=F)
a <- ggplot(dat.in.new, aes(mod, intro)) +
  geom_bin2d(bins = 40, color ="white")+
  scale_fill_gradient(low =  "blue", high = "red")+
  theme_minimal() +
  xlab("Alt. allele frq. on introgressed haplotype") + ylab("Alt. allele frq. on non-introgressed haplotype")

b <- ggplot(dat.in.sub, aes(mod, intro)) +
  geom_bin2d(bins = 40, color ="white")+
  scale_fill_gradient(low =  "blue", high = "red")+
  theme_minimal() +
  xlab("Alt. allele frq. on introgressed haplotype") + ylab("Alt. allele frq. on non-introgressed haplotype")

ggarrange(a, b, ncol = 2, nrow = 1, align="h",
          common.legend = T, labels = c("A", "B"), legend = "right")
dev.off()
