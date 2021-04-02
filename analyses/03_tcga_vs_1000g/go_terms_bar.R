library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

dat.lower <- read.table("pantherChart_lower_eur.txt",
                   sep = "\t")
dat.lower$cat <- "Depleted_TCGA"
dat.lower$prop <- dat.lower$V3 / sum(dat.lower$V3)

dat.higher <- read.table("pantherChart_higher_eur.txt",
                        sep = "\t")
dat.higher$cat <- "Enriched_TCGA"
dat.higher$prop <- dat.higher$V3 / sum(dat.higher$V3)

dat.tot <- rbind(dat.lower, dat.higher)

nb.cols <- length(unique(as.character(dat.tot$V2)))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

p1 <- ggplot(dat.tot, aes(fill=V2, y=prop, x=cat)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = mycolors) +
  labs(fill = "GO Term") +
  ggtitle("Europeans") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = .5, vjust = 0.5)) +
  ylab("Proportion") + xlab("") +
  scale_x_discrete(labels=c("Depleted in TCGA", "Enriched in TCGA"))


dat.lower <- read.table("pantherChart_lower_asn.txt",
                        sep = "\t")
dat.lower$cat <- "Depleted_TCGA"
dat.lower$prop <- dat.lower$V3 / sum(dat.lower$V3)

dat.higher <- read.table("pantherChart_higher_asn.txt",
                         sep = "\t")
dat.higher$cat <- "Enriched_TCGA"
dat.higher$prop <- dat.higher$V3 / sum(dat.higher$V3)

dat.tot <- rbind(dat.lower, dat.higher)

nb.cols <- length(unique(as.character(dat.tot$V2)))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

p2 <- ggplot(dat.tot, aes(fill=V2, y=prop, x=cat)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = mycolors) +
  labs(fill = "GO Term") +
  ggtitle("East Asians") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = .5, vjust = 0.5)) +
  ylab("Proportion") + xlab("") +
  scale_x_discrete(labels=c("Depleted in TCGA", "Enriched in TCGA"))

pdf("FigS2.pdf", width = 12, height = 12)

ggarrange(
  p1, p2,
  nrow = 2,
  labels = c("A", "B"),
  #common.legend = T,
  #legend="right",
  align = "v"
)
dev.off()
