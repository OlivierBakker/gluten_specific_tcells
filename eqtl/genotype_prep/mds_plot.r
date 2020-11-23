
library(BEDMatrix)
library(wordcloud)

kg.annot <- read.table("~/Documents/projects/pr_gs_tcells/eqtl/genotype_prep/igsr_samples_EUR.tsv", stringsAsFactors = F, sep="\t", header=T, row.names = 1)

#t    <- read.table("~/Documents/projects/pr_gs_tcells/eqtl/genotype_prep/mds_plot.mds", header=T, stringsAsFactors = F)
t    <- read.table("~/Documents/tmp/1000G_eivind_olivier_merged.mds", header=T, stringsAsFactors = F)

n.tcell <- 40
cols <- c(rep("GSTCELLS", n.tcell), kg.annot[t$IID[n.tcell:nrow(t)], 3])

cols[cols=="GSTCELLS"] <- "red"
cols[cols=="CEU"] <- "green"
cols[cols=="TSI"] <- "blue"
cols[cols=="FIN"] <- "orange"
cols[cols=="IBS"] <- "purple"
cols[cols=="GBR"] <- "pink"

plot(t$C1, t$C2, pch=20, col=cols, xlab="MDS component 1", ylab="MDS component 2")
legend("topright", c("TCELL", "CEU", "TSI", "FIN", "IBS", "GBR"), fill=c("red", "green", "blue", "orange", "purple", "pink"))
tt <- t[1:n.tcell,]
xlim <- c(min(t$C1), max(t$C1))
ylim <- c(min(t$C2), max(t$C2))

textplot(tt$C1, tt$C2, tt$IID, xlim=xlim, ylim=ylim, new=F, cex=0.75)


## Validate with PCA
geno <- BEDMatrix("~/Documents/tmp/1000G_eivind_olivier_merged.bed")

geno <- na.omit(t(as.matrix(geno)))
pc <- prcomp(t(geno))
plot(pc$x[,1], pc$x[,2], pch=20, col=c(rep("red", 40), rep("blue", 503)))
textplot(pc$x[1:40,1], pc$x[1:40,2], rownames(pc$x)[1:40], new=F)


