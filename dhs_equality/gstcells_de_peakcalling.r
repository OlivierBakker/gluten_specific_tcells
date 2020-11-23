library(DiffBind)

dhs <- dba(sampleSheet="data/de_peakcalling_samplesheet.csv")

pdf(width=10, height=10, file="output/plots/dba_dhs_correlation_heatmap.pdf")
plot(dhs)
dev.off()


dhs <- dba.count(dhs)


rpkm <- sapply(dhs$peaks, function(x){return( x[,"RPKM"])})
x <- dhs$peaks[[1]]
rownames(rpkm) <- paste0(x[,"Chr"],":", x[,"Start"], "-", x[,"End"])
colnames(rpkm) <- dhs$samples$SampleID
write.table(rpkm, file="output/dhs_consensus_counts_rpkm.tsv", quote=F, sep="\t")

counts <- sapply(dhs$peaks, function(x){return( x[,"Reads"])})
x <- dhs$peaks[[1]]
rownames(counts) <- paste0(x[,"Chr"],":", x[,"Start"], "-", x[,"End"])
colnames(counts) <- dhs$samples$SampleID
write.table(counts, file="output/dhs_consensus_counts_raw_counts.tsv", quote=F, sep="\t")

pdf(width=10, height=10, file="output/plots/dba_dhs_counts_correlation_heatmap.pdf")
plot(dhs)
dev.off()


pdf(width=10, height=10, file="output/plots/dba_dhs_counts_pca.pdf")
dba.plotPCA(dhs, label=DBA_CONDITION)
dev.off()


tmp <- dba.contrast(dhs, dhs$masks$t0, dhs$masks$t10, "t0", "t10")
tmp <- dba.contrast(tmp, dhs$masks$t10, dhs$masks$t30, "t10", "t30")
tmp <- dba.contrast(tmp, dhs$masks$t30, dhs$masks$t180, "t30", "t180")
tmp <- dba.analyze(tmp)

dhs.t10  <- dba.report(tmp, 1)
dhs.t30  <- dba.report(tmp, 2)
dhs.t180 <- dba.report(tmp, 3)

pdf(width=6, height=6, file="output/plots/dba_dhs_volcano_all_timepoints.pdf")
dba.plotVolcano(tmp, 1)
dba.plotVolcano(tmp, 2)
dba.plotVolcano(tmp, 3)
dev.off()


write.table(dhs.t180, sep="\t", quote=F, file="dhs_t30v180_de_peaks.tsv")


# Consensus peaks
library(DiffBind)
dhs <- dba(sampleSheet="data/de_peakcalling_samplesheet.csv")
dhs_consensus <- dba.peakset(dhs, consensus=c(DBA_CONDITION), minOverlap=0.66)

write.table(dhs_consensus$peaks[[21]], file="output/consensus_peaks_t10.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(dhs_consensus$peaks[[22]], file="output/consensus_peaks_t180.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(dhs_consensus$peaks[[23]], file="output/consensus_peaks_t0.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(dhs_consensus$peaks[[24]], file="output/consensus_peaks_t30.bed", quote=F, row.names=F, col.names=F, sep="\t")

dhs_consensus <- dba.peakset(dhs, consensus=c(DBA_TREATMENT), minOverlap=0.66)
write.table(dhs_consensus$peaks[[21]], file="output/consensus_peaks_all.bed", quote=F, row.names=F, col.names=F)
