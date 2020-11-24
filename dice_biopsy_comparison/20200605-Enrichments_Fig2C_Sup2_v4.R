library(plyr)
library(dplyr)
library(tidyr)
library(DOSE)
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(biomaRt)
library(stringi)
library(ComplexHeatmap)
library(gtools)
library(circlize)
library(scico)

#SymbToEntrez: Uses ENSEMBL and SYMBOL to get ENTREZID, wich improves the mapping of genes.
#Input: Dataframe with at least two columns named as "ENSEMBL" and "SYMBOL"
#Output: Original Dataframe with three more columns: "ENTREZID.x" with entrezId found using ENSEMBL, "ENTREZID.y" with entrezID found using SYMBOL, and "ENTREZID" with a merging of ENTREZID.x and ENTREZID.y.
ENSEMBLToEntrez <- function(f) {
  df <- f
  gene <- df$ENSEMBL
  gene.df <- bitr(gene, fromType = "ENSEMBL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  df <- dplyr::full_join(df, gene.df, by = "ENSEMBL")
  return(df)
}

ENSEMBLToEntrezSYMBOL <- function(f) {
  df <- f
  gene <- df$ENSEMBL
  gene.df <- bitr(gene, fromType = "ENSEMBL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  df <- dplyr::full_join(df, gene.df, by = "ENSEMBL")
  return(df)
}

#setting working directory
setwd("/Users/R_projects/gsTcells")

#reading data from gluten-specific T cells
df.de.gstcells <- read.csv(file = "/Users/R_projects/gsTcells/data/de_analysis/output/gs_tcells_deseq2_batch_and_gender_corrected_leiden_fdr_005_replication_olso.tsv", sep = "\t", header = T, stringsAsFactors = T)

df.de.gstcells <- as.data.frame(df.de.gstcells)

names(df.de.gstcells)[4] <- "ENSEMBL"

ensembl.genes <- unique(df.de.gstcells$ENSEMBL)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  values=ensembl.genes,
  mart=mart)
names(genes) <- c("ENSEMBL", "ENTREZID", "SYMBOL")
genes.original <- genes
genes <- genes[!duplicated(genes$ENSEMBL),]

df.de.gstcells <- dplyr::left_join(df.de.gstcells, genes, by = "ENSEMBL")


###### Cluster #####
#getting annotation data from 6 clusters
clusters10 <- read.csv(file = "/Users/R_projects/gsTcells/data/cluster_memberships_kmeans_nclust6.tsv", header = F, stringsAsFactors = F, sep = "\t")

#renaming colnames to match with environment
names(clusters10) <- c("ENSEMBL", "Cluster", "SYMBOL", "Function")

#making sure Cluster column is a string vector to proper use mtach() and paste()
clusters10$Cluster <- as.character(clusters10$Cluster)

#reassigning cluster number based on discussion with Iris and Olivier
clusters10$Cluster <- paste(clusters10$Cluster, "A", sep = "")
clusters10$Cluster <- gsub("4A","1", clusters10$Cluster)
clusters10$Cluster <- gsub("6A","2", clusters10$Cluster)
clusters10$Cluster <- gsub("5A","3", clusters10$Cluster)
clusters10$Cluster <- gsub("2A","4", clusters10$Cluster)
clusters10$Cluster <- gsub("1A","5", clusters10$Cluster)
clusters10$Cluster <- gsub("3A","6", clusters10$Cluster)

#making list
list.clusters10 <- split(clusters10, clusters10$Cluster)

#function to give the correct format
formatting.gene.list <- function(gene.list, deg.df){
  list <- gene.list
  list <- left_join(list, deg.df, by = "ENSEMBL")
  list <- list[c("ENTREZID","log2FoldChange")]
  return(list)
}
#formatting list
list.clusters10.merged.deg <- lapply(list.clusters10, formatting.gene.list, deg.df = df.de.gstcells)

#Verication: 2 columns named "ENTREZID", "log2FoldChange"
list.clusters10.merged.deg[1]
#verification: length should be 3509, which is number of genes in cluster file
length(intersect(df.de.gstcells$ENSEMBL, clusters10$ENSEMBL))

#New format for oing enrichment, sorting genes based on L2FC
formatting.gene.list.second <- function(gene.list.2columns){
  df <- gene.list.2columns
  geneList = df$log2FoldChange
  names(geneList) = as.character(df$ENTREZID)
  geneList = sort(geneList, decreasing = TRUE)
  return(geneList)
}
#New format to only include gene names already sorted
gene.list.only.names <- function(gene.list.vector){
  df <- gene.list.vector
  return(names(df))
}

#list with vector names as genes and L2FC
list.clusters10.merged.deg.reactomeformatted <- lapply(list.clusters10.merged.deg , formatting.gene.list.second)
#list with only gene names
list.clusters10.merged.deg.reactomeformatted.2 <- lapply(list.clusters10.merged.deg.reactomeformatted, gene.list.only.names)

#REACTOME enrichment
res.df <- compareCluster(list.clusters10.merged.deg.reactomeformatted.2, fun="enrichPathway", readable=T)
#GO: Molecular function
mf <- compareCluster(list.clusters10.merged.deg.reactomeformatted.2, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "MF") 
#GO: Biological process
bp <- compareCluster(list.clusters10.merged.deg.reactomeformatted.2, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "BP") 
#GO: Celullar component
cc <- compareCluster(list.clusters10.merged.deg.reactomeformatted.2, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "CC") 
#KEGG enrichment
kegg <- compareCluster(list.clusters10.merged.deg.reactomeformatted.2, fun="enrichKEGG") 

# Convert adjusted p to -log10 scale
levels(res.df@compareClusterResult$Cluster) <- paste0("Cluster ", 1:6)
res.df@compareClusterResult$p.adjust <- -log10(res.df@compareClusterResult$p.adjust)
#write.csv(res.df, file = "./Results/reactome.gstcells.clusters.csv")

#PLOTS
dotplot(res.df, font.size = 10, showCategory = 15)
dotplot(mf, font.size = 10, showCategory = 10)
dotplot(bp, font.size = 10, showCategory = 10)
dotplot(cc, font.size = 10, showCategory = 10)
dotplot(kegg, font.size = 10, showCategory = 10)

enrichMap(bp)

#Plots using cnetplot function
x <- enrichPathway(gene= list.clusters10.merged.deg.reactomeformatted.2[["1"]], pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory = 15)
cnetplot(x, categorySize="pvalue", foldChange = list.clusters10.merged.deg.reactomeformatted[["1"]], colorEdge = TRUE, showCategory = 10)

x <- enrichPathway(gene= list.clusters10.merged.deg.reactomeformatted.2[["2"]], pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory = 15)
cnetplot(x, categorySize="pvalue", foldChange = list.clusters10.merged.deg.reactomeformatted[["2"]], colorEdge = TRUE, showCategory = 10)

x <- enrichPathway(gene= list.clusters10.merged.deg.reactomeformatted.2[["3"]], pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory = 15)
cnetplot(x, categorySize="pvalue", foldChange = list.clusters10.merged.deg.reactomeformatted[["3"]], colorEdge = TRUE, showCategory = 10)

x <- enrichPathway(gene= list.clusters10.merged.deg.reactomeformatted.2[["4"]], pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory = 15)
cnetplot(x, categorySize="pvalue", foldChange = list.clusters10.merged.deg.reactomeformatted[["4"]], colorEdge = TRUE, showCategory = 10)

x <- enrichPathway(gene= list.clusters10.merged.deg.reactomeformatted.2[["5"]], pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory = 15)
cnetplot(x, categorySize="pvalue", foldChange = list.clusters10.merged.deg.reactomeformatted[["5"]], colorEdge = TRUE, showCategory = 10)

x <- enrichPathway(gene= list.clusters10.merged.deg.reactomeformatted.2[["6"]], pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory = 15)
cnetplot(x, categorySize="pvalue", foldChange = list.clusters10.merged.deg.reactomeformatted[["6"]], colorEdge = TRUE, showCategory = 10)