#Olink analysis
dat <- read.table('~/Documents/UMCG/projects/gsTCC_paper/FINAL_FIGURES_SCRIPTS/Olink/Olink_R-input.txt', header=T)
dat <- read.table('~/Documents/UMCG/projects/gsTCC_paper/FINAL_FIGURES_SCRIPTS/Olink/Olink_R-input.txt', header=T)
dim(dat)
rownames(dat) <- dat$NPX
dat<- dat[,-1]
dat <- t(dat)
# Some protein names don't match gene names, add ENSG id
load("/Users/irishelenejonkers/Documents/UMCG/projects/gsTCC_paper/FINAL_FIGURES_SCRIPTS/datasets/ensembl_geneName_mapping.Rdata")
ensembl_geneName_mapping$gene_id <- rownames(ensembl_geneName_mapping)
nrow(dat[row.names(dat)%in%ensembl_geneName_mapping$geneName,])
dat2 <- merge(dat, ensembl_geneName_mapping, by.x=0, by.y="geneName", all.x=T, all.y=F) # many duplicates, remove alt chromosome
chromosomes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
dat2 <- dat2[as.character(dat2$chr) %in% chromosomes,] #92

# which ones are also DE?
DE <- read.table("/Users/irishelenejonkers/Documents/UMCG/projects/gsTCC_paper/2020/cluster_memberships_kmeans_nclust6_final.txt", header=T, sep='\t')
link <- dat2[dat2$Row.names%in%DE$gene_name,]

# means
link$mean_prot_blank <- rowMeans(link[,2:3])
link$mean_prot_0hrs <- rowMeans(link[,4:6])
link$mean_prot_4hrs <- rowMeans(link[, 7:9])
# stdevs
link$stdev_prot_blank <- apply((link[,2:3]),1,sd) 
link$stdev_prot_0hrs <- apply((link[,4:6]),1,sd) 
link$stdev_prot_4hrs <- apply((link[,7:9]),1,sd) 
# blank correction
link$blank_correct_0hrs <- link$mean_prot_0hrs - link$mean_prot_blank
link$blank_correct_4hrs <- link$mean_prot_4hrs - link$mean_prot_blank
# FC
link$fc_prot <- (link$blank_correct_4hrs+0.1)/(link$blank_correct_0hrs+0.1) 
nrow(link[link$fc_prot>1,]) 
nrow(link[link$fc_prot<1,]) 

#remove genes with no protein associated
link <- link[!(link$blank_correct_4hrs==0),] #30

# FC t.test 
for (i in 1:nrow(link)) {
  z <- t.test(x = link[i,4:6], y = link[i,7:9], alternative = c("two.sided"), paired = F)
  link$p.val_fc_prot[i] <-  z$p.value
} 

# add column with the max value of the means at 0 hrs and 4 hrs
for (i in 1:nrow(link)){
  link$max.means[i] <- max(link[i, 24:25])  
}

# which are significant?
sign.prot <- link[link$p.val_fc_prot<0.05,]
# 8 in total. 

# Make plots per significant cytokine:
#plot as dotplot
library(ggplot2)
library(ggpubr)
library(reshape2)

topCyt <- sign.prot$Row.names
# only get values to plot 
df <- cbind(sign.prot$`TCC-07-1_unstim_0`-sign.prot$mean_prot_blank, sign.prot$`TCC-11-1_unstim_0`-sign.prot$mean_prot_blank, sign.prot$`TCC-19-1_unstim_0`-sign.prot$mean_prot_blank,
            sign.prot$`TCC-07-1_aCD3/CD28_4h`-sign.prot$mean_prot_blank, sign.prot$`TCC-11-1_aCD3/CD28_4h`-sign.prot$mean_prot_blank, sign.prot$`TCC-19-1_aCD3/CD28_4h`-sign.prot$mean_prot_blank)
rownames(df) <- sign.prot$Row.names
colnames(df) <- c("unstim","unstim","unstim","4hrs_stim","4hrs_stim","4hrs_stim")

# Using melt to literaly melt a wide data.frame into a long data.frame
pData <- melt(df[topCyt,])
colnames(pData) <- c("protein", "stim","value")


# get summary stats function
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# get p-values
compare_means(value ~ stim, data = pData, method = "t.test",
              group.by = "protein")

# Basic dot plot

# New attempt, with cluster colors: we have this palette for clusters
# c(“#EF6461”, “#6B6174", “#FABC2A”, “#8576B6", “#3BB273”, “#2274A5")
facet_fill_color <- c("#FABC2A","#FABC2A","#FABC2A","#FABC2A","#FABC2A","#6B6174","#FABC2A","#FABC2A")

p<-ggplot(pData, aes(x=stim, y=value)) +
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.6,
               position=position_dodge(1),color = NA) +
  facet_grid(~protein, scale="free") +
  scale_y_continuous(breaks = (seq(-8, 16, 4)), labels = seq(-8, 16, 4), limits = c(-8,17)) +
  
  theme_bw()+
  theme(panel.grid.minor = element_blank()) +
  theme(strip.background = element_rect(fill=facet_fill_color))+
  theme(strip.text.x = element_text(size = 10, color = "black"))+
  ylab("relative protein level") +
  xlab("") +
  theme(axis.text.x = element_text(color = "black", size = 10, angle=45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),  
        axis.title.y = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  stat_summary(fun.data=data_summary, color="blue", position=position_dodge(1), stroke=0) +
  stat_compare_means(label = "p.signif", method = "t.test", label.y = 17, label.x = 1.5) 

p 



