## enrichment of non-coding RNAs in cluster 1 - 6 of the DE genes in gs T cells

setwd("~/Documents/UMCG/projects/gsTCC_paper/2020")
library(ggplot2)
DE <- read.table("cluster_memberships_kmeans_nclust6_final.txt", header=T, sep='\t')

# non-coding = yes includes 3prime_overlapping_ncrna, antisense, lincRNA, miRNA, processed_transcript, sense_overlapping, sense_intronic, misc_RNA
# non-coding = no includes protein_coding, pseudogene, IG_C_pseudogene, IG_V_gene
# non-coding = blank includes snRNA, rRNA, snoRNA (the last three are clearly defined ncRNAs with specific function, thus not of interest to this class)
nrow(DE) # 3509
nrow(DE[DE$noncoding=="yes",]) #812
nrow(DE[DE$noncoding=="no",]) #2596
nrow(DE[DE$noncoding=="",]) #101

nrow(DE[DE$new_clus==1 & DE$noncoding == "yes",]) #168
nrow(DE[DE$new_clus==1 & DE$noncoding == "no",]) #174
nrow(DE[DE$new_clus==1 & DE$noncoding == "",]) #24

nrow(DE[DE$new_clus==2 & DE$noncoding == "yes",]) #77
nrow(DE[DE$new_clus==2 & DE$noncoding == "no",]) #73
nrow(DE[DE$new_clus==2 & DE$noncoding == "",]) #12

nrow(DE[DE$new_clus==3 & DE$noncoding == "yes",]) #152
nrow(DE[DE$new_clus==3 & DE$noncoding == "no",]) #821
nrow(DE[DE$new_clus==3 & DE$noncoding == "",]) #29

nrow(DE[DE$new_clus==4 & DE$noncoding == "yes",]) #62
nrow(DE[DE$new_clus==4 & DE$noncoding == "no",]) #537
nrow(DE[DE$new_clus==4 & DE$noncoding == "",]) #10

nrow(DE[DE$new_clus==5 & DE$noncoding == "yes",]) #230
nrow(DE[DE$new_clus==5 & DE$noncoding == "no",]) #333
nrow(DE[DE$new_clus==5 & DE$noncoding == "",]) #25

nrow(DE[DE$new_clus==6 & DE$noncoding == "yes",]) #123
nrow(DE[DE$new_clus==6 & DE$noncoding == "no",]) #658
nrow(DE[DE$new_clus==6 & DE$noncoding == "",]) #1

mat_clus5 <- matrix(c(230,333,812,2596),nrow=2)
mat_clus4 <- matrix(c(62,537,812,2596),nrow=2)
mat_clus6 <- matrix(c(123,658,812,2596),nrow=2)
mat_clus1 <- matrix(c(168,174,812,2596),nrow=2)
mat_clus3 <- matrix(c(152,821,812,2596),nrow=2)
mat_clus2 <- matrix(c(77,73,812,2596),nrow=2)

fisher.test(mat_clus1)
fisher.test(mat_clus2)
fisher.test(mat_clus3)
fisher.test(mat_clus4)
fisher.test(mat_clus5)
fisher.test(mat_clus6)

fisher.test(mat_clus1)

# Fisher's Exact Test for Count Data
# 
# data:  mat_clus1
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.444867 3.894572
# sample estimates:
# odds ratio 
#   3.085678 

fisher.test(mat_clus2)
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  mat_clus2
# p-value = 1.427e-12
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#     2.390722 4.756557
# sample estimates:
#     odds ratio 
# 3.370762 

fisher.test(mat_clus3)

# Fisher's Exact Test for Count Data
# 
# data:  mat_clus3
# p-value = 2.342e-08
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.4860696 0.7177112
# sample estimates:
# odds ratio 
#  0.5919666 
# 
fisher.test(mat_clus4)

# 	Fisher's Exact Test for Count Data
# 
# data:  mat_clus4
# p-value = 4.049e-15
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#     0.2758763 0.4871169
# sample estimates:
#     odds ratio 
# 0.3691935 

fisher.test(mat_clus5)

# Fisher's Exact Test for Count Data
# 
# data:  mat_clus5
# p-value = 2.542e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.824745 2.668480
# sample estimates:
# odds ratio 
#   2.207681 

fisher.test(mat_clus6)

# 	Fisher's Exact Test for Count Data
# 
# data:  mat_clus6
# p-value = 5.334e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#     0.4812970 0.7379167
# sample estimates:
#     odds ratio 
# 0.5976909 


# make a plot

tmp1 <- fisher.test(mat_clus1)
tmp2 <- fisher.test(mat_clus2)
tmp3 <- fisher.test(mat_clus3)
tmp4 <- fisher.test(mat_clus4)
tmp5 <- fisher.test(mat_clus5)
tmp6 <- fisher.test(mat_clus6)

boxLabels <- c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5","cluster 6")
dat <- data.frame(cluster = 1:length(boxLabels), 
                  pval = c(-log10(tmp1$p.value),-log10(tmp2$p.value),-log10(tmp3$p.value),-log10(tmp4$p.value), -log10(tmp5$p.value), -log10(tmp6$p.value)),
                 OR = c(log10(tmp1$estimate),log10(tmp2$estimate),log10(tmp3$estimate),log10(tmp4$estimate), log10(tmp5$estimate), log10(tmp6$estimate)))

p <- ggplot(dat, aes(x = OR, y = boxLabels)) + 
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    geom_point(aes(colour = dat$pval), size = 4.5) +
    scale_x_continuous(breaks = (seq(-0.75, 0.75, 0.25)), labels = seq(-0.75, 0.75, 0.25),
                       limits = c(-0.54,0.54)) +
    scale_y_discrete(limits = rev(boxLabels)) +
                         
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position="right") +
    scale_colour_continuous("-log10(pval)") +
    ylab("") +
    xlab("log10 Odds ratio") +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"))

    #ggtitle("enrichment non-coding RNAs per cluster")



