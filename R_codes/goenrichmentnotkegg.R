#GO and KEGG functional enrichment analyses for significantly expanded families

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(viridis)

#Genes of interest:
eggnog_data <- read.csv("../eggnog/three_leafhopper/genes_interest.csv",
                        header=T)
genes <- eggnog_data$query


#All genes

eggnog_data2<- read.csv("../eggnog/three_leafhopper/MM_h8cwolle.emapper.annotations.csv",
                        header=T)

kegg_data <- eggnog_data2[c(1,8,10)]
kegg_data$GOs <- gsub( "GO:", "", as.character(kegg_data$GOs))
kegg_data <- data.table(kegg_data)
ko_data <- kegg_data[, list(GOs = unlist(strsplit(GOs , ","))), by = c("query", "Description")]

kegg_final <- ko_data[,c(3,1,2)]
kegg_description <-  ko_data[,c(3,2)]

enr_res <- enricher(genes, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10)
dotplot(enr_res, showCategory=10, font.size = 8, label_format = 50) + ggtitle("dotplot for ALF enrichment")
barplot(enr_res, showCategory=10, font.size = 8, label_format = 50) + ggtitle("barplot for ALF enrichment")



