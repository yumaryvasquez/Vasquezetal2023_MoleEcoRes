library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(viridis)
library("stringr") 


#Genes of interest:

eggnog_data <- read.csv("../positive_selection_enrichment/four_positive.csv",
                        header=T)

genes <- eggnog_data$Ortho.Group.Label


#All genes

eggnog_data <- read.csv("../positive_selection_enrichment/four_all.csv",
                       header=T)


#get columns 
kegg_data <- eggnog_data[c(1,4,6)]

# clean up by removing the "ko:" 
kegg_data$KEGG_ko <- gsub( "ko:", "", as.character(kegg_data$KEGG.ko.Categories))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_data <- data.table(kegg_data)

ko_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = c("Ortho.Group.Label", "KEGG.Description")]


kegg_final <- ko_data[,c(3,1,2)]
kegg_description <-  ko_data[,c(3,2)]

enr_res <- enricher(genes, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 2)
dotplot(enr_res, showCategory=10)



enr_res <- enricher(genes, TERM2GENE=kegg_final, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 2)
dotplot(enr_res, showCategory=10) + scale_color_viridis(option="plasma")  + theme_bw()














#Genes of interest:

eggnog_data <- read.csv("../positive_selection_enrichment/four_positive.csv",
                        header=T)

genes <- eggnog_data$Ortho.Group.Label


#All genes

eggnog_data <- read.csv("../positive_selection_enrichment/four_all.csv",
                        header=T)


#get columns 
kegg_data <- eggnog_data[c(1,4,5)]

# clean up by removing the "ko:" 
kegg_data$KEGG_ko <- gsub( "GO:", "", as.character(kegg_data$GO.Categories))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_data <- data.table(kegg_data)

ko_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = c("Ortho.Group.Label", "KEGG.Description")]


kegg_final <- ko_data[,c(3,1,2)]
#kegg_description <-  ko_data[,c(3,2)]

#enr_res <- enricher(genes, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 10)
#dotplot(enr_res, showCategory=10)

enr_res <- enricher(genes, TERM2GENE=kegg_final, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 10)
p1<- dotplot(enr_res, showCategory=10) + scale_color_viridis(option="plasma")  + theme_bw()





#TWO#####


#Genes of interest:

eggnog_data <- read.csv("../positive_selection_enrichment/two_positive.csv",
                        header=T)

genes <- eggnog_data$Ortho.Group.Label


#All genes

eggnog_data <- read.csv("../positive_selection_enrichment/two_all.csv",
                        header=T)


#get columns 
kegg_data <- eggnog_data[c(1,4,6)]

# clean up by removing the "ko:" 
kegg_data$KEGG_ko <- gsub( "ko:", "", as.character(kegg_data$KEGG.ko.Categories))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_data <- data.table(kegg_data)

ko_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = c("Ortho.Group.Label", "KEGG.Description")]


kegg_final <- ko_data[,c(3,1,2)]
kegg_description <-  ko_data[,c(3,2)]

enr_res <- enricher(genes, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 2)
dotplot(enr_res, showCategory=10)


















#Genes of interest:

eggnog_data <- read.csv("../positive_selection_enrichment/two_positive.csv",
                        header=T)

genes <- eggnog_data$Ortho.Group.Label


#All genes

eggnog_data <- read.csv("../positive_selection_enrichment/two_all.csv",
                        header=T)


#get columns 
kegg_data <- eggnog_data[c(1,4,5)]

# clean up by removing the "ko:" 
kegg_data$KEGG_ko <- gsub( "GO:", "", as.character(kegg_data$GO.Categories))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_data <- data.table(kegg_data)

ko_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = c("Ortho.Group.Label", "KEGG.Description")]


kegg_final <- ko_data[,c(3,1,2)]
#kegg_description <-  ko_data[,c(3,2)]

#enr_res <- enricher(genes, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 10)
#dotplot(enr_res, showCategory=10)

enr_res <- enricher(genes, TERM2GENE=kegg_final, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 10)
p2 <- dotplot(enr_res, showCategory=10) + scale_color_viridis(option="plasma")  + theme_bw()


ggsave(file="../four_enrichment.svg", plot=p1, width=5, height=5, dpi=300)
ggsave(file="../two_enrichment.svg", plot=p2, width=5, height=5, dpi=300)



