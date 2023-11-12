#GO and KEGG functional enrichment analyses for significantly expanded families

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(viridis)
library("stringr") 

#Genes of interest:



eggnog_data <- read.csv("../redo_noiso/eggnog/similarities3",
                        header=F)

genes <- eggnog_data$V1

#genes <- head(genes)

#All genes

eggnog_data <- read.csv("../redo_noiso/eggnog/MM_zbuaoxef.emapper.annotations.tsv",
                        sep="\t",skip=4, header=T)


#get columns 
kegg_data <- eggnog_data[c(1,8,12)]

# clean up by removing the "ko:" 
kegg_data$KEGG_ko <- gsub( "ko:", "", as.character(kegg_data$KEGG_ko))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_data <- data.table(kegg_data)

ko_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = c("X.query", "Description")]
#go_data <- kegg_data[, list(GOs = unlist(strsplit(GOs , ","))), by = X.query]

#counts <- go_data %>% 
#  count(go_data$GOs)

#counts <- counts[with(counts,order(-counts$n)),]
#first_ten <- counts[1:10,]



#counts_kegg <- ko_data %>% 
#  count(ko_data$KEGG_ko)

#counts_kegg <- counts_kegg[with(counts_kegg,order(-counts_kegg$n)),]
#first_ten_kegg <- counts_kegg[1:11,]





#K04834: voltage-gated sodium channel type II alpha
#K04841: voltage-gated sodium channel type IX alpha
#K06767: Down syndrome cell adhesion molecule
#K12567: TTN; titin [EC:2.7.11.1]
#K08145: SLC2A8, GLUT8; MFS transporter, SP family, 
#solute carrier family 2 (facilitated glucose transporter), member 8
#K17591: RIMBP2; RIMS-binding protein 2
#K19922: TSPOAP1, BZRAP1, RIMBP1; peripheral-type benzodiazepine receptor-associated protein 1
#K02177: BRU; bruno
#K13207: CUGBP, BRUNOL, CELF; CUG-BP- and ETR3-like factor
#K10415: DYNC1I, DNCI; dynein cytoplasmic 1 intermediate chain



#ko_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = X.query]

kegg_final <- ko_data[,c(3,1,2)]
kegg_description <-  ko_data[,c(3,2)]

enr_res <- enricher(genes, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 10)
dotplot(enr_res, showCategory=10)

length(unique(kegg_final$X.query))

#splitting by expansions and contracts


expcont <- read.csv("../redo_noiso/eggnog/macexpcont2.csv", header = F)
exp <- expcont[which(expcont$V2=='Expansion'),]
cont <- expcont[which(expcont$V2=='Contraction'),]


genes_exp <- exp$V1
genes_cont <- cont$V1

enr_res_exp <- enricher(genes_exp, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 10)
dotplot(enr_res_exp, showCategory=10)

enr_res_cont <- enricher(genes_cont, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 2)
dotplot(enr_res_cont, showCategory=10)

png("../redo_noiso/eggnog/enriched_expanded.png", width = 3065, height = 2025, res = 300)
dotplot(enr_res_exp, showCategory=10)
dev.off()



t<- enr_res_exp@result
t2<- enr_res_cont@result

write.csv(t, "../redo_noiso/eggnog/contracted_results2.csv")

write.csv(t2, "../redo_noiso/eggnog/expanded_results2.csv")



f <- read.csv(file.choose())







cont_df <- read.csv("../eggnog/hemiptera/output_contractions.txt", header = FALSE)
cont_df$label <- "contraction"
df_cont<-left_join(kegg_final, cont_df, by=c('X.query'='V1'))

exps_df <- read.csv("../eggnog/hemiptera/output_expansions.txt", header = FALSE)
exps_df$label <- "expansion"
df_exp<-left_join(kegg_final, exps_df, by=c('X.query'='V1'))


kegg_final2 <- inner_join(df_cont, df_exp, by=c('KEGG_ko'='KEGG_ko', 'X.query'='X.query'))
kegg_final2$merged <- coalesce(kegg_final2$label.x, kegg_final2$label.y)
kegg_final3 <- kegg_final2[,c(1,2,3,7 )]


kegg_final3[is.na(kegg_final3)] <- "nothing"
expansions_kegg <- kegg_final3[!kegg_final3$merged == "contraction", ]
contractions_kegg <- kegg_final3[!kegg_final3$merged == "expansion", ]

kegg_description_exp <-  expansions_kegg[,c(1,3)]
kegg_description_con <-  contractions_kegg[,c(1,3)]

enr_res_exp <- enricher(genes, TERM2GENE=expansions_kegg, TERM2NAME = kegg_description_exp, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10)

enr_res_con <- enricher(genes, TERM2GENE=contractions_kegg, TERM2NAME = kegg_description_con, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10)


dotplot(enr_res_exp, showCategory=13, title = "Significantly expanded gene families in ALF") + scale_color_gradient(low = "#56B1F7", high = "#132B43")
dotplot(enr_res_con, showCategory=10, title = "Significantly contracted gene families in ALF") + scale_color_gradient(low = "#56B1F7", high = "#132B43")


p1 <- barplot(enr_res_exp, showCategory=13, font.size = 8, label_format = 50) + 
  ggtitle("Top 10 Significantly expanded gene families in ALF")+
  scale_fill_gradient(low = "#56B1F7", high = "#56B1F7")

p2 <- barplot(enr_res_con, showCategory=40, font.size = 8, label_format = 50) + 
  ggtitle("Significantly contracted gene families in ALF")+
  scale_fill_gradient(low = "#56B1F7", high = "#56B1F7")




df_top <- read.csv("../redo_noiso/eggnog/barplot_expcont.csv")

#df_top$X.1 <- factor(df_top$X.1, levels=c('Expansion', 'Contraction'))



p1<- df_top %>% arrange(Category, Count) %>% 
  mutate(Description = factor(Description, unique(Description))) %>%
  ggplot(aes(x=Description, y=Count, fill = Category)) +
  geom_bar(stat="identity" , position='dodge' ) + coord_flip() + 
  scale_fill_manual(values = c("#A3B18A", "#3A5A40"))  +        
  scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
  xlab("") + ylab("Gene Count") + theme(legend.title=element_blank()) + theme_bw(base_size = 10) +
  guides(fill= guide_legend(reverse=TRUE))  + theme(axis.title.y = element_text( color="black", 
                                                                                 face="bold",
                                                                                 angle=0))
ggsave(file="expansinos_contractions_barplot2.svg", plot=p1, width=15, height=10, dpi = 300)
















eggnog_data <- read.csv("../redo_noiso/eggnog/MM_zbuaoxef.emapper.annotations.tsv",
                        sep="\t",skip=4, header=T)


#get columns 
kegg_data <- eggnog_data[c(1,8,10)]

# clean up by removing the "ko:" 
kegg_data$GOs<- gsub( "GO:", "", as.character(kegg_data$GOs))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_data <- data.table(kegg_data)

ko_data <- kegg_data[, list(GOs = unlist(strsplit(GOs , ","))), by = c("X.query", "Description")]
#go_data <- kegg_data[, list(GOs = unlist(strsplit(GOs , ","))), by = X.query]

#counts <- go_data %>% 
#  count(go_data$GOs)

#counts <- counts[with(counts,order(-counts$n)),]
#first_ten <- counts[1:10,]



#counts_kegg <- ko_data %>% 
#  count(ko_data$KEGG_ko)

#counts_kegg <- counts_kegg[with(counts_kegg,order(-counts_kegg$n)),]
#first_ten_kegg <- counts_kegg[1:11,]





#K04834: voltage-gated sodium channel type II alpha
#K04841: voltage-gated sodium channel type IX alpha
#K06767: Down syndrome cell adhesion molecule
#K12567: TTN; titin [EC:2.7.11.1]
#K08145: SLC2A8, GLUT8; MFS transporter, SP family, 
#solute carrier family 2 (facilitated glucose transporter), member 8
#K17591: RIMBP2; RIMS-binding protein 2
#K19922: TSPOAP1, BZRAP1, RIMBP1; peripheral-type benzodiazepine receptor-associated protein 1
#K02177: BRU; bruno
#K13207: CUGBP, BRUNOL, CELF; CUG-BP- and ETR3-like factor
#K10415: DYNC1I, DNCI; dynein cytoplasmic 1 intermediate chain



#ko_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = X.query]

kegg_final <- ko_data[,c(3,1,2)]
kegg_description <-  ko_data[,c(3,2)]

enr_res <- enricher(genes_exp, TERM2GENE=kegg_final, TERM2NAME = kegg_description, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 3)
dotplot(enr_res, showCategory=20)




