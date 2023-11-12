#GO and KEGG functional enrichment analyses for significantly expanded families

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(data.table)
library(dplyr)


eggnog_data <- read.csv("../eggnog/MM_d1fzrebj.emapper.annotations.tsv",
                sep="\t",skip=4, header=T)



#get columns 
kegg_data <- eggnog_data[c(1,7,10,12)]

# clean up by removing the "ko:" 
kegg_data$KEGG_ko <- gsub( "ko:", "", as.character(kegg_data$KEGG_ko))

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_data <- data.table(kegg_data)

ko_data <- kegg_data[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = X.query]
go_data <- kegg_data[, list(GOs = unlist(strsplit(GOs , ","))), by = X.query]

counts <- go_data %>% 
          count(go_data$GOs)

counts <- counts[with(counts,order(-counts$n)),]
first_ten <- counts[1:10,]



counts_kegg <- ko_data %>% 
  count(ko_data$KEGG_ko)

counts_kegg <- counts_kegg[with(counts_kegg,order(-counts_kegg$n)),]
first_ten_kegg <- counts_kegg[1:11,]





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






