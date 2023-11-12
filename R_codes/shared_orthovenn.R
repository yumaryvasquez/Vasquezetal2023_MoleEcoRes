#orthovenn -> graphics

library(dplyr)
library(tidyverse)
library(ggbreak) 
library(patchwork)


#ALF
bio_df <- read.csv("../OrthoVenn/MAC_go/bio_process_file.txt",
                   sep="\t", header=F)
cell_df <- read.csv("../OrthoVenn/MAC_go/cell_component_file.txt",
                    sep="\t", header=F)
molec_df <- read.csv("../OrthoVenn/MAC_go/molecular_function_file.txt",
                     sep="\t", header=F)
bio_df$label <- "Biological process"
cell_df$label <- "Cellular component"
molec_df$label <- "Molecular function"
df_alf <- tibble(rbind(bio_df, cell_df, molec_df))
df_alf$species <- "ALF"
total = sum(df_alf$V3)


#GWSS
bio_df <- read.csv("../OrthoVenn/Homo_go/bio_process_file.txt",
                   sep="\t", header=F)
cell_df <- read.csv("../OrthoVenn/Homo_go/cell_component_file.txt",
                    sep="\t", header=F)
molec_df <- read.csv("../OrthoVenn/Homo_go/molecular_function_file.txt",
                     sep="\t", header=F)
bio_df$label <- "Biological process"
cell_df$label <- "Cellular component"
molec_df$label <- "Molecular function"
df_gwss <- tibble(rbind(bio_df, cell_df, molec_df))
df_gwss$species <- "GWSS"

#RGL
bio_df <- read.csv("../OrthoVenn/Nepho_go/bio_process_file.txt",
                   sep="\t", header=F)
cell_df <- read.csv("../OrthoVenn/Nepho_go/cell_component_file.txt",
                    sep="\t", header=F)
molec_df <- read.csv("../OrthoVenn/Nepho_go/molecular_function_file.txt",
                     sep="\t", header=F)
bio_df$label <- "Biological process"
cell_df$label <- "Cellular component"
molec_df$label <- "Molecular function"
df_rgl <- tibble(rbind(bio_df, cell_df, molec_df))
df_rgl$species <- "RGL"
total = sum(df_rgl$V3)



#all
df_all <- tibble(rbind(df_alf, df_gwss, df_rgl))

df_all %>% arrange(label, desc(V3)) %>% 
  mutate(V2 = factor(V2, unique(V2))) %>%
ggplot(aes(x=V2, y=V3, fill=species)) +
  geom_bar(position="stack", stat="identity" ) + coord_flip() +
  theme_minimal() +scale_fill_brewer(palette="Dark2") +
  facet_grid(~label, scales = "free_x") 



ggplot(df_all, aes(x=species, y=V3, fill=label)) +
  geom_bar(position="stack", stat="identity" ) + coord_flip() +
  theme_minimal() +scale_fill_brewer(palette="Dark2") 




#shared genes

bio_df <- read.csv("../OrthoVenn/Shared_all/bio_process_file.txt",
                   sep="\t", header=F)
cell_df <- read.csv("../OrthoVenn/Shared_all/cell_component_file.txt",
                    sep="\t", header=F)
molec_df <- read.csv("../OrthoVenn/Shared_all/molecular_function_file.txt",
                     sep="\t", header=F)
bio_df$label <- "Biological process"
cell_df$label <- "Cellular component"
molec_df$label <- "Molecular function"
bio_df2 <- bio_df[which((bio_df$V3 > 70) & (bio_df$V3 < 1000)) ,]
cell_df2 <- cell_df[which(cell_df$V3 > 3),]
molec_df2 <- molec_df[which(molec_df$V3 > 3),]

df <- tibble(rbind(bio_df2, cell_df2, molec_df2))

df %>% arrange(label, desc(V3)) %>% 
  mutate(V2 = factor(V2, unique(V2))) %>%
  ggplot( aes(x=V2, y=V3, fill=label)) +
  geom_bar(stat="identity" , position='dodge' ) + coord_flip() +
  theme_minimal() +scale_fill_brewer(palette="Dark2") +
  labs(title="Shared Genes",
       x ="Number of genes", y = "Functional categories")+ 
  theme(legend.title = element_blank(),text = element_text(size=14) )

