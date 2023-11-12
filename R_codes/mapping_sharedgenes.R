te <- read.csv("../redo_noiso/genomes/new/shared_results/shared_only.out2", header = F, sep ="\t")
y2<-read.csv("../redo_noiso/eggnog/trinity_files/t1_1.csv", skip=1)
y3<-read.csv("../redo_noiso/eggnog/trinity_files/t2_1.csv", skip=1)
y4<-read.csv("../redo_noiso/eggnog/trinity_files/t3_1.csv", skip=1)
y5<-read.csv("../redo_noiso/eggnog/trinity_files/t4_1.csv", skip=1)
y6<-read.csv("../redo_noiso/eggnog/trinity_files/t5_1.csv", skip=1)
y7<-read.csv("../redo_noiso/eggnog/trinity_files/t6_1.csv", skip=1)
pe <-  read.csv("../redo_noiso/genomes/new/shared_results/all_labels.txt", header = F, sep ="\t")


te_new <- te %>% group_by(V2) %>% slice_max(V12) %>% ungroup()
te_new2 <-te_new[!duplicated(te_new$V2), ]


large_df <- bind_rows(y2, y3, y4, y5, y6, y7)
clean_df <- large_df %>% distinct(id, .keep_all = TRUE)


te_3 <- left_join(te_new2,clean_df,  by = c("V1"="id"))
te_3$body_avg <- rowMeans(te_3[, 17:19])
te_3$nasuia_avg <- rowMeans(te_3[, 20:22])
te_3$sulcia_avg <- rowMeans(te_3[, 23:25])

te_3$max_col <- apply(te_3[, c("body_avg", "nasuia_avg", "sulcia_avg")], 1, function(x) names(x)[which.max(x)])

te_4 <- te_3[!is.na(te_3$sulcia_avg),]


te_5 <- left_join(te_4,pe,  by = c("V2"="V2"))


te_6 <- apply(te_5,2,as.character)
write.csv(te_6, "../redo_noiso/genomes/new/shared_results/shared_combined.csv")


g <- read.csv("../redo_noiso/genomes/new/shared_results/shared_barplot.csv")

mdfr <- melt(g, id.vars = "Description")



p1 <- mdfr %>% arrange(value) %>% 
  mutate(Description = factor(Description, unique(Description))) %>%
  ggplot( aes(x=reorder(Description, +value), y = value, fill =variable)) + 
  geom_bar(stat="identity" , position='dodge' ) +
  coord_flip()+ 
  scale_fill_manual(values = c("#564D65", "#3E8989",  "#D77A61", "#D4E6B5")) + 
  theme_bw(base_size = 10) + facet_wrap(~Description, scales = "free")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(file="../redo_noiso/genomes/new/shared_results/shared_trinity_barplot.svg", plot=p1, width=8, height=10, dpi = 300)


mdfr %>% arrange(value) %>% 
  mutate(Description = factor(Description, unique(Description))) %>%
  ggplot( aes(x=reorder(Description, +value), y = value, fill =variable)) +
  geom_bar(stat="identity" , position='stack' ) + coord_flip()+ 
  scale_fill_brewer(palette = "Set2")

p1 <- mdfr %>% arrange(value) %>% 
  mutate(Description = factor(Description, unique(Description))) %>%
  ggplot( aes(x=reorder(Description, +value), y = value, fill =variable)) +
  geom_bar(stat="identity" , position='fill' ) + coord_flip()+ 
  scale_fill_manual(values = c("#564D65", "#3E8989",  "#D77A61", "#D4E6B5")) + 
  theme_bw(base_size = 10)





ggsave(file="../redo_noiso/trinity/new/expansions_trinity_barplot.svg", plot=p1, width=15, height=10, dpi = 300)


#------------------------------------#

te <- read.csv("../redo_noiso/genomes/new/mac_results/mac_only.out2", header = F, sep ="\t")
y2<-read.csv("../redo_noiso/eggnog/trinity_files/t1_1.csv", skip=1)
y3<-read.csv("../redo_noiso/eggnog/trinity_files/t2_1.csv", skip=1)
y4<-read.csv("../redo_noiso/eggnog/trinity_files/t3_1.csv", skip=1)
y5<-read.csv("../redo_noiso/eggnog/trinity_files/t4_1.csv", skip=1)
y6<-read.csv("../redo_noiso/eggnog/trinity_files/t5_1.csv", skip=1)
y7<-read.csv("../redo_noiso/eggnog/trinity_files/t6_1.csv", skip=1)
pe <-  read.csv("../redo_noiso/genomes/new/mac_results/mac_label.txt", header = F, sep ="\t")


te_new <- te %>% group_by(V2) %>% slice_max(V12) %>% ungroup()
te_new2 <-te_new[!duplicated(te_new$V2), ]


large_df <- bind_rows(y2, y3, y4, y5, y6, y7)
clean_df <- large_df %>% distinct(id, .keep_all = TRUE)


te_3 <- left_join(te_new2,clean_df,  by = c("V1"="id"))
te_3$body_avg <- rowMeans(te_3[, 17:19])
te_3$nasuia_avg <- rowMeans(te_3[, 20:22])
te_3$sulcia_avg <- rowMeans(te_3[, 23:25])

te_3$max_col <- apply(te_3[, c("body_avg", "nasuia_avg", "sulcia_avg")], 1, function(x) names(x)[which.max(x)])

te_4 <- te_3[!is.na(te_3$sulcia_avg),]


te_5 <- left_join(te_4,pe,  by = c("V2"="V2"))


te_6 <- apply(te_5,2,as.character)
write.csv(te_6, "../redo_noiso/genomes/new/mac_results/mac_combined.csv")


g <- read.csv("../redo_noiso/genomes/new/mac_results/mac_barplot.csv")

mdfr <- melt(g, id.vars = "Description")



p1 <- mdfr %>% arrange(value) %>% 
  mutate(Description = factor(Description, unique(Description))) %>%
  ggplot( aes(x=reorder(Description, +value), y = value, fill =variable)) + 
  geom_bar(stat="identity" , position='dodge' ) +
  coord_flip()+ 
  scale_fill_manual(values = c("#564D65", "#3E8989",  "#D77A61", "#D4E6B5")) + 
  theme_bw(base_size = 10) + facet_wrap(~Description, scales = "free")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(file="../redo_noiso/genomes/new/mac_results/mac_trinity_barplot.svg", plot=p1, width=8, height=10, dpi = 300)



