df <- read.csv("../eggnog/hemiptera/bioprocess.csv")


df$Label <- factor(df$Label, levels=c( 'alf', 'grl','gwss','Shared_alf_grl','Shared_alf_gwss','Shared_gwss_grl','Shared_between_all'))




df2 <- df[which(df$Count > 40) ,]


ggplot(df2, aes(x=reorder(Description, Count), y=Count, fill = Label)) +
  geom_bar(stat="identity" , position='stack') +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() + coord_flip()..



df <- read.csv("../eggnog/hemiptera/enrichment.csv")





p1<-ggplot(df, aes(x=reorder(Description, Count), y=Count, fill = Type)) +
  geom_bar(stat="identity" , position='dodge') +
  scale_fill_manual(values = c("#023047", "#8ECAE6"))+
  theme_bw() + facet_wrap(~Label, scales="free", ncol=2) + coord_flip()


ggsave(file="enrichedall_barplot.svg", plot=p1, width=15, height=10, dpi = 300)











x = reorder(f.name, -age)

p1<- df_top %>% arrange(X.1, Count) %>% 
  mutate(Description = factor(Description, unique(Description))) %>%
  ggplot(aes(x=Description, y=Count, fill = X.1)) +
  geom_bar(stat="identity" , position='dodge' ) + coord_flip() + 
  scale_fill_manual(values = c("#A3B18A", "#3A5A40"))  +        
  scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
  xlab("") + ylab("Gene Count") + theme(legend.title=element_blank()) + theme_bw(base_size = 10) +
  guides(fill= guide_legend(reverse=TRUE))  + theme(axis.title.y = element_text( color="black", 
                                                                                 face="bold",
                                                                                 angle=0))
ggsave(file="expansinos_contractions_barplot.svg", plot=p1, width=15, height=10, dpi = 300)