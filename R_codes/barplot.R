library(tidyverse)
library(ggplot2)
library(svglite)

df <- read.csv("../dna_evorates/three_leafhopper/barplot.csv", header = T)



df2<- df%>% pivot_longer(cols = !Species, names_to = "status", values_to = "percentage")
df2$status[df2$status == 'X1.1.1'] <- '1:1:1'
df2$status[df2$status == 'Unassigned.Genes'] <- 'Unassigned'
df2$status[df2$status == 'Species.specific'] <- 'Species-specific'
df2$status[df2$status == 'N.N.N'] <- 'N:N:N'
df2$status[df2$status == 'Others'] <- 'Other'


df2$status <- factor(df2$status, levels=c( 'Unassigned', 'Other','Species-specific','N:N:N','1:1:1'))

df3<- df2 %>%
  arrange(percentage) %>%
  mutate(Species = factor(Species, levels=c("Homalodisca vitripennis", "Nephotettix cincticeps", "Macrosteles quadrilineatus" ))) %>%
ggplot( aes(fill = status, x = Species, y = percentage)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.3) + 
  scale_fill_brewer(palette = "Paired", direction=-1) +
  coord_flip() + theme(legend.title=element_blank()) + theme_bw(base_size = 16) +
  guides(fill= guide_legend(reverse=TRUE))  +
  ylab("Number of genes") + theme(plot.margin = margin(100, 30, 100, 30))


ggsave(file="../dna_evorates/three_leafhopper/barplot2.svg", plot=df3, width=15, height=10, dpi = 300)

