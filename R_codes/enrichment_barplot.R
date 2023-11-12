df <- read.csv("orthovenn_results.csv")

split_df <- split(df, df$Cluster)
sub_df_A <- split_df$Shared
sub_df_B <- split_df$`Sulcia-associated`
sub_df_C <- split_df$Deltocephalinae
sub_df_D <- split_df$ALF


p1 <- ggplot(sub_df_A, aes(x=reorder(Name,Count), y=Count)) +
  geom_bar(stat="identity", fill="#5e8072", width=0.5) + 
  coord_flip() + theme_minimal()

p2 <-ggplot(sub_df_B, aes(x=reorder(Name,Count), y=Count)) +
  geom_bar(stat="identity", fill="#638076", width=0.5) + 
  coord_flip() + theme_minimal() + scale_y_continuous(breaks = c(2, 4, 6))
  
  
p3 <- ggplot(sub_df_C, aes(x=reorder(Name,Count), y=Count)) +
  geom_bar(stat="identity", fill="#8d4948", width=0.5) + 
  coord_flip() + theme_minimal()
  
  
p4 <-ggplot(sub_df_D, aes(x=reorder(Name,Count), y=Count)) +
  geom_bar(stat="identity", fill="#872254", width=0.5) + 
  coord_flip() + theme_minimal()


ggsave(p1,filename ="Shared_orthovenn.svg", width = 10, height = 8)

ggsave(p2,filename ="Sulcia_orthovenn.svg", width = 10, height = 8)

ggsave(p3,filename ="Delto_orthovenn.svg", width = 10, height = 8)

ggsave(p4,filename ="ALF_orthovenn.svg", width = 10, height = 8)

