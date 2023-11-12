library(dplyr)

df3 <- read.csv("../Cafe/results/Base_branch_probabilities.tab",
                sep="\t", header=T)

new_df = subset(df3, select = c(X.FamilyID, Macrosteles_quadrilineatus.2.))


#new_df$pvalue0.05 <- with(new_df,new_df$Macrosteles_quadrilineatus.2.<0.05)
new_df$pvalue0.01 <- with(new_df,new_df$Macrosteles_quadrilineatus.2.<0.01)
#new_df$pvalue0.001 <- with(new_df,new_df$Macrosteles_quadrilineatus.2.<0.001)
summary(new_df$pvalue0.05)
  #Mode   FALSE    TRUE 
  #logical     324     255 
summary(new_df$pvalue0.01)
  #Mode   FALSE    TRUE 
  #logical     409     170 
summary(new_df$pvalue0.001)
  #Mode   FALSE    TRUE 
  #logical     467     112 

df <- new_df[new_df$pvalue0.01 == "TRUE", ]   

df2 <- read.csv("../Cafe/results/Base_change.tab",
                sep="\t", header=T)
new_df = subset(df2, select = c(FamilyID, Macrosteles_quadrilineatus.2.))

df3<-left_join(df, new_df, by=c('X.FamilyID'='FamilyID'))
df3$expcont <- ifelse(df3$Macrosteles_quadrilineatus.2..y > 0, "Expansion", "Contraction")


write.csv(df3, "../Cafe/mac.csv")
