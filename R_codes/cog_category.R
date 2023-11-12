library(tidyverse)
library(dplyr)

t1 <- read.csv("../four/cog/four_leafcog.csv")
t2 <- read.csv("../four/cog/two_leafcog.csv")
t3 <- read.csv("../four/cog/four_cog_all.csv")
t4 <- read.csv("../four/cog/two_cog_all.csv")



expanded_df <- t1 %>% separate_rows(COG.category, sep = "")

expanded_df2 <- t2 %>% separate_rows(COG.category, sep = "")

expanded_df3 <- t3 %>% separate_rows(COG.category, sep = "")

expanded_df4 <- t4 %>% separate_rows(COG.category, sep = "")

string_counts <- as.data.frame(table(expanded_df$COG.category))
string_counts2 <- as.data.frame(table(expanded_df2$COG.category))
string_counts3 <- as.data.frame(table(expanded_df3$COG.category, expanded_df3$State))
string_counts4 <- as.data.frame(table(expanded_df4$COG.category, expanded_df4$State))

filtered_values <- string_counts[string_counts$Freq < 1000, ]
filtered_values2 <- string_counts2[string_counts2$Freq < 1000, ]
filtered_values3 <- string_counts3[string_counts3$Freq < 1000, ]
filtered_values4 <- string_counts4[string_counts4$Freq < 1000, ]

ggplot(data=filtered_values, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="#377EB8")+
  theme_minimal() + labs(x = "COG Category", y = "Count", title = "Positive Selection in COG Categories",
                         subtitle = "Comparison between four leafhoppers") 


ggplot(data=filtered_values2, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="#377EB8")+
  theme_minimal() + labs(x = "COG Category", y = "Count", title = "Positive Selection in COG Categories",
                         subtitle = "Comparison between two leafhoppers") 


p1 <- ggplot(data=filtered_values3, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal() + labs(x = "COG Category", y = "Count", 
                         title = "Selection in COG Categories",
                         subtitle = "Comparison between four leafhoppers",
                         fill = "Selection") + 
  scale_fill_brewer(palette = "Set1")


p2 <- ggplot(data=filtered_values4, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal() + labs(x = "COG Category", y = "Count", 
                         title = "Selection in COG Categories",
                         subtitle = "Comparison between two leafhoppers",
                         fill = "Selection") + 
  scale_fill_brewer(palette = "Set1")

ggsave(file="../four/cog/four.svg", plot=p1, width=10, height=10, dpi=300)
ggsave(file="../four/cog/two.svg", plot=p2, width=10, height=10, dpi=300)



expanded_df3 <- t3 %>% separate_rows(COG.category, sep = "")
expanded_df3$COG.category <- na_if(expanded_df3$COG.category, '')
expanded_df3 <- na.omit(expanded_df3)

string_counts3 <- as.data.frame(table(expanded_df3$COG.category, expanded_df3$State))
positive_selection_df <- subset(string_counts3, Var2 == "Positive Selection")
sum_freq_positive_selection <- sum(positive_selection_df$Freq)
negative_selection_df <- subset(string_counts3, Var2 == "Negative Neutral Selection")
sum_freq_negative_selection <- sum(negative_selection_df$Freq)






results <- list()

# Get unique COG categories
unique_categories <- unique(expanded_df3$COG.category)

for (category in unique_categories) {
  # Subset the data for the current COG category
  subset_data <- subset(expanded_df3, COG.category == category)
  
  # Create a contingency table
  contingency_table <- table(subset_data$COG.category,subset_data$State)
  temp <- as.data.frame(contingency_table)
  neg_sum <- 2141 - temp$Freq[1]
  pos_sum <- 1062 - temp$Freq[2]
  contingency_table <- rbind(contingency_table, c(neg_sum, pos_sum)) 
  contingency_table <- contingency_table[, ncol(contingency_table):1]
  
  # Perform Fisher's exact test
  fisher_result <- chisq.test(contingency_table, simulate.p.value = TRUE)
  
  # Store the contingency table and Fisher's exact test result
  results[[category]] <- list(contingency_table = contingency_table, fisher_result = fisher_result)
}



for (category in unique_categories) {
  cat("Contingency table for COG category", category, ":\n")
  print(results[[category]]$contingency_table)
  
  cat("Fisher's exact test result for COG category", category, ":\n")
  print(results[[category]]$fisher_result)
  cat("\n")
}













# Perform Fisher's exact test for each COG category
results <- lapply(rownames(contingency_table), function(category) {
  subset_table <- contingency_table[category, ]
  fisher.test(subset_table)
})

# Extract p-values and determine enrichment
enrichment <- sapply(results, function(x) {
  p_value <- x$p.value
  enrichment <- ifelse(p_value < 0.05, "Enriched", "Not Enriched")
  list(p_value = p_value, enrichment = enrichment)
})








results <- by(expanded_df3, expanded_df3$COG.category, function(x) {
  table <- table(x$State)
  chisq.test(table)
})

enrichment <- sapply(results, function(x) {
  if (length(x) > 0) {
    p_value <- x$p.value
    enrichment <- ifelse(p_value < 0.05, "Enriched", "Not Enriched")
    list(p_value = p_value, enrichment = enrichment)
  } else {
    list(p_value = NA, enrichment = NA)
  }
})



# Convert results to a dataframe
enrichment_df <- as.data.frame(t(enrichment))
 
enrichment_df<- filtered_values3 %>% left_join(rownames_to_column(enrichment_df), by=c("Var1" = "rowname"))


# Plot the p-values using ggplot2
ggplot(enrichment_df, aes(x = Var1, y = Freq, fill=Var2)) +
  geom_bar(stat = "identity",position=position_dodge() ) +
  geom_text(aes(label = ifelse(p_value < 0.05, "*", "")),  color = "black", size = 3.5) +
  labs(x = "Category", y = "Count", title = "Enrichment Analysis") +
  theme_minimal()




temp_total <- as.data.frame(table(expanded_df3$COG.category))
tt <- temp_total[temp_total$Freq < 1300, ]

math <- filtered_values %>% left_join(tt,by=c("Var1" = "Var1"))
math$percent <- (math$Freq.x/math$Freq.y)



temp_total <- as.data.frame(table(expanded_df4$COG.category))
tt <- temp_total[temp_total$Freq < 1300, ]

math2 <- filtered_values2 %>% left_join(tt,by=c("Var1" = "Var1"))
math2$percent <- (math2$Freq.x/math2$Freq.y)
