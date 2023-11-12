library(tidyverse)



likeRT <- function(x, y) {
  -2 * (x-y) 
}





island <- read.csv("../dna_evorates/three_leafhopper/M7_M8.csv")

likelihoodratio_n <- c()
for (i in 1:nrow(island)) {
  likelihoodratio_n <- c(likelihoodratio_n, 
                         likeRT(island$M7[i], 
                                island$M8[i]))
}

island <- island %>%
  mutate(likelihoodratio = likelihoodratio_n)

island$likelihoodratio <- likelihoodratio_n



chi_pval_n <- c()
for (i in 1:nrow(island)) {
  chi_pval_n <- c(chi_pval_n,
                  pchisq(island$likelihoodratio[i],
                         df=2, 
                         lower.tail = FALSE))
}

island <- island %>% 
  mutate(Chi_pvalue = chi_pval_n)


island$Chi_pvalue <- chi_pval_n


write.csv(island, "../dna_evorates/three_leafhopper/M7_M8.csv")

