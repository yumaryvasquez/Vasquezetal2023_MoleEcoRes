island <- read.csv("../four/codeml/two/Book4.csv")

likelihoodratio_n_M7M8 <- c()
for (i in 1:nrow(island)) {
  likelihoodratio_n_M7M8 <- c(likelihoodratio_n_M7M8, 
                         likeRT(island$M7[i], 
                                island$M8[i]))
}

island <- island %>%
  mutate(likelihoodratio_M7M8 = likelihoodratio_n_M7M8)

island$likelihoodratio_M7M8 <- likelihoodratio_n_M7M8



chi_pval_n_M7M8 <- c()
for (i in 1:nrow(island)) {
  chi_pval_n_M7M8 <- c(chi_pval_n_M7M8,
                  pchisq(island$likelihoodratio_M7M8[i],
                         df=2, 
                         lower.tail = FALSE))
}

island <- island %>% 
  mutate(Chi_pvalue_M7M8 = chi_pval_n_M7M8)


island$Chi_pvalue_M7M8 <- chi_pval_n_M7M8











likelihoodratio_n_M1M2 <- c()
for (i in 1:nrow(island)) {
  likelihoodratio_n_M1M2 <- c(likelihoodratio_n_M1M2, 
                              likeRT(island$M1[i], 
                                     island$M2[i]))
}

island <- island %>%
  mutate(likelihoodratio_M1M2 = likelihoodratio_n_M1M2)

island$likelihoodratio_M1M2 <- likelihoodratio_n_M1M2



chi_pval_n_M1M2 <- c()
for (i in 1:nrow(island)) {
  chi_pval_n_M1M2 <- c(chi_pval_n_M1M2,
                       pchisq(island$likelihoodratio_M1M2[i],
                              df=2, 
                              lower.tail = FALSE))
}

island <- island %>% 
  mutate(Chi_pvalue_M1M2 = chi_pval_n_M1M2)


island$Chi_pvalue_M1M2 <- chi_pval_n_M1M2



write.csv(island, "../four/codeml/two/Book4_results.csv")
