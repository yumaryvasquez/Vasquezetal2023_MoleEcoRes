library(tidyverse)



likeRT <- function(x, y) {
  -2 * (x-y) 
}





island <- read.csv("../redo_noiso/codeml/three/m1_m2_m7_m8.csv",)

likelihoodratio_m1m2 <- c()
for (i in 1:nrow(island)) {
  likelihoodratio_m1m2 <- c(likelihoodratio_m1m2, 
                         likeRT(island$M1[i], 
                                island$M2[i]))
}

island <- island %>%
  mutate(likelihoodratio_m1m2 = likelihoodratio_m1m2)

island$likelihoodratio_m1m2 <- likelihoodratio_m1m2



chi_pval_m1m2<- c()
for (i in 1:nrow(island)) {
  chi_pval_m1m2 <- c(chi_pval_m1m2,
                  pchisq(island$likelihoodratio_m1m2[i],
                         df=2, 
                         lower.tail = FALSE))
}

island <- island %>% 
  mutate(Chi_pvalue_m1m2 = chi_pval_m1m2)


island$Chi_pvalue_m1m2 <- chi_pval_m1m2

############

likelihoodratio_m7m8 <- c()
for (i in 1:nrow(island)) {
  likelihoodratio_m7m8 <- c(likelihoodratio_m7m8, 
                            likeRT(island$M7[i], 
                                   island$M8[i]))
}

island <- island %>%
  mutate(likelihoodratio_m7m8 = likelihoodratio_m7m8)

island$likelihoodratio_m7m8 <- likelihoodratio_m7m8



chi_pval_m7m8<- c()
for (i in 1:nrow(island)) {
  chi_pval_m7m8 <- c(chi_pval_m7m8,
                     pchisq(island$likelihoodratio_m7m8[i],
                            df=2, 
                            lower.tail = FALSE))
}

island <- island %>% 
  mutate(Chi_pvalue_m7m8 = chi_pval_m7m8)


island$Chi_pvalue_m7m8 <- chi_pval_m7m8


################

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
























write.csv(island, "../redo_noiso/codeml/three/redo_m1m2m7m8_chisquared.csv")










df <- read.csv("../codeml_results/eggnog_results.csv")
df2 <- read.csv("../codeml_results/list.txt", header = F)

df3 <- left_join(df, df2, by = c("query"="V2"))
write.csv(df3, "../codeml_results/egg_nog_ortho.csv")
