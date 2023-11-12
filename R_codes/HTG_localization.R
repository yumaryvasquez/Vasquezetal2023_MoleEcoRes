#localization 

df  <- read.csv("../all.txt",
                        sep="\t",skip=4, header=F)
df2 <- df %>% group_by(V1, V2) %>% summarise(max = max(V12, na.rm=TRUE))
df3<-left_join(df, df2, by=c('V1'='V1', 'V2'='V2', 'V12'='max'))
df3<-left_join(df2, df, by=c('V1'='V1', 'V2'='V2', 'max'='V12'))



#ideogram
require(RIdeogram)

chromosome <-read.table("../gene_ideogram/chromosome_length.csv", sep = "\t", header = T, stringsAsFactors = F)
gene_density2 <- read.table("../gene_ideogram/all2.csv", sep = "\t", header = T, stringsAsFactors = F)
HTG_localization <- read.table("../gene_ideogram/htg_places2.csv", sep = "\t", header = T, stringsAsFactors = F)
#chromosome$CE_start <- 0
#chromosome$CE_end <-chromosome$End
HTG_localization$Chr <- as.character(HTG_localization$Chr)
#ideogram(karyotype = chromosome)
#convertSVG("chromosome.svg", device = "png") 
#ideogram(karyotype = chromosome, overlaid = gene_density2)
#convertSVG("chromosome.svg", device = "png") 
ideogram(karyotype = chromosome, overlaid = gene_density2, 
         label = HTG_localization, label_type = "marker",
         colorset1 = c("#DAD7CD", "#588157", "#344E41"))
convertSVG("chromosome.svg", device = "png") 
  





chromosome <-read.table("../gene_ideogram/chromosome_length.csv", sep = "\t", header = T, stringsAsFactors = F)
gene_density2 <- read.table("../gene_ideogram/all2.csv", sep = "\t", header = T, stringsAsFactors = F)
HTG_localization <- read.table("../gene_ideogram/mito_support2.csv", sep = "\t", header = T, stringsAsFactors = F)
#chromosome$CE_start <- 0
#chromosome$CE_end <-chromosome$End
HTG_localization$Chr <- as.character(HTG_localization$Chr)
#ideogram(karyotype = chromosome)
#convertSVG("chromosome.svg", device = "png") 
#ideogram(karyotype = chromosome, overlaid = gene_density2)
#convertSVG("chromosome.svg", device = "png") 
ideogram(karyotype = chromosome, overlaid = gene_density2, 
         label = HTG_localization, label_type = "marker",
         colorset1 = c("#DAD7CD", "#588157", "#344E41"))
convertSVG("chromosome.svg", device = "png") 




















  
  #Data processing:
df  <- read.csv("../gene_ideogram/unique_chromosome9.txt",
                  sep="\t",skip=2, header=F)
df2 <- df %>% arrange(V1)  
df2 <- df2 %>% 
  rename(
    start = V1,
    end = V2
  )
df <- read.csv("../gene_ideogram/chr9_local.csv", header = T)
df1 <- df[,-1]
df1 <- df1 %>% 
  rename(
    start = Start,
    end = End
  )
df1 <- na.omit(df1)


# Function to check if a value falls within a range
in_range <- function(value, range) {
  return(value >= range[1] && value <= range[2])
}

# Initialize count column in df1
df1$count <- 0

# Loop through each row of df1
for (i in 1:nrow(df1)) {
  # Initialize count for this row
  count <- 0
  # Loop through each row of df2
  for (j in 1:nrow(df2)) {
    # Check if the range in df2 overlaps with the range in df1
    if (in_range(df2$start[j], c(df1$start[i], df1$end[i])) || 
        in_range(df2$end[j], c(df1$start[i], df1$end[i])) ||
        (df2$start[j] < df1$start[i] && df2$end[j] > df1$end[i])) {
      # Increment count
      count <- count + 1
    }
  }
  # Update count column in df1
  df1$count[i] <- count
}

write.csv(df1,"../gene_ideogram/chr9_local.csv" )


#next step

df1 <- read.csv("../gene_ideogram/chr1_local.csv", 
               header = T,  row.names = 1)
df2 <- read.csv("../gene_ideogram/chr2_local.csv", 
               header = T,  row.names = 1)
df3 <- read.csv("../gene_ideogram/chr3_local.csv", 
                header = T,  row.names = 1)
df4 <- read.csv("../gene_ideogram/chr4_local.csv", 
                header = T,  row.names = 1)
df5 <- read.csv("../gene_ideogram/chr5_local.csv", 
                header = T,  row.names = 1)
df6 <- read.csv("../gene_ideogram/chr6_local.csv", 
                header = T,  row.names = 1)
df7 <- read.csv("../gene_ideogram/chr7_local.csv", 
                header = T,  row.names = 1)
df8 <- read.csv("../gene_ideogram/chr8_local.csv", 
                header = T,  row.names = 1)
df9 <- read.csv("../gene_ideogram/chr9_local.csv", 
                header = T,  row.names = 1)
df1 <- cbind(Chr = 1, df1)
df2 <- cbind(Chr = 2, df2)
df3 <- cbind(Chr = 3, df3)
df4 <- cbind(Chr = 4, df4)
df5 <- cbind(Chr = 5, df5)
df6 <- cbind(Chr = 6, df6)
df7 <- cbind(Chr = 7, df7)
df8 <- cbind(Chr = 8, df8)
df9 <- cbind(Chr = 9, df9)

df_all <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9)


write.csv(df_all, "../gene_ideogram/all.csv")










#DO NOT USE!!!

for (i in 1:nrow(ranges)){
  
  r1 <- ranges[i,2]
  r2 <- range[i, 3]
  
}


df1 <- tibble(start = c(1, 10, 20), end = c(5, 15, 25))
df2 <- tibble(min = c(2, 7, 12, 17), max = c(8, 13, 18, 23))

# Create a function to count the number of overlaps
count_overlaps <- function(start, end, values) {
  sum(between(values, start, end))
}

# Use a nested loop to iterate over each row in df1 and df2 and count the overlaps
df1 %>%
  mutate(count = map_dbl(seq_along(start), function(i) {
    values <- df2 %>% 
      filter(max >= start[i], min <= end[i]) %>% 
      pull(min:max) 
    count_overlaps(start[i], end[i], values)
  }))


# Generate sample data
df1 <- data.frame(start = c(1, 6, 11,21), end = c(5, 10, 20,30))
df2 <- data.frame(start = c(1, 2,1, 12, 22), end = c(2, 6,25, 16, 26))

