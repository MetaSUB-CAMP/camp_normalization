#data <- read.csv2(snakemake@input[[1]], sep=",")
data <- read.table(snakemake@input[[1]], sep = ",",
                   header = TRUE, row.names = 1)

# remove the all zero counts
data <- data[rowSums(data) > 0, ]
# remove samples with only 1 non-zero features
data <- data[colSums(data > 0) >= 2]

write.table(data, snakemake@output[[1]], sep = ",")

