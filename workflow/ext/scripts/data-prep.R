#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

# TODO: improve
# args = c(snakemake@input[[1]], snakemake@output[[1]])

#data <- read.csv2(snakemake@input[[1]], sep=",")
data <- read.table(args[1], sep = ",",
                   header = TRUE, row.names = 1)

# remove samples with only 1 non-zero features
data <- data[colSums(data > 0) >= 2]
# remove the all zero counts
data <- data[rowSums(data) > 0, ]

write.table(data, args[2], sep = ",")

