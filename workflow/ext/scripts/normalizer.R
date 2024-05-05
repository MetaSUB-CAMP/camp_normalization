#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

normalize <- function(sample_prepped, sample_out, method, batch_factor=NULL){
  
  ftable <- read.table(sample_prepped, sep = ",", header = TRUE, row.names = 1)
  
  # Total Sum Scaling
  if(method == "tss"){
    ftable_norm <- sweep(ftable, 2, colSums(ftable), "/")
  }
  
  # Upper Quantile
  if(method == "uq"){
    uppquant <- quantile(ftable[ftable>0])["75%"]
    ftable_norm = sweep(ftable, 2, uppquant, "/")
  }
  
  # Median 
  if(method == "med"){
    med <- median(ftable[ftable>0])
    ftable_norm = sweep(ftable, 2, med, "/")
  }
  
  # Cumulative Sum Scaling
  if(method == "css"){
    require(metagenomeSeq)
    meta_ftable <- newMRexperiment(ftable)
    meta_ftable_css <- cumNorm(meta_ftable, p=cumNormStatFast(meta_ftable))
    ftable_norm <- data.frame(MRcounts(meta_ftable_css, norm=TRUE, log=TRUE))
  }
  
  # Centered log-ratio
  if(method == "clr"){
    require("compositions")
    ftable_norm <- data.frame(clr(ftable))
  }
  
  # Additive log-ratio
  if(method == "alr"){
    require("compositions")
    ftable_norm <- data.frame(alr(ftable))
  }

  # Isometric log-ratio
  if(method == "ilr"){
    require("compositions")
    ftable_norm <- data.frame(alr(ftable))
  }
  
  # Blom
  if(method == "blom"){
    # TSS normalization
    ftable_tss <- sweep(ftable, 2, colSums(ftable), "/")
    # function for blom transformation
    blom.func <- function(tab){
      tab <- as.matrix(tab)
      # a small noise term is added before data transformation to handle the ties
      noise <- matrix(rnorm(nrow(tab)*ncol(tab),mean=0,sd=10^-10),nrow=nrow(tab),ncol=ncol(tab))
      tab_noised <- tab+noise
      c <- 3/8
      tab_trans <- as.data.frame(t(apply(tab_noised,1,function(x) qnorm((rank(x)-c)/(ncol(tab_noised)-2*c+1)))))
      as.data.frame(tab_trans)
    }
    # blom transformation
    ftable_norm <- blom.func(ftable_tss)
  }
  
  # Non-paranormal normalization
  if(method=="npn"){
    require("huge")
    ftable_tss <- sweep(ftable, 2, colSums(ftable), "/")
    # NPN transformation
    ftable_norm <- as.data.frame(t(huge.npn(t(ftable_tss),npn.func="truncation")))
  }
  
  #Trimmed mean of M-values
  if(method == "tmm"){
    require(edgeR)
    # function for TMM normalization
    tmm.func <- function(tab){
      tab <- as.matrix(tab)
      tab_dge <- DGEList(counts=tab)
      tab_tmm_dge <- calcNormFactors(tab_dge, method="TMM")
      tab_norm <- cpm(tab_tmm_dge)
      as.data.frame(tab_norm)
    }
    # TMM normalization
    ftable_norm <- tmm.func(ftable)
  }
  
  # ComBat
  if(method=="combat"){
    require(sva)
    # TSS normalization
    ftable_tss <- sweep(ftable, 2, colSums(ftable), "/")
    ftable_tss[ftable_tss==0] <- min(ftable_tss[ftable_tss!=0])*0.65
    # log transformation
    trans_ftable <- log(ftable_tss)
    # combat
    ftable_norm <- as.data.frame(ComBat(trans_ftable, batch=batch_factor))
  }

  # limma
  if(method=="limma"){
    require(limma)
    # TSS normalization
    ftable_tss <- sweep(ftable, 2, colSums(ftable), "/")
    ftable_tss[ftable_tss==0] <- min(ftable_tss[ftable_tss!=0])*0.65
    # log transformation
    trans_ftable <- log(ftable_tss)
    # combat
    ftable_norm <- as.data.frame(removeBatchEffect(trans_ftable, batch=batch_factor))
  }
  
  write.csv(ftable_norm, sample_out)
  
}

# TODO: poor logic, improve
# args = c(snakemake@input[[1]], snakemake@input[[2]], snakemake@params[[1]], snakemake@output[[1]])

sample_prepped <- args[1]
sample_out <- args[4]
method <- args[3]
pathtometa <- args[2]

batchcorr_methods <- c('conqur', 'combat', 'limma')
if(method %in%  batchcorr_methods){
  meta <- read.table(pathtometa, sep = ",", header = TRUE, row.names = 1)
  batch_factor <- factor(meta$batch)
  normalize(sample_prepped, sample_out, method, batch_factor)
} else{
  normalize(sample_prepped, sample_out, method)
}
  
