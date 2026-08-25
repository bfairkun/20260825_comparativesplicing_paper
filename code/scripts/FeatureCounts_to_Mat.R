#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : FeatureCounts_to_Mat
# @created     : Friday Sep 06, 2024 11:03:38 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "../code/featureCounts/Rabbit_UCSC.oryCun2_ensv101/AllSamplesUnstrandedCounting.Counts.txt ../code/scratch/test.mat.tsv.gz ../code/scratch/test.mat.rpkm.tsv.gz ../code/scratch/test.mat.tpm.tsv.gz ../code/scratch/test.mat.cpm.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(data.table)
library(edgeR)

f_in <- args[1]
f_out_RawCounts <- args[2]
f_out_log2RPKM <- args[3]
f_out_TPM <- args[4]
f_out_CPM <- args[5]

dat <- fread(f_in) %>%
  rename_at(-c(1:6), ~str_replace(.x, "^.+/(.+?)/Aligned.sortedByCoord.out.bam$", "\\1"))

counts <- dat %>%
  dplyr::select(-c(2:6)) %>%
  column_to_rownames("Geneid") %>%
  DGEList()

counts$counts %>%
  as.data.frame() %>%
  rownames_to_column("Geneid") %>%
  dplyr::select(Geneid, everything()) %>%
  write_tsv(f_out_RawCounts)


rpkm.mat <- calcNormFactors(counts, method='TMMwsp') %>%
  rpkm(log=T, gene.length=dat$Length, prior.count=0.1)

rpkm.mat %>%
  round(3) %>%
  as.data.frame() %>%
  rownames_to_column("Geneid") %>%
  dplyr::select(Geneid, everything()) %>%
  write_tsv(f_out_log2RPKM)

cpm.mat <- calcNormFactors(counts, method='TMMwsp') %>%
  cpm(log=T, prior.count=0)

cpm.mat %>%
  round(3) %>%
  as.data.frame() %>%
  rownames_to_column("Geneid") %>%
  dplyr::select(Geneid, everything()) %>%
  write_tsv(f_out_CPM)

tpm.mat <- log2(t( t(2**rpkm.mat) / colSums(2**rpkm.mat) ) * 1e6)


tpm.mat %>%
  round(3) %>%
  as.data.frame() %>%
  rownames_to_column("Geneid") %>%
  dplyr::select(Geneid, everything()) %>%
  write_tsv(f_out_TPM)
