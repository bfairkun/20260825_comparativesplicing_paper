#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : SpearmanCor_Mazin_LeafcutterPSI
# @created     : Thursday Aug 22, 2024 15:32:08 CDT
#
# @description : Spearman correlation test from PSI table of Mazin leafcutter processed data
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "../code/rna-seq/SplicingAnalysis/leafcutter/Chicken_UCSC.galGal6_ensv101/juncTableBeds/PSI_ByMax.sorted.bed.gz ../code/config/Cordoso_Moreira_SampleList.tsv ../code/scratch/test.out.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

PSI_table <- args[1]
Samples_table <- args[2]
Out_table <- args[3]

library(data.table)
library(tidyverse)
library(qvalue)

samples <- read_tsv(Samples_table)

leafcutter.PSI <- fread(PSI_table) %>%
  dplyr::select(-c(1:3,5,6)) %>%
  pivot_longer(names_to = "sample", values_to = "PSI", -junc) %>%
  inner_join(
    samples %>%
      filter(`Used library?` %in% c("Yes", "yes")) %>%
      dplyr::select(Tissue = Tissue_ForDevelopementalAnalysis, sample=ID, Ordinal_stage, Species=ID_Species) %>%
      separate_rows(Tissue, sep=",")
  ) %>%
  mutate(NonZeroPSI = PSI>0) %>%
  mutate(group = paste(junc, Species, Tissue, sep=';'))


options(warn=-1)

leafcutter.PSI.filtered <- leafcutter.PSI %>%
  # head(1E6) %>%
  group_by(junc, Species, Tissue) %>%
  mutate(NumberNonZeroPSI_Samples = sum(NonZeroPSI, na.rm=T)) %>%
  ungroup() %>%
  filter(NumberNonZeroPSI_Samples > 5)

spearman.coefs <- lapply(split(leafcutter.PSI.filtered, leafcutter.PSI.filtered$group), function(x){cor.test(x$PSI, x$Ordinal_stage, method = "spearman", exact=F)$estimate})
spearman.pvals <- lapply(split(leafcutter.PSI.filtered, leafcutter.PSI.filtered$group), function(x){cor.test(x$PSI, x$Ordinal_stage, method = "spearman", exact=F)$p.value})


# Bring it together
results <- data.frame(group = names(spearman.coefs), corr = unlist(spearman.coefs), P=unlist(spearman.pvals), row.names = NULL) %>%
  separate(group, into=c("junc", "Species", "Tissue"), sep=";") %>%
  separate(junc, into=c("start", "stop", "chrom", "cluster"), sep=":", remove = F, convert=T) %>%
  group_by(Tissue, Species) %>%
  mutate(q = qvalue(P)$qvalue)

write_tsv(results, Out_table)
