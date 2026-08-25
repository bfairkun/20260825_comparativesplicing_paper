#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79
# @file        :
# @created     : 2024-10-23 14:28:02
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text="../code/MazinLeafcutterAnalysis/ClassifyJuncs/Rat_UCSC.rn7_RefSeqv108.AllObserved._junction_classifications.txt scratch/test.out.tsv.gz ../code/rna-seq/SplicingAnalysis/juncfiles/Rat_Testis_18_1824.junccounts.tsv.gz rna-seq/SplicingAnalysis/juncfiles/Rat_Testis_0dpb_1750.junccounts.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(data.table)

Annotations_fn_in <- args[1]
output_fn <- args[2]
input_fns <- args[-c(1:2)]

Annotations <- fread(Annotations_fn_in) %>%
  separate(Intron_coord, into=c("chrom", "start", "stop"), sep="[:-]", convert=T, remove=F) %>%
  mutate(junc = str_glue("{chrom}_{start}_{stop+1}_{Strand}"))

Counts <- input_fns %>%
  setNames(nm=.) %>%
  lapply(fread) %>%
  bind_rows(.id="fn") %>%
  dplyr::rename("junc"="V1", "Count"="V2")

# Counts %>% head(10)
#
# Annotations %>%
#   filter(chrom == "chr1" & start == 471892)
# "chr1_155194_326511_-", "stop+1"
# "chr1_471892_512060_+", "stop+1"

Counts %>%
  inner_join(Annotations) %>%
  group_by(fn, Annot, Coding, UTR, GencodePC) %>%
  summarise(TotalCounts = sum(Count)) %>%
  ungroup() %>%
  write_tsv(output_fn, col_names = F)
