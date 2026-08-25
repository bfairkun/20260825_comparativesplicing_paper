#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PrepAllJuncsFor_JunctionClassifier
# @created     : Thursday Aug 29, 2024 15:14:13 CDT
#
# @description :
# rule ConvertLeacfutterJuncNames_To_True_Bed_Coords:
#     input:
#         junclist = "rna-seq/SplicingAnalysis/leafcutter/{GenomeName}/clustering/leafcutter_perind.counts.gz",
#     output:
#         "MazinLeafcutterAnalysis/ClassifyJuncs/{GenomeName}.juncList.tsv.gz"
#     shell:
#         """
#         zcat {input} | awk 'NR==1 {{ print $1 }} NR>1 {{split($1, a, ":"); print a[1]":"a[2]":"a[3]-1":"a[4] }}' | gzip - > {output}
#         """
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "rna-seq/SplicingAnalysis/leafcutter/Human_UCSC.hg38_GencodeComprehensive46/clustering/leafcutter_perind.counts.gz rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/Human_UCSC.hg38_GencodeComprehensive46.uniq.annotated.tsv.gz MazinLeafcutterAnalysis/ClassifyJuncs//Human_UCSC.hg38_GencodeComprehensive46.Clustered.juncList.tsv.gz MazinLeafcutterAnalysis/ClassifyJuncs//Human_UCSC.hg38_GencodeComprehensive46.AllObserved.juncList.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

fn_in.leaf <- args[1]
fn_in.AllObs <- args[2]
fn_out <- args[3]
fn_out.AllObs <- args[4]

library(tidyverse)
library(data.table)


dat.in <- fread(fn_in.leaf, select = 1, sep=' ') %>%
  dplyr::select(leafName = chrom) %>%
  separate(leafName, into=c("chrom", "start", "end", "cluName"), sep=":", convert=T, remove=F) %>%
  mutate(end = as.integer(end - 1), start=as.integer(start)) %>%
  mutate(leafName = str_glue("{chrom}:{start}:{end}:{cluName}"))

## Check that dat.in is proper bed coords
# dat.in %>%
#   # head(10000) %>%
#   mutate(strand = str_extract(cluName, "[+-]$"), score=".") %>%
#   arrange(chrom, start, end) %>%
#   head(10000) %>%
#   dplyr::select(chrom, start, end, leafName, score, strand) %>%
#   write_tsv("scratch/leaf.bed", col_names = F)

all.juncs <- read_tsv(fn_in.AllObs) %>%
  mutate(end = end - 1)

all.juncs %>%
  dplyr::select(1:4) %>%
  head()

all.juncs %>%
  left_join(dat.in) %>%
  mutate(IsClustered = !is.na(cluName)) %>%
  count(IsClustered, strand)

dat.in %>% dplyr::select(cluName)

dat.in %>%
    dplyr::select(chrom = leafName) %>%
    mutate(chrom = as.character(chrom)) %>%
    write_tsv(fn_out)

all.juncs %>%
  left_join(dat.in) %>%
  mutate(start = as.character(as.integer(start)), end=as.character(as.integer(end))) %>%
  mutate(leafName = case_when(
    is.na(cluName) ~ str_glue("{chrom}:{start}:{end}:clu_NA_{strand}"),
    TRUE ~ leafName
  )) %>%
  arrange(chrom, start, end) %>%
  dplyr::select(chrom = leafName) %>%
  write_tsv(fn_out.AllObs)
