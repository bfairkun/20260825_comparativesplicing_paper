#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79
# @file        :
# @created     : 2024-10-14 18:59:00
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text="MazinLeafcutterAnalysis/ContrastGroupFiles/BeforeB_vs_AferS.Brain.Human.txt MazinLeafcutterAnalysis/ClassifyJuncs/Chicken_UCSC.galGal6_ensv101.AllObserved._junction_classifications.txt MazinLeafcutterAnalysis/Contrasts_ds/BeforeB_vs_AferS.Brain.Human_cluster_significance.txt MazinLeafcutterAnalysis/Contrasts_ds/BeforeB_vs_AferS.Brain.Human_effect_sizes.txt MazinLeafcutterAnalysis/Expression/Human_UCSC.hg38_GencodeComprehensive46.counts.tsv.gz MazinLeafcutterAnalysis/Contrasts_ds_tidy/BeforeB_vs_AferS.Brain.Human.de.tsv.gz MazinLeafcutterAnalysis/Contrasts_ds_tidy/BeforeB_vs_AferS.Brain.Human.joined.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}



fn_in_groups <- args[1]
fn_in_juncclassifications <- args[2]
fn_in_clustersig <- args[3]
fn_in_leafeffect <- args[4]
fn_in_counts <- args[5]
fn_out_de <- args[6]
fn_out_tidy <- args[7]

library(tidyverse)
library(data.table)
library(edgeR)

contrast.groups <- read_delim(fn_in_groups, col_names = c("ID", "Group"), delim = ' ')

Counts <- fread(fn_in_counts) %>%
  dplyr::select(Geneid, contrast.groups$ID) %>%
  column_to_rownames("Geneid")

y <- DGEList(counts=Counts, group=contrast.groups$Group)

y$samples

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

dim(y)

results <- y %>%
  calcNormFactors() %>%
  estimateDisp() %>%
  exactTest() %>%
  topTags(n=Inf) %>%
  as.data.frame() %>%
  rownames_to_column("Geneid")



cluster.sig <- read_tsv(fn_in_clustersig)
effect.sizes <- read_tsv(fn_in_leafeffect)
junc.classifications <- read_tsv(fn_in_juncclassifications)


results %>%
  write_tsv(fn_out_de)

effect.sizes %>%
  separate(intron, into=c("chrom", "start", "stop", "cluster"), sep=":", convert=T) %>%
  mutate(Intron_coord = str_glue("{chrom}:{start}-{stop-1}")) %>%
  unite(cluster, chrom, cluster, sep = ":") %>%
  inner_join(cluster.sig) %>%
  left_join(junc.classifications) %>%
  left_join(
    results,
    by=c("Gene_name"="Geneid")
  ) %>%
  write_tsv(fn_out_tidy)




