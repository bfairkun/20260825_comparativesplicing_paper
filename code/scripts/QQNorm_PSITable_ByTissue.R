#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79
# @file        :
# @created     : 2024-10-01 09:43:04
#
# @description : standardize and qqnorm PSI table by tissue
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text="../code/rna-seq/SplicingAnalysis/leafcutter/Human_UCSC.hg38_GencodeComprehensive46/juncTableBeds/PSI.sorted.bed.gz ../code/config/Cordoso_Moreira_SampleList.tsv ../code/scratch/PSI.norm. ../code/scratch/test.qqnorm.PSI.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(data.table)
library(broom)

f_in <- args[1]
sample_list_in <- args[2]
f_out_tables_prefix <- args[3]
f_out <- args[4]


PSI.table <- read_tsv(f_in)

CordosoSamples <- read_tsv(sample_list_in)

PSI.table.tidy <- PSI.table %>%
  pivot_longer(names_to = "ID", values_to = "PSI", -c(1:6)) %>%
  inner_join(
    CordosoSamples %>%
      filter(`Used library?` %in% c("yes", "Yes")) %>%
      dplyr::select(ID:ID_Stage, Ordinal_stage, Tissue_ForDevelopementalAnalysis) %>%
      separate_rows(Tissue_ForDevelopementalAnalysis, sep=",")
  ) %>%
  group_by(junc, Tissue_ForDevelopementalAnalysis) %>%
  mutate(PSI.normalized = scale(PSI)) %>%
  mutate(var.PSI = var(PSI, na.rm = T)) %>%
  ungroup() %>%
  group_by(ID, Tissue_ForDevelopementalAnalysis) %>%
  mutate(PSI.normalized = qnorm(dense_rank(PSI.normalized)/(sum(!is.na(PSI.normalized)) + 1))) %>%
  ungroup()

# write out normalized psi tables by tissue
PSI.tables.list <- PSI.table.tidy %>%
  mutate(PSI.normalized = round(PSI.normalized, 4)) %>%
  dplyr::select(1:6,  ID, PSI.normalized) %>%
  split(paste(PSI.table.tidy$Tissue_ForDevelopementalAnalysis)) %>%
  lapply(function(x) pivot_wider(x, names_from = "ID", values_from = "PSI.normalized"))


lapply(seq_along(PSI.tables.list), function(i) write_tsv(PSI.tables.list[[i]], paste0(f_out_tables_prefix, names(PSI.tables.list)[[i]],  ".tsv.gz") ))

sample_n_of <- function(data, size, ...) {
  dots <- quos(...)

  group_ids <- data %>%
    group_by(!!! dots) %>%
    group_indices()

  sampled_groups <- sample(unique(group_ids), size)

  data %>%
    filter(group_ids %in% sampled_groups)
}


# Get spearman coef
# and linear model summary stats
lm.results <- PSI.table.tidy %>%
  # head(1E6) %>%
  filter(var.PSI > 0) %>%
  mutate(Ordinal_stage_scaled = scale(Ordinal_stage)) %>%
  group_by(Tissue_ForDevelopementalAnalysis, junc) %>%
  do(tidy(lm(PSI.normalized ~ Ordinal_stage_scaled, data = .))) %>%
  filter(term == "Ordinal_stage_scaled") %>%
  dplyr::select(Tissue_ForDevelopementalAnalysis, junc, lm_estimate = estimate, lm_p = p.value)

spearman.results <- PSI.table.tidy %>%
  # sample_n_of(100, junc) %>%
  # head(1E6) %>%
  filter(var.PSI > 0) %>%
  group_by(Tissue_ForDevelopementalAnalysis, junc) %>%
  filter(sum(is.finite(PSI.normalized)) >= 5) %>%
  summarize(
    # Spearman correlation results
    cor_test = list(cor.test(Ordinal_stage, PSI.normalized, method = 'spearman')),
    spearman_estimate = map_dbl(cor_test, "estimate"),
    spearman_p = map_dbl(cor_test, "p.value")
  ) %>%
  dplyr::select(Tissue_ForDevelopementalAnalysis, junc, spearman_estimate, spearman_p)

full_join(lm.results, spearman.results) %>%
  ggplot(aes(x=lm_estimate, y=spearman_estimate)) +
  geom_point(alpha=0.2) +
  coord_cartesian(xlim=c(-2,2))

full_join(lm.results, spearman.results) %>%
  write_tsv(f_out)

