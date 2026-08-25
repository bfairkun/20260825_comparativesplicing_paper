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
                 "../code/MazinLeafcutterAnalysis/Expression/Chicken_UCSC.galGal6_ensv101.log2rpkm.tsv.gz ../code/config/Cordoso_Moreira_SampleList.tsv ../code/scratch/test.out.tsv.gz ../code/scratch/test.out.extended.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

RPKM_table <- args[1]
Samples_table <- args[2]
Out_table <- args[3]
Out_extended_table <- args[4]

library(data.table)
library(tidyverse)
library(qvalue)
library(broom)

samples <- read_tsv(Samples_table)

RPKM <- fread(RPKM_table) %>%
  pivot_longer(names_to = "sample", values_to = "log2rpkm", -Geneid) %>%
  inner_join(
    samples %>%
      filter(`Used library?` %in% c("Yes", "yes")) %>%
      dplyr::select(Tissue = Tissue_ForDevelopementalAnalysis, sample=ID, Ordinal_stage, Species=ID_Species) %>%
      separate_rows(Tissue, sep=",")
  ) %>%
  mutate(NonZeroPSI = log2rpkm>-5) %>%
  mutate(group = paste(Geneid, Species, Tissue, sep=';'))


options(warn=-1)

RPKM.filtered <- RPKM %>%
  # head(1E6) %>%
  group_by(Geneid, Species, Tissue) %>%
  mutate(NumberNonZeroPSI_Samples = sum(NonZeroPSI, na.rm=T)) %>%
  ungroup() %>%
  filter(NumberNonZeroPSI_Samples > 5)

dim(RPKM)
dim(RPKM.filtered)

spearman.coefs <- lapply(split(RPKM.filtered, RPKM.filtered$group), function(x){cor.test(x$log2rpkm, x$Ordinal_stage, method = "spearman", exact=F)$estimate})
spearman.pvals <- lapply(split(RPKM.filtered, RPKM.filtered$group), function(x){cor.test(x$log2rpkm, x$Ordinal_stage, method = "spearman", exact=F)$p.value})


# Bring it together
results <- data.frame(group = names(spearman.coefs), corr = unlist(spearman.coefs), P=unlist(spearman.pvals), row.names = NULL) %>%
  separate(group, into=c("Geneid", "Species", "Tissue"), sep=";") %>%
  group_by(Tissue, Species) %>%
  mutate(q = qvalue(P)$qvalue)

# results %>%
#   ggplot(aes(x=P)) +
#   geom_histogram() +
#   facet_wrap(~Tissue)

write_tsv(results, Out_table)

### calculated and write out extended output table
expression.table.tidy <- fread(RPKM_table) %>%
  pivot_longer(names_to = "ID", values_to = "log2expression", -Geneid) %>%
  inner_join(
    samples %>%
      filter(`Used library?` %in% c("yes", "Yes")) %>%
      dplyr::select(ID:ID_Stage, Ordinal_stage, Tissue_ForDevelopementalAnalysis) %>%
      separate_rows(Tissue_ForDevelopementalAnalysis, sep=",")
  )

expression.table.tidy.qqnormed <- expression.table.tidy %>%
  group_by(Geneid, Tissue_ForDevelopementalAnalysis) %>%
  mutate(log2expression.normalized = scale(log2expression)) %>%
  mutate(var.log2expression = var(log2expression, na.rm = T)) %>%
  mutate(maxExpression = max(log2expression, na.rm=T)) %>%
  ungroup() %>%
  group_by(ID, Tissue_ForDevelopementalAnalysis) %>%
  mutate(log2expression.normalized = qnorm(dense_rank(log2expression.normalized)/(sum(!is.na(log2expression.normalized)) + 1))) %>%
  ungroup()


lm.results <- expression.table.tidy.qqnormed %>%
  # head(1E6) %>%
  filter(var.log2expression > 0) %>%
  mutate(Ordinal_stage_scaled = scale(Ordinal_stage)) %>%
  group_by(Tissue_ForDevelopementalAnalysis, Geneid) %>%
  do(tidy(lm(log2expression.normalized ~ Ordinal_stage_scaled, data = .))) %>%
  filter(term == "Ordinal_stage_scaled") %>%
  dplyr::select(Tissue_ForDevelopementalAnalysis, Geneid, lm_estimate = estimate, lm_p = p.value)

spearman.results <- expression.table.tidy.qqnormed %>%
  # head(1E6) %>%
  filter(var.log2expression > 0) %>%
  group_by(Tissue_ForDevelopementalAnalysis, Geneid) %>%
  filter(sum(is.finite(log2expression.normalized)) >= 5) %>%
  summarize(
    # Spearman correlation results
    cor_test = list(cor.test(Ordinal_stage, log2expression.normalized, method = 'spearman')),
    spearman_estimate = map_dbl(cor_test, "estimate"),
    spearman_p = map_dbl(cor_test, "p.value")
  ) %>%
  dplyr::select(Tissue_ForDevelopementalAnalysis, Geneid, spearman_estimate, spearman_p)

# full_join(lm.results, spearman.results) %>%
#   ggplot(aes(x=lm_estimate, y=spearman_estimate)) +
#   geom_point(alpha=0.05) +
#   coord_cartesian(xlim=c(-2,2))

full_join(lm.results, spearman.results) %>%
  group_by(Tissue_ForDevelopementalAnalysis) %>%
  mutate(lm_q = qvalue(lm_p)$qvalue,
         spearman_q = qvalue(spearman_p)$qvalue) %>%
  ungroup() %>%
  dplyr::select(Tissue = Tissue_ForDevelopementalAnalysis, everything() ) %>%
  full_join(
    results %>%
      dplyr::rename(UnNormSpearman_estimate=corr, UnNormSpearman_p=P, UnNormSpearman_q=q) %>%
      dplyr::select(-Species)
  ) %>%
  full_join(
    expression.table.tidy.qqnormed %>%
      dplyr::select(Geneid, maxExpression, Tissue = Tissue_ForDevelopementalAnalysis) %>%
      distinct()
  ) %>%
  write_tsv(Out_extended_table)

