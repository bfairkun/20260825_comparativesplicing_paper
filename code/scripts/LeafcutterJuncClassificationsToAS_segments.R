#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : test
# @created     : Wednesday Jul 03, 2024 10:00:17 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
  args <- scan(text=
                 "../code/kaessmanAnalysis/JuncListsForLeafcutter2/Chicken_ensemblv84._junction_classifications.txt ../code/kaessmanAnalysis/NMDFinderB_ASSegments/Chicken_ensemblv84.tsv.gz ../code/GenomeFiles/Chicken_ensemblv84/Reannotated.B.bed.gz ../code/scratch/OutFull.tsv.gz ../code/scratch/OutCollapsed.tsv.gz", what='character')
} else{
  args <- commandArgs(trailingOnly=TRUE)
}

Input_leaf2JunctionClassifications_fn <- args[1]
Input_NMDFinderBTsv_fn <- args[2]
Input_TranscriptsBed12_fn <- args[3]
Output_FullJoinedTable_fn <- args[4]
Output_CollapsedTable_fn <- args[5]

library(tidyverse)

JunctionsAnnotations <- read_tsv(Input_leaf2JunctionClassifications_fn)

AS_segmentAnnotations <- read_tsv(Input_NMDFinderBTsv_fn)

transcripts_to_genes <- read_tsv(Input_TranscriptsBed12_fn, col_names=F) %>%
  dplyr::select(4,13) %>%
  dplyr::rename("transcript"="X4", "Gene_name"="X13") %>%
  distinct()

AS_segmentAnnotations.tidy <- AS_segmentAnnotations %>%
  filter(!is.na(WhichIsoformIsAnnotated)) %>%
  distinct() %>%
  separate_rows(LongIsoform_UniqueJuncs, ShortIsoform_UniqueJuncs, sep=",") %>%
  pivot_longer(cols=LongIsoform_UniqueJuncs:ShortIsoform_NMDFinderB, names_pattern="^(.+?)_(.+)$", names_to = c("isoform", ".value")) %>%
  distinct() %>%
  mutate(UniqueJuncs = str_replace(UniqueJuncs, "^(.+)_([0-9]*)_([0-9]*)_([-+])$", "\\1;\\2;\\3;\\4"))


JoinedAnnotations <- AS_segmentAnnotations.tidy %>%
  separate(UniqueJuncs, into=c("junc_chrom", "junc_start", "junc_end", "junc_strand"), sep=";", convert=T, remove=F) %>%
  inner_join(transcripts_to_genes) %>%
  mutate(Intron_coord = str_glue("{junc_chrom}:{junc_start}-{junc_end+1}")) %>%
  inner_join(JunctionsAnnotations)

AS_segment_AnnotationsFromJuncs <- JoinedAnnotations %>%
  group_by(isoform, AS_segment, Gene_name) %>%
  summarise(PercentCoding = sum(Coding)/n(), n=n()) %>%
  ungroup() %>%
  mutate(isoform = recode(isoform, "LongIsoform"="AS_seg_included", "ShortIsoform"="AS_seg_excluded"))

write_tsv(JoinedAnnotations, Output_FullJoinedTable_fn)
write_tsv(AS_segment_AnnotationsFromJuncs, Output_CollapsedTable_fn)

