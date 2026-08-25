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
                 "../code/kaessman_AS_dat/Supplementary_Data/Supplementary_Data_9.csv ../code/kaessman_AS_dat/AS_segment_lists/", what='character')
} else{
  args <- commandArgs(trailingOnly=TRUE)
}

Input_fn <- args[1]
OutputPrefix <- args[2]

library(tidyverse)
library(data.table)

AllSegmentEffects <- fread(Input_fn, nrows = Inf, select=1:8) %>%
  rename("segment.name"="V1") %>%
  mutate(species = str_replace(segment.name, "^(.+?)\\..+$", "\\1")) %>%
  mutate(species = recode(species, "hum"="Human_ensemblv75", "chi"="Chicken_ensemblv84", "mac"="Macaque_ensemblv84","mou"="Mouse_ensemblv84", "opo"="Opossum_ensemblv84", "rab"="Rabbit_ensemblv_84", "rat"="Rat_ensemblv_79")) %>%
  mutate(color = recode(as.type, "AA"="#7570b3","AD"="#e7298a","CE"="#1b9e77","RI"="#d95f02","complex"="#636363")) %>%
  mutate(strand = recode(strand, `-1`="-", `1`="+"))


AllSegmentEffects %>% head()

AllSegmentEffects %>%
  mutate(start = start - 1) %>%
  mutate(score = ".", thickStart = start, thickStop=stop) %>%
  dplyr::select(species, chr=chr_id, start, stop,segment.name, score, strand, thickStart, thickStop, color) %>%
  arrange(species, chr, start, stop) %>%
  group_by(species) %>%
  group_walk(~ write_tsv(.x, paste0(OutputPrefix, .y$species, ".bed"), col_names = F))




