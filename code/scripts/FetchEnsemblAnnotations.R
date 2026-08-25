#!/usr/bin/env Rscript
# Gene/transcript biotypes and human orthologs for the seven study species, from Ensembl.
#
# Usage: Rscript FetchEnsemblAnnotations.R <transcript_info.tsv.gz> <human_homologs.tsv.gz>
#
# Two notes on versions:
#   * UCSC's "Ensembl" tracks are snapshots. For mouse, UCSC's mm39 Ensembl track matches
#     Ensembl v110, so ENSEMBL_VERSION below is pinned rather than left at "current" --
#     otherwise the join in FixMouseEnsemblBed.py drops transcripts.
#   * biomaRt talks to a live service, so this rule is network-dependent and its output is
#     an intermediate in scratch/ rather than a shipped artifact.

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(readr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
out_transcripts <- args[1]
out_homologs <- args[2]

ENSEMBL_VERSION <- 110

DATASETS <- c(
  Human   = "hsapiens_gene_ensembl",
  Macaque = "mmulatta_gene_ensembl",
  Mouse   = "mmusculus_gene_ensembl",
  Rat     = "rnorvegicus_gene_ensembl",
  Rabbit  = "ocuniculus_gene_ensembl",
  Opossum = "mdomestica_gene_ensembl",
  Chicken = "ggallus_gene_ensembl"
)

TRANSCRIPT_ATTRS <- c(
  "ensembl_gene_id", "ensembl_gene_id_version", "gene_biotype",
  "ensembl_transcript_id", "ensembl_transcript_id_version", "transcript_biotype",
  "external_gene_name"
)

HOMOLOG_ATTRS <- c(
  "ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_perc_id",
  "hsapiens_homolog_perc_id_r1", "hsapiens_homolog_orthology_confidence",
  "hsapiens_homolog_orthology_type"
)

get_mart <- function(dataset) {
  useEnsembl(biomart = "genes", dataset = dataset, version = ENSEMBL_VERSION)
}

message("fetching transcript info ...")
transcripts <- imap_dfr(DATASETS, function(dataset, species) {
  message("  ", species)
  getBM(attributes = TRANSCRIPT_ATTRS, mart = get_mart(dataset)) %>%
    mutate(Species = species)
}) %>%
  select(Species, all_of(TRANSCRIPT_ATTRS))
write_tsv(transcripts, out_transcripts)
message("wrote ", nrow(transcripts), " transcript rows -> ", out_transcripts)

message("fetching human orthologs ...")
homologs <- imap_dfr(DATASETS[names(DATASETS) != "Human"], function(dataset, species) {
  message("  ", species)
  getBM(attributes = HOMOLOG_ATTRS, mart = get_mart(dataset)) %>%
    mutate(Species = species)
}) %>%
  filter(hsapiens_homolog_ensembl_gene != "") %>%
  left_join(
    transcripts %>% distinct(Species, ensembl_gene_id, external_gene_name),
    by = c("Species", "ensembl_gene_id")
  ) %>%
  select(Species, all_of(HOMOLOG_ATTRS), external_gene_name)
write_tsv(homologs, out_homologs)
message("wrote ", nrow(homologs), " ortholog rows -> ", out_homologs)
