import filecmp
import re
import os
import pathlib
import sys

import pandas as pd

# ---------------------------------------------------------------- configuration
# Everything read here is shipped with the project: sample metadata, genome metadata and
# contrast definitions ARE the workflow's configuration and are not derived from anything.

samples = pd.read_csv(config["samples"], sep="\t", index_col=0)
STAR_genomes = pd.read_csv(config["STAR_genomes"], sep="\t", index_col=0)

MazinGenomes = [g for g in STAR_genomes.index if g not in ["Lamprey_ensemblv_112"]]
SampleGenomes = samples["STARGenomeName"].unique()
CordosMoreiraGenomes = list(SampleGenomes)
CordosMoreiraGenomes_dict = dict(zip([g.split("_")[0] for g in CordosMoreiraGenomes],
                                     CordosMoreiraGenomes))

CordosoMoreira_df = pd.read_csv("config/Cordoso_Moreira_SampleList.tsv", sep="\t",
                                index_col=0)
CordosoMoreira_df["JuncFileGroupForParralelization"] = (
    CordosoMoreira_df.groupby("ID_Species").cumcount() % 10 + 1)

Extra_Gtfs = pd.read_csv("config/CordosoGenomes_Extra_Gtfs.tsv", sep="\t", index_col=0)
GenomesWithRefSeq = Extra_Gtfs[Extra_Gtfs["RefSeqGtf"].notna()].index
GenomesWithEnsembl = Extra_Gtfs[Extra_Gtfs["EnsemblGtf"].notna()].index

ChainFiles = pd.read_csv("config/ChainFiles.tsv", sep="\t", index_col=0)
Cordoso_contrasts = pd.read_csv("config/CordosoTimeSeriesContrasts.tsv", sep="\t",
                                index_col=0)
Cordoso_TissueContrasts = pd.read_csv("config/CordosoTissueContrasts.tsv", sep="\t",
                                      index_col=0)

# The supplementary contrast table is the authority for which contrasts the paper used and
# which libraries are in each arm. The leafcutter group files are derived from it, and the
# differential-splicing targets are restricted to it -- `CordosoTimeSeriesContrasts.tsv` also
# defines four alternative stage splits (BeforeM_vs_AfterS, BeforeB_vs_AfterB,
# BeforeB_vs_AferS, FirstHalf_vs_SecondHalf) that were exploratory and are not in the paper.
CONTRAST_TABLE = "../data/20260820_SupplementalTable_Contrasts.tsv"
_ct = pd.read_csv(CONTRAST_TABLE, sep="\t")


def _contrast_stem(name, species):
    """Contrast name + species -> the leafcutter group/result file stem."""
    m = re.match(r"^Adult (\w+) v Adult (\w+)$", name)
    if m:
        return "adult", f"{species}_{m.group(1)}-{m.group(2)}"
    m = re.match(r"^(?:Neonate|Juvenile) (\w+) v Adult \w+$", name)
    if m:
        return "dev", f"Neonate_vs_Adult.{m.group(1)}.{species}"
    m = re.match(r"^Early embryo (\w+) v Adult \w+$", name)
    if m:
        return "dev", f"First10_vs_Second10.{m.group(1)}.{species}"
    raise ValueError(f"unrecognised contrast name: {name!r}")


_stems = [_contrast_stem(r.ContrastName, r.Species) for _, r in _ct.iterrows()]
PAPER_DEV_CONTRASTS = sorted({s for kind, s in _stems if kind == "dev"})
PAPER_ADULT_CONTRASTS = sorted({s for kind, s in _stems if kind == "adult"})

maf_filelist = pd.read_csv("../data/multiz100way_maf_md5sum.txt", sep=r"\s+",
                           names=["md5sum", "fn"])
maf_filelist["fn_base"] = maf_filelist["fn"].str.replace(r"\.maf\.gz$", "", regex=True)
maf_filelist["link"] = maf_filelist["fn"].apply(
    lambda fn: f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/{fn}")

# Written by rule Export_NoEmbryoSampleFlags and shipped, so a fresh clone parses; the
# guard keeps the workflow parseable if someone deletes it before regenerating.
_noembryo = "../data/Cordoso_SampleList.WhichAreNonHumanEmbryo.tsv.gz"
Cordoso_SamplesNoHumanEmbryo = (pd.read_csv(_noembryo, sep="\t", index_col=0)
                                if os.path.exists(_noembryo) else None)

# ---------------------------------------------------------------- derived paths
# Anything generated at parse time goes under scratch/, never config/, so that a
# read-only clone can still be parsed and `git status` stays clean.
DERIVED = "scratch/derived_config"
os.makedirs(DERIVED, exist_ok=True)

# Reference genomes live inside the bundle by default, but the path must be ABSOLUTE: rules
# imported from the rna_seq module have their relative paths rewritten under `prefix:`, so a
# relative prefix would send the module looking under rna-seq/resources/. Resolving it here
# keeps the config portable while staying exempt from that rewrite. Point `GenomesPrefix` at
# a shared directory if you already have the indexes built.
GENOMES_PREFIX = os.path.abspath(
    config.get("GenomesPrefix", "resources/ReferenceGenomes")) + "/"
config["GenomesPrefix"] = GENOMES_PREFIX
config.setdefault("rna_seq", {})
config["rna_seq"]["GenomesPrefix"] = GENOMES_PREFIX

# FASTQs are downloaded into this directory. samples.tsv stores R1 relative to the
# workdir so the config is portable, but rules imported from the rna_seq module get their
# relative inputs rewritten under `prefix:`, so hand that module absolute paths.
FASTQ_DIR = os.path.abspath(config.get("fastq_dir", "CordosoMoreira_Fastq"))
_abs_samples = os.path.join(DERIVED, "samples.abs.tsv")
_s = samples.copy()
_s["R1"] = [os.path.join(FASTQ_DIR, f"{i}.fastq.gz") for i in _s.index]
_s.to_csv(_abs_samples, sep="\t")
config.setdefault("rna_seq", {})
config["rna_seq"]["samples"] = _abs_samples

# aspera batch manifest for the bulk download helper
CordosoMoreira_df["formatted_column"] = CordosoMoreira_df.apply(
    lambda row: f"{str(row['fastq_aspera']).replace('fasp.sra.ebi.ac.uk:', '')} "
                f"{FASTQ_DIR}/{row.name}.fastq.gz", axis=1)
CordosoMoreira_df["formatted_column"].to_csv(os.path.join(DERIVED, "batch.txt"),
                                             index=False, header=False)


# ---------------------------------------------------------------- helpers

def has_differences(dcmp):
    """https://stackoverflow.com/questions/4187564"""
    try:
        differences = dcmp.left_only + dcmp.right_only + dcmp.diff_files
        if differences:
            return True
        return any([has_differences(subdcmp) for subdcmp in dcmp.subdirs.values()])
    except NotADirectoryError:
        return True


def CreateSymlinksOfDir1ContentsIntoDir2(Dir1, Dir2):
    """Symlink a module's scripts into the main workdir.

    Rules imported from a module workflow may invoke scripts by a path relative to that
    module, e.g. `Rscript scripts/Foo.R`. Snakemake runs them from the main workdir, so the
    scripts have to be reachable from here too.
    """
    Dir1_sanitized = Dir1.rstrip("/") + "/"
    Dir2_sanitized = Dir2.rstrip("/") + "/"
    for filepath in pathlib.Path(Dir1_sanitized).glob("[!.]*"):
        module_script_file = os.path.abspath(filepath)
        new_script_link = Dir2_sanitized + os.path.basename(filepath)
        try:
            os.symlink(module_script_file, new_script_link)
            print(f"Making link: {new_script_link}->{module_script_file}", file=sys.stderr)
        except FileExistsError:
            if os.path.realpath(new_script_link) != module_script_file:
                print(f"not making link, fix clashing file names: {new_script_link}",
                      file=sys.stderr)


TRUTHY = ["T", "True", "TRUE", True]


def _big_chroms(genome):
    """Genomes with a chromosome over 512 Mbp need a CSI index rather than TBI."""
    return STAR_genomes.loc[genome]["ChromLargerThan512Mbp"] in TRUTHY


def GetIndexingParams(wildcards):
    return "--csi" if _big_chroms(wildcards.GenomeName) else ""


def GetIndexSuffix(wildcards):
    return "csi" if _big_chroms(wildcards.GenomeName) else "tbi"


def _expand_samples_by_flag(FormattedString, flag):
    """Expand a template over the libraries of one genome that carry a given usage flag."""
    def InputFunctionToReturn(wildcards):
        df = Cordoso_SamplesNoHumanEmbryo
        keep = df.loc[(df["GenomeName"] == wildcards.GenomeName)
                      & (df[flag].isin(TRUTHY))].index.unique()
        return expand(FormattedString, sample=keep)
    return InputFunctionToReturn


def ExpandAllSamplesInFormatStringFromGenomeNameWildcard_NoHumanEmbryo(FormattedString):
    return _expand_samples_by_flag(FormattedString, "NoHumanEmbryo")


def ExpandAllSamplesInFormatStringFromGenomeNameWildcard_NoEmbryo(FormattedString):
    return _expand_samples_by_flag(FormattedString, "NoEmbryo")
