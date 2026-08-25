#!/usr/bin/env python3
"""Reduce pipeline intermediates to the small plot-ready tables the figure notebook reads.

Each subcommand corresponds to one Snakemake rule in rules/PlotData.smk and writes exactly
one file under output/. The reductions collapse ~9 GB of intermediates into <1 MB, so the
figures can be regenerated without re-running the pipeline.

Every subcommand takes its inputs explicitly on the command line rather than globbing, so
Snakemake owns the dependency graph.
"""
import argparse
import gc
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

SPECIES = ["Human", "Macaque", "Mouse", "Rat", "Rabbit", "Opossum", "Chicken"]
SRC_ORDER = ["Gencode", "Ensembl", "RefSeq"]
RPM_THRESH = 0.02
DS_COLS = ["cluster", "deltapsi", "Intron_coord", "p.adjust", "Gene_name",
           "Annot", "Coding", "UTR", "GencodePC", "logFC", "FDR"]
GZ = dict(sep="\t", index=False, compression={"method": "gzip", "mtime": 0})


# ---------------------------------------------------------------- shared helpers

def junc_name(intron_coord, strand):
    """leafcutter2 Intron_coord + Strand -> the regtools junction name."""
    parts = intron_coord.str.rsplit(":", n=1, expand=True)
    rng = parts[1].str.rsplit("-", n=1, expand=True)
    return parts[0] + "_" + rng[0] + "_" + rng[1] + "_" + strand


def productivity_label(utr, coding, annot, gencode_pc):
    """leafcutter2 flags -> NE / UP / PR.

    Mirrors the R case_when: NA flags fall through to "PR".
    """
    flag = 1 * utr + 2 * coding + 4 * annot + 8 * gencode_pc
    return pd.Series(
        np.select([flag.isin([1, 5]), flag.isin([0, 4])], ["NE", "UP"], default="PR"),
        index=flag.index), flag


def _as_float_flags(d, cols=("Annot", "Coding", "UTR", "GencodePC")):
    return {c: d[c].map({True: 1.0, False: 0.0, "True": 1.0, "False": 0.0}) for c in cols}


def species_of(genome):
    return genome.split("_")[0]


def one_event_per_cluster(d):
    """Filter one differential-splicing contrast to one representative event per cluster.

    Keeps the unproductive junction with the largest |dPSI| if the cluster has one, else the
    productive junction with the largest |dPSI|. Requires p.adjust < 0.05, |dPSI| > 0.05.
    Each input file is exactly one species:tissue contrast, so grouping never crosses files.
    """
    f = _as_float_flags(d)
    pl, _ = productivity_label(f["UTR"], f["Coding"], f["Annot"], f["GencodePC"])
    d = d.assign(PL=pl)
    d = d[d.PL != "NE"]
    d = d[(d["p.adjust"] < 0.05) & (d.deltapsi.abs() > 0.05)].copy()
    if not len(d):
        return d.assign(Prod=[], ContainsNC=[])
    d["Prod"] = d.PL == "PR"
    d["ContainsNC"] = d.groupby("cluster")["Prod"].transform(lambda s: (~s).any())
    d["_negabs"] = -d.deltapsi.abs()
    return (d.sort_values(["cluster", "ContainsNC", "Prod", "_negabs"])
              .drop_duplicates("cluster").drop(columns=["_negabs"]))


def spearman(x, y):
    from scipy.stats import spearmanr
    ok = ~(np.isnan(x) | np.isnan(y))
    if ok.sum() <= 2:
        return np.nan, np.nan
    return spearmanr(x[ok], y[ok])


# ---------------------------------------------------------------- reductions

def junc_class_tallies(args):
    """1.8 GB of per-junction classifications + 230 MB of counts -> ~100 tallies."""
    read_rows, uniq_rows = [], []
    for path in args.classifications:
        # .../ClassifyJuncs/<Source>/<Genome>.AllObserved._junction_classifications.txt
        source = os.path.basename(os.path.dirname(path))
        genome = os.path.basename(path).replace(
            ".AllObserved._junction_classifications.txt", "")
        d = pd.read_csv(path, sep="\t")
        for c in ["Annot", "Coding", "UTR", "GencodePC"]:
            d[c] = d[c].astype(bool)
        pl, _ = productivity_label(1.0 * d.UTR, 1.0 * d.Coding, 1.0 * d.Annot, 1.0 * d.GencodePC)
        d["PL"] = pl.values
        d["name"] = junc_name(d.Intron_coord, d.Strand)
        d = d[["name", "Annot", "PL"]].drop_duplicates("name")

        counts = pd.read_csv(os.path.join(args.count_dir, f"{genome}.uniq.junc.gz"),
                             sep="\t", header=None, usecols=[3, 4], names=["name", "score"])
        d = d.merge(counts, on="name", how="inner")
        total = d.score.sum()          # includes NE: the RPM floor is against the whole library
        d["JuncRPM"] = d.score / total * 1e6
        nn = d[d.PL != "NE"].copy()
        nn["IsProd"] = nn.PL == "PR"
        a = nn.groupby(["Annot", "IsProd"])["score"].sum().rename("n").reset_index()
        b = (nn[nn.JuncRPM > RPM_THRESH].groupby(["Annot", "IsProd"])
             .size().rename("n").reset_index())
        for df, sink, metric in ((a, read_rows, "JunctionReads"),
                                 (b, uniq_rows, "UniqueJunctions")):
            df["Metric"], df["AnnotationSource"], df["OriginGenome"] = metric, source, genome
            sink.append(df)

    out = pd.concat(read_rows + uniq_rows, ignore_index=True)
    out["Species"] = out.OriginGenome.map(species_of)
    out["Percent"] = (out.groupby(["Metric", "AnnotationSource", "Species"])["n"]
                      .transform(lambda s: s / s.sum() * 100))
    out[["Metric", "AnnotationSource", "OriginGenome", "Species", "Annot", "IsProd",
         "n", "Percent"]].to_csv(args.out, **GZ)


def transcripts_per_gene(args):
    """12 MB of transcript-biotype tables -> per-gene histogram plus per-panel totals."""
    rows = []
    for path in args.transcript_types:
        source = os.path.basename(os.path.dirname(path))
        genome = os.path.basename(path).replace(".TranscriptTypes.tsv.gz", "")
        d = pd.read_csv(path, sep="\t")
        d["AnnotationSource"], d["Species"] = source, species_of(genome)
        d["OriginGenome"] = genome
        rows.append(d)
    tt = pd.concat(rows, ignore_index=True)

    # UCSC's mouse Ensembl GTF stores transcript ids in the gene_id field; substitute the
    # corrected BED for that one annotation x genome combination.
    if args.mouse_ensembl_fix:
        mfix = pd.read_csv(args.mouse_ensembl_fix, sep="\t", header=None)[[13, 12, 15, 14]]
        mfix.columns = ["gname", "tname", "gtype", "ttype"]
        mfix["AnnotationSource"], mfix["Species"] = "Ensembl", "Mouse"
        mfix["OriginGenome"] = "Mouse_UCSC.mm39_GencodeComprehensive46"
        tt = pd.concat([tt[~((tt.AnnotationSource == "Ensembl") &
                             (tt.OriginGenome == "Mouse_UCSC.mm39_GencodeComprehensive46"))],
                        mfix], ignore_index=True)

    tt = tt[tt.gtype == "protein_coding"].copy()
    tt["IsPC"] = tt.ttype == "protein_coding"

    totals = (tt.groupby(["AnnotationSource", "Species"])
              .agg(NumCoding=("IsPC", "sum"), Total=("IsPC", "size"),
                   nGenes=("gname", "nunique")).reset_index())
    totals["NumNoncoding"] = totals.Total - totals.NumCoding
    totals[["AnnotationSource", "Species", "NumCoding", "NumNoncoding", "Total",
            "nGenes"]].to_csv(args.out_totals, **GZ)

    pg = (tt.groupby(["AnnotationSource", "Species", "gname"])["IsPC"]
          .agg(NumCoding="sum", Total="size").reset_index())
    pg["NumNoncoding"] = pg.Total - pg.NumCoding
    pg["Grp"] = pg.Total.clip(upper=10)         # ">10" bucket
    hist = (pg.groupby(["AnnotationSource", "Species", "Grp"])
            .agg(sumNC=("NumNoncoding", "sum"), sumTot=("Total", "sum"),
                 NumGenes=("gname", "size")).reset_index())
    hist["BarNoncoding"] = hist.sumNC / hist.sumTot * hist.NumGenes
    hist["BarCoding"] = hist.NumGenes - hist.BarNoncoding
    hist[["AnnotationSource", "Species", "Grp", "NumGenes", "BarCoding",
          "BarNoncoding"]].to_csv(args.out, **GZ)


def contrast_summaries(args):
    """2.7 GB across 404 contrast files -> one Spearman per contrast x productivity class."""
    rows = []
    for path in args.joined:
        base = os.path.basename(path)
        m = re.match(r"(Neonate_vs_Adult|First10_vs_Second10)\.(.+?)\.(.+?)\.joined\.tsv\.gz",
                     base)
        if m:
            fam, tissue, sp = m.groups()
            group = ("Neonate v adult within tissue" if fam == "Neonate_vs_Adult"
                     else "Early embryo v adult within tissue")
            if group.startswith("Early embryo") and sp == "Human":
                continue            # no human embryonic samples were used
            label = (f"Neonate {tissue} v Adult {tissue}" if fam == "Neonate_vs_Adult"
                     else f"Early embryo {tissue} v Adult {tissue}")
        elif re.match(r"^\w+_[\w]+-[\w]+\.joined\.tsv\.gz$", base):
            sp, pair = base.replace(".joined.tsv.gz", "").split("_", 1)
            tA, tB = pair.split("-")
            group, label = "Adult tissues", f"Adult {tA} v Adult {tB}"
        else:
            # an alternative stage split that the figures do not use
            continue

        d = one_event_per_cluster(pd.read_csv(path, sep="\t", usecols=DS_COLS,
                                              low_memory=False))
        if not len(d):
            continue
        for prod in (True, False):
            s = d[d.Prod == prod]
            rho, p = spearman(s.deltapsi.values, s.logFC.values)
            rows.append(dict(ContrastGroup=group, Tissue=label, Species=sp,
                             Coding="Productive" if prod else "Unproductive",
                             n=len(s), spearman=rho, p=p))
    out = pd.DataFrame(rows)
    out["Tissue"] = out.Tissue.replace(
        {"Neonate Testis v Adult Testis": "Juvenile Testis vs Adult Testis"})
    out.to_csv(args.out, **GZ)


def betabeta_points(args):
    """One adult tissue pair, per-cluster scatter points for the main beta-beta panel."""
    rows = []
    for path in args.joined:
        sp = os.path.basename(path).split("_")[0]
        d = one_event_per_cluster(pd.read_csv(path, sep="\t", usecols=DS_COLS,
                                             low_memory=False))
        rows.append(d.assign(Species=sp)[["Species", "Prod", "cluster", "deltapsi", "logFC"]])
    pd.concat(rows, ignore_index=True).to_csv(args.out, **GZ)


def junc_class_by_sample(args):
    """Per-library junction-productivity read counts, consolidated across species."""
    rows = []
    for path in args.summaries:
        d = pd.read_csv(path, sep="\t", header=None,
                        names=["fn", "Annot", "Coding", "UTR", "GencodePC", "TotalCounts"])
        d["ID"] = d.fn.str.replace(
            r"rna-seq/SplicingAnalysis/juncfiles/(.+?)\.junccounts\.tsv\.gz", r"\1", regex=True)
        rows.append(d.drop(columns=["fn"]))
    pd.concat(rows, ignore_index=True)[
        ["ID", "Annot", "Coding", "UTR", "GencodePC", "TotalCounts"]].to_csv(args.out, **GZ)


def testis_percent_up(args):
    """Percent unproductive reads per library (the Fig 3c metric)."""
    d = pd.read_csv(args.counts, sep="\t")
    pl, _ = productivity_label(1.0 * d.UTR, 1.0 * d.Coding, 1.0 * d.Annot, 1.0 * d.GencodePC)
    d["PL"] = pl.values
    d = d[d.PL != "NE"]
    tot = d.groupby(["ID", "PL"])["TotalCounts"].sum().unstack(fill_value=0)
    tot["Percent"] = tot["UP"] / (tot["UP"] + tot["PR"]) * 100
    tot.reset_index()[["ID", "UP", "PR", "Percent"]].to_csv(args.out, **GZ)


def scn8a_phylop(args):
    """Per-base PhyloP restricted to the three plotted SCN8A exons."""
    bg = pd.read_csv(args.bedgraph, sep="\t", names=["chrom", "start", "stop", "PhyloP"])
    trios = pd.read_csv(args.dnds, sep="\t")
    mphy = _mean_phylop(args.mean_phylop)

    out = []
    for label, region in zip(args.exon_labels, args.exon_regions):
        chrom, rng = region.split(":")
        lo, hi = (int(v) for v in rng.split("-"))
        w = bg[(bg.stop >= lo) & (bg.stop <= hi)].copy()
        w["Exon"] = label
        w["CodonPos"] = (w.stop - lo) % 3
        mp = mphy.loc[mphy.region == region, "MeanPhyloP"]
        w["MeanPhyloP_bigwig"] = float(mp.iloc[0]) if len(mp) else np.nan
        dn = _dnds_for_region(trios, region)
        w["dNdS"] = dn
        out.append(w[["Exon", "chrom", "stop", "PhyloP", "CodonPos",
                      "MeanPhyloP_bigwig", "dNdS"]])
    pd.concat(out, ignore_index=True).to_csv(args.out, **GZ)


def _mean_phylop(path):
    m = pd.read_csv(path, sep="\t", names=["chrom", "bedstart", "bedstop", "MeanPhyloP"])
    m["region"] = (m.chrom + ":" + (m.bedstart + 1).astype(str) + "-"
                   + m.bedstop.astype(str))
    return m[["region", "MeanPhyloP"]]


def _dnds_for_region(trios, region):
    """dN/dS for a region, whichever exon slot of the trio it occupies."""
    for col, val in [("cassette_exon", "dn_ds_ratio_ce"),
                     ("upstream_exon", "dn_ds_ratio_us"),
                     ("downstream_exon", "dn_ds_ratio_ds")]:
        hit = trios.loc[trios[col] == region, val].dropna()
        if len(hit):
            return float(hit.iloc[0])
    return np.nan


def dnds_phylop_trios(args):
    """dN/dS and mean PhyloP per exon of each cassette-exon trio, one row per exon."""
    tidy = _dnds_tidy(args)
    parts = [
        tidy.assign(PhyloP=tidy.MeanPhyloP_exon_ce, dnds=tidy.dn_ds_ratio_ce,
                    Type="CassetteExon", BroadType="Cassette exon"),
        tidy.assign(PhyloP=tidy.MeanPhyloP_exon_us, dnds=tidy.dn_ds_ratio_us,
                    Type="Downstream", BroadType="Flanking exons"),
        tidy.assign(PhyloP=tidy.MeanPhyloP_exon_ds, dnds=tidy.dn_ds_ratio_ds,
                    Type="UpstreamExon", BroadType="Flanking exons"),
    ]
    out = pd.concat(parts, ignore_index=True)
    out["Productivity"] = out.Productivity.map({"PR": "Productive", "UP": "Unproductive"})
    names = pd.read_csv(args.gene_names, sep="\t")
    names = (names[names.Species == "Human"]
             [["ensembl_gene_id_version", "external_gene_name"]].drop_duplicates())
    names.columns = ["gene_name", "external_gene_name"]
    out = out.merge(names, on="gene_name", how="left")
    out[["rn", "gene_name", "external_gene_name", "BroadType", "Type", "Productivity",
         "DifferentiallySpliced", "PhyloP", "dnds"]].to_csv(args.out, **GZ)


def _dnds_tidy(args):
    """One row per human cluster x productivity, with dN/dS and mean PhyloP per exon slot."""
    trios = pd.read_csv(args.dnds, sep="\t")
    mphy = _mean_phylop(args.mean_phylop)
    with np.errstate(divide="ignore"):
        trios["MeanDnDs.Flank"] = np.exp(
            (np.log(trios.dn_ds_ratio_us) + np.log(trios.dn_ds_ratio_ds)) / 2)
    # prefer the differentially spliced entry, then the least-constrained flanks
    trios = (trios.sort_values(["cassette_exon", "DifferentiallySpliced", "MeanDnDs.Flank"],
                               ascending=[True, False, True])
                  .drop_duplicates("cassette_exon"))
    long = []
    for rt, rcol in [("exon_us", "upstream_exon"), ("exon_ds", "downstream_exon"),
                     ("exon_ce", "cassette_exon")]:
        t = trios[[rcol]].copy()
        t.columns = ["region"]
        t["region_type"] = rt
        long.append(t)
    long = pd.concat(long).reset_index().rename(columns={"index": "_i"})
    long = long.merge(mphy, on="region", how="inner")
    wide = long.pivot(index="_i", columns="region_type", values=["region", "MeanPhyloP"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    tidy = trios.join(wide, how="inner")
    tidy["MeanPhyloP.Flank"] = (tidy.MeanPhyloP_exon_us + tidy.MeanPhyloP_exon_ds) / 2
    tidy = (tidy.sort_values(["human.cluster", "Productivity", "MeanPhyloP_exon_ce"],
                             ascending=[True, True, False])
                .drop_duplicates(["human.cluster", "Productivity"]))
    tidy["rn"] = np.arange(len(tidy))
    return tidy


def cassette_exon_supplement(args):
    """Supplementary table of cassette-exon trios with conservation and constraint.

    `IsInBoxedRegion` marks the quadrant highlighted in the dN/dS vs PhyloP panel:
    unproductive cassette exons that are nucleotide-conserved but not amino-acid constrained.
    """
    tidy = _dnds_tidy(args)
    names = pd.read_csv(args.gene_names, sep="\t")
    names = (names[names.Species == "Human"]
             [["ensembl_gene_id_version", "external_gene_name"]].drop_duplicates())
    names.columns = ["gene_name", "external_gene_name"]
    out = tidy.merge(names, on="gene_name", how="left")
    out["IsInBoxedRegion"] = ((out.Productivity == "UP") & (out.MeanPhyloP_exon_ce > 4)
                              & (out.dn_ds_ratio_ce > 0.5) & (out.dn_ds_ratio_ce < 50))
    out = out.rename(columns={"rn": "row_number", "flanking_int_chrom": "chrom",
                              "dn_ds_ratio_us": "dnds_us", "dn_ds_ratio_ds": "dnds_ds",
                              "dn_ds_ratio_ce": "dnds_ce"})
    out = out.drop(columns=[c for c in ["MatchedFlankingJuncDs", "upstream_exon",
                                        "cassette_exon", "downstream_exon"]
                            if c in out.columns])
    lead = ["row_number", "gene_name", "external_gene_name", "Productivity",
            "human.cluster", "chrom", "strand"]
    lead = [c for c in lead if c in out.columns]
    out[lead + [c for c in out.columns if c not in lead]].to_csv(args.out, **GZ)


def arhgap17_flanks(args):
    """Mean PhyloP over the 50 bp flanking each ARHGAP17 exon (bed -> tsv.gz)."""
    pd.read_csv(args.bed, sep="\t",
                names=["chrom", "start", "stop", "meanPhyloP", "score", "strand"]
                ).to_csv(args.out, **GZ)


def nmd_expression(args):
    """119 MB of expression + 12 MB of ortholog tables -> core NMD gene deltas per library."""
    CORE = ["UPF3A", "UPF3B", "SMG7", "SMG6", "SMG5", "SMG1", "UPF2", "UPF1"]
    orth = pd.read_csv(args.orthologs, sep="\t")
    genes = pd.read_csv(args.transcript_info, sep="\t")

    core_hs = (genes[(genes.Species == "Human") & genes.external_gene_name.isin(CORE)]
               .drop_duplicates("ensembl_gene_id_version")
               .rename(columns={"ensembl_gene_id": "hsapiens_homolog_ensembl_gene"}))
    o = orth.merge(core_hs[["hsapiens_homolog_ensembl_gene"]].drop_duplicates(),
                   on="hsapiens_homolog_ensembl_gene")
    o = (o.sort_values(["Species", "hsapiens_homolog_ensembl_gene",
                        "hsapiens_homolog_perc_id"])
          .drop_duplicates(["Species", "hsapiens_homolog_ensembl_gene"]))
    hs_self = core_hs[["hsapiens_homolog_ensembl_gene", "external_gene_name"]].copy()
    hs_self["Species"] = "Human"
    hs_self["ensembl_gene_id"] = hs_self.hsapiens_homolog_ensembl_gene
    core = pd.concat([o, hs_self], ignore_index=True).merge(
        genes[["ensembl_gene_id", "ensembl_gene_id_version"]].drop_duplicates(),
        on="ensembl_gene_id")

    # Rat's RefSeq expression matrix is keyed on gene symbol, so no Ensembl id matches.
    RAT = ["Upf1", "Upf2", "Upf3a", "Upf3b", "Smg1", "Smg5", "Smg6", "Smg7"]
    rat = pd.DataFrame({"ensembl_gene_id": RAT, "ensembl_gene_id_version": RAT,
                        "external_gene_name": [s.upper() for s in RAT], "Species": "Rat"})
    rat = rat.merge(core_hs[["hsapiens_homolog_ensembl_gene", "external_gene_name"]]
                    .drop_duplicates(), on="external_gene_name")
    hs_name = (core_hs[["hsapiens_homolog_ensembl_gene", "external_gene_name"]]
               .drop_duplicates().rename(columns={"external_gene_name": "human.name"}))
    core = pd.concat([core[core.Species != "Rat"], rat], ignore_index=True).merge(
        hs_name, on="hsapiens_homolog_ensembl_gene")

    wanted = set(core.ensembl_gene_id_version)
    rows = []
    for path in args.expression:
        genome = os.path.basename(path).replace(".log2rpkm.tsv.gz", "")
        d = pd.read_csv(path, sep="\t")
        d = d[d.Geneid.isin(wanted)].melt(id_vars=["Geneid"], var_name="ID",
                                          value_name="log2RPKM")
        d["OriginGenome"] = genome
        rows.append(d)
    rpkm = pd.concat(rows, ignore_index=True).merge(
        core, left_on="Geneid", right_on="ensembl_gene_id_version")

    su = _samples_by_tissue(args.sample_stages)
    dat = rpkm.merge(su[["ID", "Species", "Tissue", "PercentRank_OrdinalStage"]],
                     on=["ID", "Species"], how="inner")
    t0 = (dat.assign(_m=dat.groupby(["OriginGenome", "Geneid", "Tissue"])
                            ["PercentRank_OrdinalStage"].transform("min"))
             .query("PercentRank_OrdinalStage == _m")
             .groupby(["OriginGenome", "Geneid", "Tissue"])["log2RPKM"].mean()
             .rename("MeanT0").reset_index())
    dat = dat.merge(t0, on=["OriginGenome", "Geneid", "Tissue"])
    dat["Delta"] = dat.log2RPKM - dat.MeanT0
    dat[["ID", "Species", "Tissue", "human.name", "PercentRank_OrdinalStage",
         "log2RPKM", "Delta"]].to_csv(args.out, **GZ)


def _samples_by_tissue(sample_stages):
    sw = pd.read_csv(sample_stages, sep="\t").rename(columns={"ID_Species": "Species"})
    sw["PercentRank_OrdinalStage"] = (sw.groupby("Species")["Ordinal_stage"]
                                      .transform(lambda s: s.rank(method="min").sub(1)
                                                 / (len(s) - 1)))
    su = sw[sw.NoHumanEmbryo].copy()
    su = su.assign(Tissue=su.Tissue_ForDevelopementalAnalysis.astype(str)
                   .str.split(",")).explode("Tissue")
    su["Tissue"] = su.Tissue.str.strip()
    return su[su.Tissue.ne("nan") & su.Tissue.notna()]


def developmental_markers(args):
    """Percent-rank position of each developmental landmark, per species."""
    sw = pd.read_csv(args.sample_stages, sep="\t").rename(columns={"ID_Species": "Species"})
    sw["PercentRank_OrdinalStage"] = (sw.groupby("Species")["Ordinal_stage"]
                                      .transform(lambda s: s.rank(method="min").sub(1)
                                                 / (len(s) - 1)))
    stages = pd.read_csv(args.stages, sep="\t")
    sw = sw.merge(stages[["Ordinal_stage", "Species", "Marker"]].dropna(subset=["Marker"])
                        .drop_duplicates(), on=["Ordinal_stage", "Species"], how="left")
    out = (sw[sw.Marker.notna()].groupby(["Species", "Marker"])
           ["PercentRank_OrdinalStage"].min().unstack()
           .rename(columns={"Onset meiosis spermatogenesis": "OnsetMeiosis"}).reset_index())
    out.to_csv(args.out, **GZ)


def noembryo_samples(args):
    """Flag which libraries are used for which developmental analyses."""
    samples = pd.read_csv(args.samples, sep="\t")
    stages = pd.read_csv(args.stages, sep="\t")
    birth = (stages[stages.Marker == "Birth"][["Species", "Ordinal_stage"]]
             .drop_duplicates().rename(columns={"Species": "ID_Species",
                                                "Ordinal_stage": "BirthStage"}))
    keep = ["ID", "ID_Species", "Tissue_ForDevelopementalAnalysis", "Ordinal_stage",
            "Used library?"]
    d = samples[keep].merge(birth, on="ID_Species", how="left")
    used = d["Used library?"].isin(["Yes", "yes"])
    d["NoEmbryo"] = used & (d.Ordinal_stage >= d.BirthStage)
    d["NoHumanEmbryo"] = used & (((d.Ordinal_stage >= d.BirthStage)
                                  & (d.ID_Species == "Human"))
                                 | (d.ID_Species != "Human"))
    allsamples = pd.read_csv(args.star_samples, sep="\t")[["SampleName", "STARGenomeName"]]
    d = d.merge(allsamples.rename(columns={"SampleName": "ID",
                                           "STARGenomeName": "GenomeName"}),
                on="ID", how="left")
    d.to_csv(args.out, **GZ)


# ---------------------------------------------------------------- CLI

def conserved_as_examples(args):
    """PSI, cluster depth and host-gene expression at the plotted example loci.

    Slices ~1.5 GB of PSI tables down to a few thousand rows. Rows are streamed with
    `zcat | awk` rather than tabix because the tabix indexes are stale for several species.
    """
    import io
    import subprocess

    loci = pd.read_csv(args.loci, sep="\t")
    groups = pd.read_csv(args.sample_groups, sep="\t")

    def slice_junc(genome, junccoord, fname):
        chrom, rng = junccoord.split(":")
        start, end = (int(v) for v in rng.split("-"))
        path = os.path.join(args.psi_dir, genome, "juncTableBeds", fname)
        cmd = (f"zcat '{path}' | awk -F'\\t' "
               f"'NR==1 || ($1==\"{chrom}\" && $2=={start} && $3=={end + 1})'")
        txt = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True,
                             check=True).stdout
        d = pd.read_csv(io.StringIO(txt), sep="\t")
        d = d.rename(columns={d.columns[0]: "chrom"})
        return d.melt(id_vars=list(d.columns[:6]), var_name="ID",
                      value_name="value")[["ID", "value"]]

    rows = []
    for _, r in loci.iterrows():
        psi = slice_junc(r.OriginGenome, r.JuncCoord,
                         "PSI.sorted.bed.gz").rename(columns={"value": "PSI"})
        den = slice_junc(r.OriginGenome, r.JuncCoord,
                         "PSIDenom.sorted.bed.gz").rename(columns={"value": "ReadsInCluster"})
        expr = pd.read_csv(os.path.join(args.expression_dir,
                                        f"{r.OriginGenome}.log2rpkm.tsv.gz"), sep="\t")
        expr = expr[expr.Geneid == r.Geneid].melt(id_vars=["Geneid"], var_name="ID",
                                                  value_name="log2RPKM")
        d = psi.merge(den, on="ID").merge(expr[["ID", "log2RPKM"]], on="ID", how="left")
        d["gene_name"], d["SpeciesShort"] = r.gene_name, r.SpeciesShort
        d["Species"], d["JuncCoord"], d["Geneid"] = r.OriginGenome, r.JuncCoord, r.Geneid
        rows.append(d)

    out = pd.concat(rows, ignore_index=True)
    # the grouping is gene-specific: brain genes are grouped against brain samples, testis
    # genes against testis samples, so join on gene as well as library
    out = out.merge(groups, on=["gene_name", "SpeciesShort", "ID"], how="inner")
    out["ReadsInClusterGroup"] = out.ReadsInCluster >= args.min_reads
    out = out.rename(columns={"SpeciesShort": "Species.short"})
    out[["JuncCoord", "ID", "PSI", "Species", "Species.short", "ReadsInCluster",
         "ReadsInClusterGroup", "NewerGroup", "gene_name", "Geneid",
         "log2RPKM"]].to_csv(args.out, **GZ)


def contrast_group_files(args):
    """Write one leafcutter group file per contrast, from the supplementary contrast table.

    The table is the authority for which contrasts the paper used and which libraries are in
    each arm, so the group files leafcutter consumes are derived from it rather than shipped
    independently. Each file is `<sample> <group>`, space separated, no header.

    Contrast names map onto filenames as:
      Neonate <T> v Adult <T>        -> ContrastGroupFiles/Neonate_vs_Adult.<T>.<Sp>.txt
      Juvenile Testis v Adult Testis -> ContrastGroupFiles/Neonate_vs_Adult.Testis.<Sp>.txt
      Early embryo <T> v Adult <T>   -> ContrastGroupFiles/First10_vs_Second10.<T>.<Sp>.txt
      Adult <A> v Adult <B>          -> ContrastAdultTissueGroupFiles/<Sp>_<A>-<B>.txt
    """
    tab = pd.read_csv(args.table, sep="\t")
    os.makedirs(args.dev_dir, exist_ok=True)
    os.makedirs(args.adult_dir, exist_ok=True)

    written = []
    for _, r in tab.iterrows():
        name, sp = r.ContrastName, r.Species
        m_adult = re.match(r"^Adult (\w+) v Adult (\w+)$", name)
        m_neo = re.match(r"^Neonate (\w+) v Adult \w+$", name)
        m_juv = re.match(r"^Juvenile (\w+) v Adult \w+$", name)
        m_emb = re.match(r"^Early embryo (\w+) v Adult \w+$", name)
        if m_adult:
            tA, tB = m_adult.groups()
            path = os.path.join(args.adult_dir, f"{sp}_{tA}-{tB}.txt")
            labels = (tA, tB)
        elif m_neo or m_juv:
            tissue = (m_neo or m_juv).group(1)
            path = os.path.join(args.dev_dir, f"Neonate_vs_Adult.{tissue}.{sp}.txt")
            labels = ("Early", "Late")
        elif m_emb:
            tissue = m_emb.group(1)
            path = os.path.join(args.dev_dir, f"First10_vs_Second10.{tissue}.{sp}.txt")
            labels = ("Early", "Late")
        else:
            sys.exit(f"unrecognised contrast name: {name!r}")

        lines = []
        for col, label in zip(["SamplesInGroup1", "SamplesInGroup2"], labels):
            if isinstance(r[col], str):
                lines += [f"{s} {label}" for s in r[col].split(", ")]
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        written.append(path)
    print(f"wrote {len(written)} contrast group files", file=sys.stderr)


def join_by_human_coords(args):
    """Join every species' differential-splicing results onto human by lifted coordinates.

    For each contrast, a non-human junction is matched to the human junction it lifts over
    to, and the two sets of statistics are placed side by side. Only the columns the
    downstream cascade needs are kept, which keeps this intermediate to tens of MB rather
    than the several hundred it would otherwise be.
    """
    COLS = ["cluster", "deltapsi", "Intron_coord", "p", "p.adjust",
            "Annot", "Coding", "UTR", "GencodePC"]

    def load(path):
        base = os.path.basename(path)
        m = re.match(r"(Neonate_vs_Adult|First10_vs_Second10)\.(.+?)\.(.+?)\.joined\.tsv\.gz$",
                     base)
        if m:
            fam, tissue, sp = m.groups()
            group = ("Neonate v adult within tissue" if fam == "Neonate_vs_Adult"
                     else "Early embryo v adult within tissue")
            label = tissue
        elif re.match(r"^\w+_[\w]+-[\w]+\.joined\.tsv\.gz$", base):
            sp, pair = base.replace(".joined.tsv.gz", "").split("_", 1)
            group, label = "Adult tissues", pair
        else:
            return None
        d = pd.read_csv(path, sep="\t", usecols=COLS, low_memory=False)
        f = _as_float_flags(d)
        pl, _ = productivity_label(f["UTR"], f["Coding"], f["Annot"], f["GencodePC"])
        d["ProductivityClass"] = pl
        d["Species"], d["ContrastGroup"], d["Tissue"] = sp, group, label
        # the junction name, so the lifted coordinates can be matched
        d["name"] = junc_name(d.Intron_coord,
                              d.cluster.str[-1] if d.cluster.dtype == object else "+")
        return d[["Species", "ContrastGroup", "Tissue", "cluster", "Intron_coord", "name",
                  "deltapsi", "p", "p.adjust", "ProductivityClass"]]

    frames = [f for f in (load(p) for p in args.joined) if f is not None]
    allds = pd.concat(frames, ignore_index=True)
    del frames; gc.collect()

    lift = []
    BEDCOLS = ["chrom", "start", "end", "name", "score", "strand", "thickStart",
               "thickStop", "color", "nBlocks", "blockSizes", "blockStarts"]
    for path in args.lifted:
        d = pd.read_csv(path, sep="\t", header=None, names=BEDCOLS,
                        usecols=["chrom", "start", "end", "name", "strand"])
        d["human.name"] = (d.chrom + "_" + (d.start + 1).astype(str) + "_"
                           + (d.end - 1).astype(str) + "_" + d.strand)
        lift.append(d[["name", "human.name"]])
    lift = pd.concat(lift, ignore_index=True).drop_duplicates()

    human = allds[allds.Species == "Human"].copy()
    other = allds[allds.Species != "Human"].merge(lift, on="name", how="inner")
    del allds; gc.collect()

    out = other.merge(human.rename(columns={"name": "human.name"}),
                      on=["ContrastGroup", "Tissue", "human.name"],
                      suffixes=(".OtherSpecies", ".Human"), how="inner")
    out = out.rename(columns={
        "Species.OtherSpecies": "QuerySpecies",
        "deltapsi.OtherSpecies": "QuerySpecies.leaf_deltapsi",
        "p.OtherSpecies": "QuerySpecies.leaf_nominal_p",
        "ProductivityClass.OtherSpecies": "QuerySpecies.ProductivityClass",
        "Intron_coord.Human": "Human.Intron_coordinates",
        "cluster.Human": "Human.leaf_cluster_name",
        "deltapsi.Human": "Human.leaf_deltapsi",
        "p.adjust.Human": "Human.leaf_adjusted_p",
        "ProductivityClass.Human": "Human.ProductivityClass"})
    keep = ["QuerySpecies", "ContrastGroup", "Tissue", "QuerySpecies.leaf_deltapsi",
            "QuerySpecies.leaf_nominal_p", "QuerySpecies.ProductivityClass",
            "Human.Intron_coordinates", "Human.leaf_cluster_name", "Human.leaf_deltapsi",
            "Human.leaf_adjusted_p", "Human.ProductivityClass"]
    out[[c for c in keep if c in out.columns]].to_csv(args.out, **GZ)
    print(f"joined {len(out)} rows", file=sys.stderr)


def _cascade_frames(args):
    """Shared inputs for the cascade counts and the conserved-junction supplement.

    Returns the lifted+classified junction table, the clustered junction list, and the
    human-only view of it. Two details matter for the counts: classification rows are NOT
    deduplicated (a junction assigned to several overlapping genes contributes several rows),
    and junctions must also appear in the observed-junction set.
    """
    BEDCOLS = ["chrom", "start", "end", "name", "score", "strand", "thickStart",
               "thickStop", "color", "nBlocks", "blockSizes", "blockStarts"]
    lift = []
    for path in args.lifted:
        genome = os.path.basename(path).replace(".Lifted.bed.gz", "")
        d = pd.read_csv(path, sep="\t", header=None, names=BEDCOLS,
                        usecols=["chrom", "start", "end", "name", "strand"])
        d["OriginGenome"], d["Species"] = genome, species_of(genome)
        d["human.name"] = (d.chrom + "_" + (d.start + 1).astype(str) + "_"
                           + (d.end - 1).astype(str) + "_" + d.strand)
        lift.append(d[["OriginGenome", "Species", "name", "human.name"]])
    lift = pd.concat(lift, ignore_index=True)

    jp = {}
    for path in args.classifications:
        genome = os.path.basename(path).replace(
            ".AllObserved._junction_classifications.txt", "")
        d = pd.read_csv(path, sep="\t")
        for c in ["Coding", "GencodePC"]:
            d[c] = d[c].astype(bool)
        # coarser than the UP/PR/NE label: near-UTR junctions count as unproductive here
        d["Coding"] = d.Coding | d.GencodePC
        d["name"] = junc_name(d.Intron_coord, d.Strand)
        obs = set(pd.read_csv(os.path.join(args.observed_dir,
                                           f"{genome}.uniq.annotated.tsv.gz"),
                              sep="\t", usecols=["name"]).name)
        jp[species_of(genome)] = d.loc[d.name.isin(obs), ["name", "Coding"]]

    human = jp["Human"].rename(columns={"name": "human.name", "Coding": "Coding.Human"})
    other = pd.concat([v.assign(Species=k).rename(columns={"Coding": "Coding.OtherSpecies"})
                       for k, v in jp.items() if k != "Human"], ignore_index=True)
    lift_up = (lift.merge(other, on=["Species", "name"], how="inner")
                   .merge(human, on="human.name", how="inner"))

    clus = []
    for path in args.clustered:
        genome = os.path.basename(path).replace(".Clustered.juncList.tsv.gz", "")
        d = pd.read_csv(path, sep="\t")
        ex = d.chrom.str.extract(r"^(.+?):(.+?):(.+?):(.+_[+-])$")
        d["name"] = ex[0] + "_" + ex[1] + "_" + ex[2] + "_" + ex[3].str[-1]
        d["cluster"] = ex[3] + "_"
        d["OriginGenome"] = genome
        clus.append(d[["OriginGenome", "name", "cluster", "chrom"]])
    clus = pd.concat(clus, ignore_index=True)
    hum_clu = (clus[clus.OriginGenome.str.startswith("Human_")]
               .rename(columns={"name": "human.name", "cluster": "human.cluster",
                                "chrom": "human.leaf2Name"})
               .drop(columns=["OriginGenome"]))
    return {"lift_up": lift_up, "clus": clus.drop(columns=["chrom"]), "hum_clu": hum_clu}


def _conserved_events(path):
    """The cross-species join, reduced to ONE row per (species, human cluster)."""
    j = pd.read_csv(path, sep="\t", low_memory=False, usecols=[
        "QuerySpecies", "QuerySpecies.leaf_deltapsi", "QuerySpecies.leaf_nominal_p",
        "QuerySpecies.ProductivityClass", "Human.Intron_coordinates",
        "Human.leaf_cluster_name", "Human.leaf_deltapsi", "Human.leaf_adjusted_p",
        "Human.ProductivityClass"]).rename(columns={
        "QuerySpecies": "Species", "QuerySpecies.leaf_deltapsi": "dpsi.Other",
        "QuerySpecies.leaf_nominal_p": "p.Other",
        "QuerySpecies.ProductivityClass": "PC.Other",
        "Human.leaf_cluster_name": "cluster.Human",
        "Human.leaf_deltapsi": "dpsi.Human", "Human.leaf_adjusted_p": "q.Human",
        "Human.ProductivityClass": "PC.Human"})
    j["Coding.Other"] = j["PC.Other"] == "PR"
    j["Coding.Human"] = j["PC.Human"] == "PR"
    j = j[(j["q.Human"] < 0.05) & (j["dpsi.Human"].abs() > 0.05)
          & (j["p.Other"] < 0.05)].copy()
    ic = j["Human.Intron_coordinates"].str.rsplit(":", n=1, expand=True)
    rng = ic[1].str.rsplit("-", n=1, expand=True)
    j["human.name"] = (ic[0] + "_" + rng[0] + "_" + rng[1] + "_"
                       + j["cluster.Human"].str[-1])
    j["UP_in_both"] = (~j["Coding.Other"]) & (~j["Coding.Human"])
    j["IsMatchingSigns"] = np.sign(j["dpsi.Human"]) == np.sign(j["dpsi.Other"])
    return (j.sort_values(["Species", "cluster.Human", "IsMatchingSigns", "UP_in_both",
                          "Coding.Other"], ascending=[True, True, False, False, True])
             .drop_duplicates(["Species", "cluster.Human"]))


def conserved_event_cascade(args):
    """Conserved-event counts after each successive filter, unproductive and productive.

    Collapses ~2.3 GB into 70 numbers. `*Any` counts distinct human features across all
    non-macaque species, so it is not the sum of the per-species bars.
    """
    cas = _cascade_frames(args)
    lift_up, clus, hum_clu = cas["lift_up"], cas["clus"], cas["hum_clu"]
    j = _conserved_events(args.joined_by_human)

    def counts(df, key, dedup_species=False):
        """Per-species counts plus `*Any`.

        dedup_species=True collapses to one row per (key, Species) first, which the
        cluster-level step needs; the junction-level steps count rows directly.
        """
        base = df.drop_duplicates([key, "Species"]) if dedup_species else df
        out = base.groupby("Species").size().to_dict()
        out["*Any"] = df[df.Species != "Macaque"].drop_duplicates(key).shape[0]
        return out

    rows = []
    f1 = counts(lift_up, "human.name")
    for path_name, is_prod in [("Unproductive", False), ("Productive", True)]:
        s2 = lift_up[lift_up["Coding.OtherSpecies"] == is_prod]
        s3 = s2[s2["Coding.Human"] == is_prod]
        s4 = (s3.merge(clus, on=["OriginGenome", "name"], how="inner")
                .merge(hum_clu, on="human.name", how="inner"))
        ev = (j[j.UP_in_both & j.IsMatchingSigns] if not is_prod
              else j[j["Coding.Other"] & j["Coding.Human"] & j.IsMatchingSigns])
        steps = [f1, counts(s2, "human.name"), counts(s3, "human.name"),
                 counts(s4, "human.cluster", dedup_species=True),
                 counts(ev, "cluster.Human")]
        for i, st in enumerate(steps, 1):
            for sp, n in st.items():
                rows.append(dict(Path=path_name, FilterStep=i, Species=sp, n=n))
    pd.DataFrame(rows).to_csv(args.out, **GZ)


def conserved_juncs_past_macaque(args):
    """Supplementary table: conserved junctions shared with a non-macaque species.

    One row per human junction x productivity class, flagged for whether the event is
    differentially spliced with a consistent direction in both species.
    """
    cas = _cascade_frames(args)
    lift_up, clus, hum_clu = cas["lift_up"], cas["clus"], cas["hum_clu"]

    rows = []
    for label, is_prod in [("PR", True), ("UP", False)]:
        s = lift_up[(lift_up["Coding.OtherSpecies"] == is_prod)
                    & (lift_up["Coding.Human"] == is_prod)]
        s = (s.merge(clus, on=["OriginGenome", "name"], how="inner")
              .merge(hum_clu, on="human.name", how="inner"))
        s = s[s.Species != "Macaque"]
        s = s[["human.leaf2Name", "human.name", "human.cluster"]].drop_duplicates()
        s["Productivity"] = label
        rows.append(s)
    out = pd.concat(rows, ignore_index=True)

    # the flag comes from the conserved-event table, which is already one row per
    # (species, human cluster) -- dedup first, then ask about matching signs
    j = _conserved_events(args.joined_by_human)
    j = j[j.IsMatchingSigns & (j.Species != "Macaque")]
    ds = pd.DataFrame({"human.name": j["human.name"].unique()})
    ds["DifferentiallySpliced"] = True

    out = out.merge(ds, on="human.name", how="left")
    out["DifferentiallySpliced"] = out.DifferentiallySpliced.fillna(False)
    out[["human.name", "human.leaf2Name", "human.cluster", "Productivity",
         "DifferentiallySpliced"]].to_csv(args.out, **GZ)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("junc-class-tallies")
    p.add_argument("--classifications", nargs="+", required=True)
    p.add_argument("--count-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=junc_class_tallies)

    p = sub.add_parser("transcripts-per-gene")
    p.add_argument("--transcript-types", nargs="+", required=True)
    p.add_argument("--mouse-ensembl-fix")
    p.add_argument("--out", required=True)
    p.add_argument("--out-totals", required=True)
    p.set_defaults(func=transcripts_per_gene)

    p = sub.add_parser("contrast-summaries")
    p.add_argument("--joined", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=contrast_summaries)

    p = sub.add_parser("betabeta-points")
    p.add_argument("--joined", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=betabeta_points)

    p = sub.add_parser("junc-class-by-sample")
    p.add_argument("--summaries", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=junc_class_by_sample)

    p = sub.add_parser("testis-percent-up")
    p.add_argument("--counts", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=testis_percent_up)

    p = sub.add_parser("scn8a-phylop")
    p.add_argument("--bedgraph", required=True)
    p.add_argument("--dnds", required=True)
    p.add_argument("--mean-phylop", required=True)
    p.add_argument("--exon-labels", nargs="+", required=True)
    p.add_argument("--exon-regions", nargs="+", required=True)
    p.add_argument("--gene", default="SCN8A")
    p.add_argument("--out", required=True)
    p.set_defaults(func=scn8a_phylop)

    p = sub.add_parser("dnds-phylop-trios")
    p.add_argument("--dnds", required=True)
    p.add_argument("--mean-phylop", required=True)
    p.add_argument("--gene-names", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=dnds_phylop_trios)

    p = sub.add_parser("conserved-as-examples")
    p.add_argument("--loci", required=True)
    p.add_argument("--sample-groups", required=True)
    p.add_argument("--psi-dir", required=True)
    p.add_argument("--expression-dir", required=True)
    p.add_argument("--min-reads", type=int, default=10)
    p.add_argument("--out", required=True)
    p.set_defaults(func=conserved_as_examples)

    p = sub.add_parser("contrast-group-files")
    p.add_argument("--table", required=True)
    p.add_argument("--dev-dir", required=True)
    p.add_argument("--adult-dir", required=True)
    p.set_defaults(func=contrast_group_files)

    p = sub.add_parser("join-by-human-coords")
    p.add_argument("--joined", nargs="+", required=True)
    p.add_argument("--lifted", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=join_by_human_coords)

    p = sub.add_parser("conserved-event-cascade")
    p.add_argument("--lifted", nargs="+", required=True)
    p.add_argument("--classifications", nargs="+", required=True)
    p.add_argument("--clustered", nargs="+", required=True)
    p.add_argument("--observed-dir", required=True)
    p.add_argument("--joined-by-human", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=conserved_event_cascade)

    p = sub.add_parser("conserved-juncs-past-macaque")
    p.add_argument("--lifted", nargs="+", required=True)
    p.add_argument("--classifications", nargs="+", required=True)
    p.add_argument("--clustered", nargs="+", required=True)
    p.add_argument("--observed-dir", required=True)
    p.add_argument("--joined-by-human", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=conserved_juncs_past_macaque)

    p = sub.add_parser("cassette-exon-supplement")
    p.add_argument("--dnds", required=True)
    p.add_argument("--mean-phylop", required=True)
    p.add_argument("--gene-names", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cassette_exon_supplement)

    p = sub.add_parser("arhgap17-flanks")
    p.add_argument("--bed", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=arhgap17_flanks)

    p = sub.add_parser("nmd-expression")
    p.add_argument("--expression", nargs="+", required=True)
    p.add_argument("--orthologs", required=True)
    p.add_argument("--transcript-info", required=True)
    p.add_argument("--sample-stages", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=nmd_expression)

    p = sub.add_parser("developmental-markers")
    p.add_argument("--sample-stages", required=True)
    p.add_argument("--stages", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=developmental_markers)

    p = sub.add_parser("noembryo-samples")
    p.add_argument("--samples", required=True)
    p.add_argument("--stages", required=True)
    p.add_argument("--star-samples", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=noembryo_samples)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
