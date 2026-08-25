#!/usr/bin/env python3
"""Repair UCSC's mouse Ensembl BED12, which stores transcript ids in the gene_id field.

The reformatted BED carries a transcript id where a gene id belongs, so grouping
transcripts by gene gives one transcript per "gene". This joins the BED against Ensembl
transcript metadata on the transcript id and writes the gene id and both biotypes back into
the trailing columns.

UCSC's mm39 Ensembl track is a snapshot of Ensembl v110, so the transcript table must be
fetched at that version or the join silently drops transcripts. The join is reported as a
rate so a version mismatch is visible rather than quiet.

Usage:
    FixMouseEnsemblBed.py --bed IN.bed --transcript-info INFO.tsv.gz --out OUT.bed
"""
import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--transcript-info", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--species", default="Mouse")
    ap.add_argument("--min-join-rate", type=float, default=0.90)
    args = ap.parse_args()

    # BED12 plus trailing annotation columns; everything read as text so ids stay intact
    bed = pd.read_csv(args.bed, sep="\t", header=None, dtype=str)
    bed = bed.iloc[:, :16]
    bed.columns = [f"X{i}" for i in range(1, bed.shape[1] + 1)]

    info = pd.read_csv(args.transcript_info, sep="\t", dtype=str)
    info = info[info.Species == args.species][
        ["ensembl_transcript_id_version", "ensembl_gene_id_version",
         "transcript_biotype", "gene_biotype"]].drop_duplicates(
        "ensembl_transcript_id_version")

    merged = bed.merge(info, left_on="X13",
                       right_on="ensembl_transcript_id_version", how="inner")
    rate = len(merged) / len(bed) if len(bed) else 0.0
    print(f"joined {len(merged)}/{len(bed)} BED rows ({rate:.1%}) on transcript id",
          file=sys.stderr)
    if rate < args.min_join_rate:
        sys.exit(f"ERROR: join rate {rate:.1%} below {args.min_join_rate:.0%}. The "
                 f"transcript table is probably from a different Ensembl release than the "
                 f"UCSC track (expected v110 for mm39).")

    cols = [f"X{i}" for i in range(1, 13)] + [
        "X13", "ensembl_gene_id_version", "transcript_biotype", "gene_biotype"]
    merged[cols].to_csv(args.out, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
