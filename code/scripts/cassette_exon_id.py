"""
This script takes leafcutter 2 input and GTF and outputs potential cassette exons. 
The script tries it match 5' and 3' end of non_coding + UTR False introns from leafcutter 2 with the 3' and 5' end of coding exons (CDS) from the gtf.
If both those conditons are met and there is a [3'-5'] distance of <500 between 2 non-coding introns then the intervening region is called as a 
cassette exon
Author: Dylan Stermer 
Date: January 30, 2025
"""

import argparse
import gzip
import sys
import logging

def setup_logging(verbose):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Identify potential cassette exons by analyzing Leafcutter 2 introns and GTF coding regions."
    )
    parser.add_argument("--gtf", required=True, help="Path to the reference GTF file.")
    parser.add_argument("--leafcutter2", required=True, help="Path to the Leafcutter 2 junction classification file. Coordinates should be 0-start, half-open, like bed file. If not, may need to change offset parameter")
    parser.add_argument("--output", required=True, help="Output BED file for detected cassette exons. Each feature is cassette exon")
    parser.add_argument("--bed12", required=False, help="Output BED12 file for detected cassette exons. Each feature is cassette exon and it's flanking introns.")
    parser.add_argument("--offset", type=int, choices=[0, 1], help="Offset for coordinate conversion between junctions and GTF.", default=0)
    parser.add_argument("--FlankingIntronLeaf2Flags", type=int, nargs='+', default=[0, 4, 5, 8, 12], help="Which leaf2 flags to consider for flanking introns to identify exons. Modify these to instead look for productive cassette exons (that are flanked by 'productive' splice junctions), vs poison exons (that are flanked by 'non-productive' splice junctions). The leaf2 flags are a binary to decimal conversion of the T/F (1/0) flags in the leaf2 file: UTR (first, rightmost digit), coding, annotation, GencodePC (last, leftmost digit). For example, UTR=False, coding=False, annotation=True, GenocdePC=False would be 0100 which is 4")
    parser.add_argument("--SkippingIntronLeaf2Flags", type=int, nargs='+', default=[12, 13, 14, 15], help="Which leaf2 flags to consider for the exon-skipping junction to identify exons. Similar to FlankingIntronLeaf2Flags, but for the exon-skipping junction.")
    parser.add_argument("--MaxExonSize", type=int, default=500, help="Maximum size of cassette exon to consider")
    parser.add_argument("--LongOutputNames", action='store_true', default=False, help="Long, more descriptive, name for bed column4")
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase verbosity')
    return parser.parse_args(args)

def parse_gtf(gtf_path):
    dict_cds = {}
    with open(gtf_path) as gtf_file:
        for line in gtf_file:
            if line.startswith("#"): 
                continue
            chr, source, ann_type, start, stop, dot, strand, frame, info = line.strip().split("\t")
            if ann_type != "CDS":
                continue
            gene_name = info.split('gene_name "')[1].split('"')[0]
            dict_cds.setdefault(gene_name, set()).add((chr, int(start), int(stop), strand, frame))
    return dict_cds

def BooleanListToDecimalFlag(bool_list):
    return sum([2**i for i, x in enumerate(bool_list) if x])

def parse_leafcutter2(leafcutter_path, valid_decimal_flags=[0,4,8,12]):
    dict_leafcutter2 = {}
    with open(leafcutter_path) as leafcutter_file:
        for line in leafcutter_file:
            if line.strip().startswith("Gene_name"):
                continue
            gene_name, intron_cord, strand, annotation, coding, UTR, GencodePC = line.strip().split("\t")
            FlagsAsList = [item == "True" for item in (UTR, coding, annotation, GencodePC)]
            DecimalFlag = BooleanListToDecimalFlag(FlagsAsList)
            chrom, positions = intron_cord.split(":")
            start, stop = positions.split("-")
            if DecimalFlag in valid_decimal_flags:
                dict_leafcutter2.setdefault(gene_name, set()).add((chrom, int(start), int(stop), DecimalFlag))
    return dict_leafcutter2

def identify_cassette_exon_junctions(dict_cds, dict_leafcutter2, offset):
    dict_cassette_5_junction = {}
    dict_cassette_3_junction = {}
    
    for gene_name in dict_leafcutter2:
        if gene_name not in dict_cds:
            continue
        
        coding_exons = dict_cds[gene_name]
        non_code_introns = dict_leafcutter2[gene_name]
        
        for non_code_chr, non_code_start, non_code_stop, non_code_annote in non_code_introns:
            for coding_chr, coding_start, coding_stop, coding_strand, frame in coding_exons:
                if non_code_chr != coding_chr:
                    continue
                if non_code_start == coding_stop:
                    dict_cassette_5_junction.setdefault(gene_name, []).append(
                        (non_code_stop, non_code_start, non_code_chr, f"{coding_chr}:{coding_start}-{coding_stop}", coding_strand, frame)
                    )
                if non_code_stop == coding_start - offset:
                    dict_cassette_3_junction.setdefault(gene_name, []).append(
                        (non_code_start, non_code_stop, non_code_chr, f"{coding_chr}:{coding_start}-{coding_stop}", coding_strand, frame)
                    )
    return dict_cassette_5_junction, dict_cassette_3_junction

def detect_cassette_exons(dict_cassette_5_junction, dict_cassette_3_junction, dict_leafcutter2_potential_skipping_juncs, MaxExonSize=500):
    dict_cassette_exons = {}
    for gene_name in dict_cassette_5_junction:
        if gene_name not in dict_cassette_3_junction or gene_name not in dict_leafcutter2_potential_skipping_juncs:
            continue
        
        non_coding_upstream_intron = sorted(dict_cassette_5_junction[gene_name])
        non_coding_downstream_intron = sorted(dict_cassette_3_junction[gene_name])
        for up_end, up_start, up_chr, up_coding_exon, coding_strand, up_exon_frame in non_coding_upstream_intron:
            for down_start, down_end, down_chr, down_coding_exon, coding_strand, down_exon_frame in non_coding_downstream_intron:
                if up_chr != down_chr:
                    continue
                if int(up_end+1) == int(down_end-1) or int(up_start) == int(down_start) or int(up_end) > int(down_start):
                    continue
                if abs(up_end - down_start) <= MaxExonSize and any((up_chr, up_start, down_end) == (j[0], j[1], j[2]) for j in dict_leafcutter2_potential_skipping_juncs[gene_name]):
                    dict_cassette_exons.setdefault(gene_name, []).append(
                        (up_end, down_start, abs(up_end - down_start), up_chr,
                         f"{up_chr}:{up_start}-{up_end+1}", f"{down_chr}:{down_start}-{down_end-1}", up_coding_exon, up_exon_frame, down_coding_exon, down_exon_frame, coding_strand)
                    )
    return dict_cassette_exons

def write_bed(cassette_exons, output_file, write_as_blocked_bed=True, LongNames=False):
    """
    Write cassette exon regions to a BED file.

    If write_as_blocked_bed is True, the output will use BED12 format where:
    - The upstream and downstream **exons** are blocks.
    - The cassette exon is a block.
    - Thin intron lines will connect the three exon blocks.

    If write_as_blocked_bed is False, it will write a standard BED6 file.
    """
    with open(output_file, "w") as bed:
        bed.write("#track name='CassetteExons' description='Predicted Cassette Exons' useScore=0\n")
        if LongNames:
            name_header = "gene_name;upstream_exon;upstream_frame;upstream_intron;cassette_exon;downstream_intron;downstream_exon;downstream_frame"
        else:
            name_header = "gene_name"
        if write_as_blocked_bed:
            bed.write(f"#chrom\tstart\tend\t{name_header}\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\n")
        else:
            bed.write(f"#chrom\tstart\tend\t{name_header}\tscore\tstrand\n")        
        for gene_name, exons in cassette_exons.items():
            for up_end, down_start, size, chrom, upstream_intron, downstream_intron, up_coding_exon, up_exon_frame, down_coding_exon, down_exon_frame, strand in exons:
                # Parse upstream and downstream intron coordinates
                up_chr, up_range = upstream_intron.split(":")
                up_start, up_end = map(int, up_range.split("-"))

                down_chr, down_range = downstream_intron.split(":")
                down_start, down_end = map(int, down_range.split("-"))

                # Parse upstream and downstream exon coordinates
                up_chr, up_exon_range = up_coding_exon.split(":")
                up_exon_start, up_exon_end = map(int, up_exon_range.split("-"))
                up_exon_start = up_exon_start - 1

                down_chr, down_exon_range = down_coding_exon.split(":")
                down_exon_start, down_exon_end = map(int, down_exon_range.split("-"))
                down_exon_start = down_exon_start - 1

                # BED start and end should span the entire feature (flanking exons + introns)
                bed_start = up_exon_start  # Start at the beginning of the upstream exon
                bed_end = down_exon_end    # End at the end of the downstream exon

                # Exon block coordinates (cassette exon)
                exon_start = up_end - 1
                exon_end = down_start

                # Compute block sizes
                up_exon_size = up_exon_end - up_exon_start
                cassette_exon_size = exon_end - exon_start
                down_exon_size = down_exon_end - down_exon_start

                block_sizes = f"{up_exon_size},{cassette_exon_size},{down_exon_size}"

                # Compute block starts relative to bed_start
                block_starts = f"0,{exon_start - bed_start},{down_exon_start - bed_start}"

                # Write BED file
                if LongNames:
                    name = f"{gene_name};{up_coding_exon};{up_exon_frame};{up_chr}:{up_start+1}-{up_end-1};{chrom}:{exon_start+1}-{exon_end};{down_chr}:{down_start+1}-{down_end+1};{down_coding_exon};{down_exon_frame}"
                else:
                    name = f"{gene_name}"
                if write_as_blocked_bed:
                    bed.write(
                        f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t0\t{strand}\t"
                        f"{bed_start}\t{bed_end}\t0\t3\t{block_sizes}\t{block_starts}\n"
                    )
                else:
                    # Standard BED6 format (only the cassette exon)
                    bed.write(
                        f"{chrom}\t{exon_start}\t{exon_end}\t{name}\t0\t{strand}\n"
                    )

def main(args=None):
    args = parse_args(args)
    setup_logging(args.verbose)
    logging.debug("Parsing GTF file...")
    dict_cds = parse_gtf(args.gtf)
    logging.debug("Reading leafcutter 2 annotation file...")
    dict_leafcutter2_potential_flanking_juncs = parse_leafcutter2(args.leafcutter2, valid_decimal_flags=args.FlankingIntronLeaf2Flags)
    dict_leafcutter2_potential_skipping_juncs = parse_leafcutter2(args.leafcutter2, valid_decimal_flags=args.SkippingIntronLeaf2Flags)
    logging.debug("Identifying cassette exons...")
    dict_cassette_5_junction, dict_cassette_3_junction = identify_cassette_exon_junctions(dict_cds, dict_leafcutter2_potential_flanking_juncs, args.offset)
    dict_cassette_exons = detect_cassette_exons(dict_cassette_5_junction, dict_cassette_3_junction, dict_leafcutter2_potential_skipping_juncs, MaxExonSize=args.MaxExonSize)
    # breakpoint()
    logging.debug("Writing output...")
    write_bed(dict_cassette_exons, args.output, write_as_blocked_bed=False, LongNames=args.LongOutputNames)
    if args.bed12:
        write_bed(dict_cassette_exons, args.bed12, write_as_blocked_bed=True, LongNames=args.LongOutputNames)

if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        # some tests
        # Find productive cassette exons
        main("--gtf MazinLeafcutterAnalysis/ReformatedGTFs/Human_UCSC.hg38_GencodeComprehensive46.gtf --leafcutter2 MazinLeafcutterAnalysis/ClassifyJuncs/Human_UCSC.hg38_GencodeComprehensive46.Clustered._junction_classifications.txt --output scratch/productive_ints.bed --bed12 scratch/productive_CE.withflanks.bed -v --offset 1 --FlankingIntronLeaf2Flags 14 15".split(' '))
        # Find poison exons
        main("--gtf MazinLeafcutterAnalysis/ReformatedGTFs/Human_UCSC.hg38_GencodeComprehensive46.gtf --leafcutter2 MazinLeafcutterAnalysis/ClassifyJuncs/Human_UCSC.hg38_GencodeComprehensive46.Clustered._junction_classifications.txt --output scratch/poison_id.bed --bed12 scratch/poison_id.withflanks.bed -v --offset 1 --LongOutputNames".split(' '))
    main()