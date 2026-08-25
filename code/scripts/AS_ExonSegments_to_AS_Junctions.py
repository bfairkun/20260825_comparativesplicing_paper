#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : AS_ExonSegments_to_AS_Junctions.py
# @created     : Friday Jun 21, 2024 12:26:19 CDT
#
# @description : 
######################################################################

import sys
import bedparse
import argparse
import gzip
import re
import pysam
import logging
from io import StringIO
import pybedtools
import pandas as pd
from collections import Counter


def get_NMD_detective_B_classification(sequence):
    """
    sequence should be marked with '^' for start, '*' for stop, and '|' for splice juncs
    """
    CDS = re.search(r"\^(\w+)\*", sequence.replace("|", ""))
    InternalStopExon = re.search(r"(^|\|)([\^ACGTNacgtn]*\*[ACGTNacgtn]*)\|", sequence)
    if "^" not in sequence or CDS == None:
        return "No CDS"
    elif re.search(r"\^[\w|]+$", sequence):
        return "No stop"
    elif re.search(r"(^|\|)[\^ACGTNacgtn]*\*[ACGTNacgtn]*$", sequence):
        return "Last exon"
    elif len(CDS.group(1)) <= 125:
        return "Start proximal"
    elif len(InternalStopExon.group(2)) >= 407:
        return "Long exon"
    elif re.search(r"\*[ACGTNacgtn]{0,50}\|[ACGTNacgtn]+$", sequence):
        return "50 nt rule"
    else:
        return "Trigger NMD"


def count_bars_until_n_position(sequence, n):
    bar_count = 0
    non_bar_count = 0
    for i,char in enumerate(sequence):
        if non_bar_count == n:
            break
        if char == "|":
            bar_count += 1
        else:
            non_bar_count += 1
    return bar_count


def AddORF_Marks(sequence, StartMarkBasePosition=None, StopMarkBasePosition=None):
    """
    returns sequence with '^' at StartmarkBasePosition and '*' at StopMarkBasePosition. '|' strings are ignored in base position. 
    """
    StartMark = '^'
    StopMark = '*'
    if StartMarkBasePosition == None:
        StartMark = ''
        StartMarkBasePosition = 0
    if StopMarkBasePosition == None:
        StopMark = ''
        StopMarkBasePosition = 0
    OffsetToStartPosition = count_bars_until_n_position(sequence, StartMarkBasePosition)
    OffsetToStopPosition = count_bars_until_n_position(sequence, StopMarkBasePosition)
    return sequence[0:StartMarkBasePosition + OffsetToStartPosition] + StartMark + sequence[(StartMarkBasePosition + OffsetToStartPosition):(StopMarkBasePosition+OffsetToStopPosition)] + StopMark + sequence[(StopMarkBasePosition+OffsetToStopPosition):]
AddORF_Marks("||ATG|G|ATAGG", 1, 5)

def insert_marks_for_first_ORF(sequence, require_STOP=True, min_ORF_len=0):
    """
    return sequence with "^" to mark start codon, "*" to mark stop for longest ORF. If not require_STOP, translation does not need stop codon at end. "|" characters (which I use mark splice junctions) are ignored, allowing codons to cross exon junctions. min_ORF_len in codons after start codon.
    """
    if require_STOP:
        regex = r"^.*?(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN]){" + str(min_ORF_len) + r",}((?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A))).*$"
    else:
        regex = r"^.*?(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN]){" + str(min_ORF_len) + r",}((?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A))?).*$"
    first_ORF_match = re.match(regex, sequence, flags=re.IGNORECASE)
    if first_ORF_match:
        start_codon_pos = first_ORF_match.start(1)
        if first_ORF_match.group(2):
            stop_codon_pos = first_ORF_match.end(1)
            return sequence[0:start_codon_pos] + "^" + first_ORF_match.group(1) + "*" + sequence[stop_codon_pos:]
        else:
            return sequence[0:start_codon_pos] + "^" + sequence[start_codon_pos:]
    else:
        return sequence

def extract_sequence(self, fasta_obj, AddMarksForORF=False):
    """
    ...method to be used on bedparse.bedline object. will monkey patch it into bedparse.bedline class
    will also add '|' between blocks to mark exon junctions
    """
    # Extract information from the BED entry
    chrom = self.chr
    start = self.start
    end = self.end
    strand = self.strand if hasattr(self, 'strand') else '+'
    block_sizes = list(map(int, self.exLengths.split(','))) if hasattr(self, 'exLengths') else [self.end - self.start]
    block_starts = list(map(int, self.exStarts.split(','))) if hasattr(self, 'exStarts') else [0]

    sequence = ''
    # Retrieve the sequence for each block and concatenate them
    if block_sizes and block_starts:
        for block_start, block_size in zip(block_starts, block_sizes):
            block_start_genomic = start + block_start
            block_end_genomic = block_start_genomic + block_size
            block_seq = fasta_obj.fetch(chrom, block_start_genomic, block_end_genomic)
            sequence += "|" + block_seq
        sequence = sequence.rstrip('|').lstrip('|')
    else:
        sequence = fasta_obj.fetch(chrom, start, end)
    # Reverse complement if on the negative strand
    if strand == '-':
        sequence = reverse_complement(sequence)
        # Add ORF markers if requested
    if AddMarksForORF and self.cds():
        if self.utr(which=5):
            cds_relative_start = sum([int(i) for i in self.utr(which=5).exLengths.split(',')])
        else:
            cds_relative_start = 0
        cds_relative_end = cds_relative_start + sum([int(i) for i in self.cds().exLengths.split(',')])
        sequence = AddORF_Marks(sequence, cds_relative_start, cds_relative_end)
    return sequence

def reverse_complement(seq):
    # Helper function to reverse complement a DNA sequence
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]

# Monkey-patch the BEDLine class
bedparse.bedline.extract_sequence = extract_sequence
bedparse.bedline.reverse_complement = staticmethod(reverse_complement)

def get_absolute_pos(bedline, coord):
    """
    Similar to bedline.tx2genome(coord, stranded=True), but I think that function is buggy and doesn't work for getting last base in tx
    """
    cumulative_block_size = 0
    exLengths = [int(i) for i in bedline.exLengths.split(',')]
    exStarts = [int(i) for i in bedline.exStarts.split(',')]
    if bedline.strand == '+':
        for i, block_size in enumerate(exLengths):
            remainder = coord - cumulative_block_size
            cumulative_block_size += block_size
            if cumulative_block_size >= coord:
                absolute_pos = bedline.start + remainder + exStarts[i]
                break
    elif bedline.strand == '-':
        for i, block_size in enumerate(list(reversed(exLengths))):
            remainder = coord - cumulative_block_size
            cumulative_block_size += block_size
            if cumulative_block_size >= coord:
                absolute_pos = bedline.start + list(reversed(exStarts))[i] + block_size - remainder
                break
    return absolute_pos

def adjust_start(interval):
    """Decrement the start position by 1."""
    interval.start = int(interval.start) - 1
    return interval

def convert_bed6_to_bed12(bed6_list, name=None, score=0, thickStart=None, thickEnd=None, color="0"):
    bed6_list = sorted(bed6_list, key=lambda x: x.start)
    chr = bed6_list[0].chr
    start = min(bed6_list, key=lambda x: x.start).start
    end = max(bed6_list, key=lambda x: x.end).end
    strand = bed6_list[0].strand
    blockCount = len(bed6_list)
    blockSizes = ",".join([str(exon.end - exon.start) for exon in bed6_list])
    blockStarts = ",".join([str(exon.start - start) for exon in bed6_list])
    
    if thickStart is None:
        thickStart = start
    if thickEnd is None:
        thickEnd = end
    if name is None:
        name = "isoform"

    bed12_fields = [
        chr, start, end, name, score, strand, thickStart, thickEnd, color, blockCount, blockSizes, blockStarts
    ]
    return bedparse.bedline(bed12_fields)

def Get_bedline_coords(bedline):
    return bedline.chr, bedline.start, bedline.end, bedline.strand

def process_ce_segment(AS_segment, transcript):
    exons = transcript.bed12tobed6()
    if Get_bedline_coords(AS_segment) in [Get_bedline_coords(exon) for exon in exons]: 
        long_isoform_exons = [exon for exon in exons]
        short_isoform_exons = [exon for exon in exons if not Get_bedline_coords(AS_segment) == Get_bedline_coords(exon)]
    else:
        long_isoform_exons = [exon for exon in exons] + [AS_segment]
        short_isoform_exons = [exon for exon in exons]
    long_isoform = convert_bed6_to_bed12(long_isoform_exons)
    short_isoform = convert_bed6_to_bed12(short_isoform_exons)
    return long_isoform, short_isoform

def process_aa_ad_segment(AS_segment, transcript):
    exons = transcript.bed12tobed6()
    short_isoform_exons = []
    long_isoform_exons = []
    for exon in exons:
        if exon.start == AS_segment.start:
            if AS_segment.end < exon.end:
                short_isoform_exons.append(bedparse.bedline([exon.chr, AS_segment.end, exon.end, ".", ".", exon.strand]))
                long_isoform_exons.append(bedparse.bedline([exon.chr, exon.start, exon.end, ".", ".", exon.strand]))
            else:
                return None, None
        elif exon.end == AS_segment.end:
            if AS_segment.start > exon.start:
                short_isoform_exons.append(bedparse.bedline([exon.chr, exon.start, AS_segment.start, ".", ".", exon.strand]))
                long_isoform_exons.append(bedparse.bedline([exon.chr, exon.start, exon.end, ".", ".", exon.strand]))
            else:
                return None, None
        elif exon.start == AS_segment.end:
            short_isoform_exons.append(bedparse.bedline([exon.chr, exon.start, exon.end, ".", ".", exon.strand]))
            long_isoform_exons.append(bedparse.bedline([exon.chr, AS_segment.start, exon.end, ".", ".", exon.strand]))
        elif exon.end == AS_segment.start:
            short_isoform_exons.append(bedparse.bedline([exon.chr, exon.start, exon.end, ".", ".", exon.strand]))
            long_isoform_exons.append(bedparse.bedline([exon.chr, exon.start, AS_segment.end, ".", ".", exon.strand]))
        else:
            short_isoform_exons.append(exon)
            long_isoform_exons.append(exon)
    long_isoform = convert_bed6_to_bed12(long_isoform_exons)
    short_isoform = convert_bed6_to_bed12(short_isoform_exons)
    if [Get_bedline_coords(exon) for exon in short_isoform_exons] == [Get_bedline_coords(exon) for exon in long_isoform_exons]:
        return None, None
    else:
        return long_isoform, short_isoform

def bedlines_are_equal(bedline1, bedline2):
    # Check basic attributes: chr, start, end, strand
    if (bedline1.chr != bedline2.chr or 
        bedline1.start != bedline2.start or 
        bedline1.end != bedline2.end or 
        bedline1.strand != bedline2.strand):
        return False

    # Check block attributes if both bedlines have blocks
    if bedline1.bedType >= 6 and bedline2.bedType >= 6:
        if bedline1.nEx != bedline2.nEx:
            return False
        
        for exon1, exon2 in zip(bedline1.bed12tobed6(), bedline2.bed12tobed6()):
            if exon1.start != exon2.start or exon1.end != exon2.end:
                return False
    return True

def GetLongIsoformBlockForAS_Event(long_isoform, short_isoform):
    long_isoform_exons = long_isoform.bed12tobed6()
    short_isoform_exons = short_isoform.bed12tobed6()
    while len(long_isoform_exons) > len(short_isoform_exons):
        short_isoform_exons.append(short_isoform_exons[0])
    for i, (long_isoform_exon, short_isoform_exon) in enumerate(zip(long_isoform_exons, short_isoform_exons)):
        if Get_bedline_coords(long_isoform_exon) != Get_bedline_coords(short_isoform_exon):
            return i

def insert_marks_for_defined_ORF(sequence, start_codon_pos=None):
    """
    return sequence with "^" to mark start codon at provided string position, "*" to mark stop for longest ORF. If no stop codon is found, then no * will be inserted.
    """
    if start_codon_pos != None and start_codon_pos < len(sequence):
        # orf_match = re.search(r"^\|?(A\|?T\|?G\|?(?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])*?)\|?(?:T\|?A\|?A|T\|?A\|?G|T\|?G\|?A)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        # orf_match_nostop = re.search(r"^\|?(A\|?T\|?G\|?(?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])+)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        orf_match = re.search(r"^\|?((?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])*?)\|?(?:T\|?A\|?A|T\|?A\|?G|T\|?G\|?A)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        orf_match_nostop = re.search(r"^\|?((?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])+)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        if orf_match:
            # orf match with start and stop codon
            return sequence[:start_codon_pos] + "^" + orf_match.group(0) + "*" + sequence[orf_match.span(0)[1] + start_codon_pos:]
        elif orf_match_nostop:
            # orf match with no stop
            return sequence[:start_codon_pos] + "^" + sequence[start_codon_pos:]
        else:
            return sequence
    else:
        return sequence

def translate_from_cds_startCodon(bedline, cdsStart, cdsEnd, fasta_obj):
    bedline.cdsStart, bedline.cdsEnd = cdsStart, cdsEnd
    try:
        sequence = bedline.extract_sequence(fasta_obj, AddMarksForORF=True).replace("*", "")
    except bedparse.BEDexception:
        sequence = bedline.extract_sequence(fasta_obj, AddMarksForORF=False).replace("*", "")
    if "^|" not in sequence and sequence.find("^") != -1:
        remarked_sequence = insert_marks_for_defined_ORF(sequence.replace("^", ""), start_codon_pos=sequence.find("^"))
    else:
        remarked_sequence = insert_marks_for_first_ORF(sequence.replace("*", "").replace("^", ""), require_STOP=False, min_ORF_len=42)
    return remarked_sequence

def parse_args(args=None):
    parser = argparse.ArgumentParser(description='Process bed12+ (mRNA_end_NF and other relevant tags in col 22, transcript type in col 18), and AS_segment bed9 (AS segments like from Mazin et al) where col9 color is coded as "CE"="#1b9e77", "RI"="#d95f02", "AA"="#7570b3", "AD"="#e7298a". Output tsv file describing affect of inclusion/exclusion on protein_coding transcripts for each transcript-intersecting AS_segment. Will translate transcript with out without the AS_segment, and determine NMD status with NMDFinderB algorithm. Also output leafcutter_perind.counts.gz file containing just junc names for processing the relevant junctions with Yangs SpliceJunctionClassifier.py script.')
    parser.add_argument('-AnnotatedTranscriptsIn', dest='transcripts_fn', required=True, help='bed12+ file for each transcript. Only protein_coding transcripts (with "protein_coding" as value for col18 and no "_NF"-containing string in col22) will be intersected with AS_segments')
    parser.add_argument('-AS_segmentsIn', dest='ASsegments_fn', required=True, help='bed9 for each AS segment, which is either CE (cassette exon) with #1b9e77 in color column, RI (retained intron) with #d95f02 in color column, alternative donor (AD, #e7298a), or alternative acceptor (AA, #7570b3)')
    parser.add_argument('-fa', dest='fasta_in', required=True, help='Input FASTA file.')
    parser.add_argument('-TsvOut', dest='tsv_out', required=True, help='tsv file out')
    parser.add_argument('-JuncListOut', dest='junclist_out', required=True, help='leafcutter_perind.counts.gz file containing just junc names for processing the relevant junctions with Yangs SpliceJunctionClassifier.py script.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase verbosity')
    return parser.parse_args(args)


def setup_logging(verbose):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

def main(args=None):
    args = parse_args(args)
    setup_logging(args.verbose)

    logging.debug(f'{args}')
    logging.info("Reading in files...")
    # Input files
    fasta_obj = pysam.FastaFile(args.fasta_in)

    # Read the input BED files
    # AS type coded as color column:
    # "CE"="#1b9e77", "RI"="#d95f02", "AA"="#7570b3", "AD"="#e7298a"
    adjusted_orth_segs = pybedtools.BedTool(args.ASsegments_fn).saveas()
    logging.debug(f'read in {len(adjusted_orth_segs)} AS segments')

    # Get protein_coding transcripts
    protein_coding_bed = pybedtools.BedTool(args.transcripts_fn).filter(lambda x: x.fields[16] == "protein_coding" and "NF" not in x.fields[20]).each(lambda x: x[0:12]).saveas()
    logging.debug(f'read in {len(protein_coding_bed)} full length protein coding transcripts...')

    # convert protein_coding transcripts to dictionary of bedparse.bedline objects, becasue some weird bug in bedtools that modifies the transcript blocks after intersection
    protein_coding_bedline_dict = dict()
    for entry in protein_coding_bed:
        protein_coding_bedline_dict[entry[3]] = bedparse.bedline(entry[0:12])

    # Convert BED12 to BED6 to get protein_coding transcript exons
    protein_coding_exons = protein_coding_bed.bed12tobed6().sort()

    # Intersect the exons with the adjusted orth_segs_bed
    intersected = protein_coding_exons.intersect(adjusted_orth_segs, s=True, wo=True)

    # Filter the intersected results to remove lines where the start and end positions are the same. I want to exclude these later
    SegmentsToExclude1 = intersected.filter(lambda x: not (x.fields[1] == x.fields[7] and x.fields[2] == x.fields[8]))

    # I also want to exclude AS_segments where that totally enclose an annotated exon
    SegmentsToExclude2 = protein_coding_exons.intersect(adjusted_orth_segs, s=True, wo=True, f=1).filter(lambda x: (int(x.fields[8]) - int(x.fields[7])) > int(x.fields[15]))
    # SegmentsToExclude2_temp = protein_coding_exons.intersect(adjusted_orth_segs, s=True, wo=True, f=1).filter(lambda x: (int(x.fields[8]) - int(x.fields[7])) > int(x.fields[15])).saveas()
    # for segment in SegmentsToExclude2_temp:
    #     if segment.fields[3] == "chi.868.s2":
    #         print(segment)

    SegmentsToExcludeUnion = set([i[9] for i in SegmentsToExclude1] + [i[9] for i in SegmentsToExclude2])
    logging.info(f"Excluding {len(SegmentsToExcludeUnion)} segments that partially overlap exons...")

    Filtered_AS_Segments = adjusted_orth_segs.filter(lambda x: (x[3] not in SegmentsToExcludeUnion) or  (x[8] not in ["#7570b3", "#e7298a", "#1b9e77"])).saveas()
    logging.debug(f"{len(Filtered_AS_Segments)} Filtered AS segments will be intersected with transcripts...")
    logging.debug(f'types of AS_segments to intersect: {Counter([l[8] for l in Filtered_AS_Segments])}')

    logging.debug("Intersecting segments with protein_coding transcripts...")
    # Intersect exons with protein_coding_transcripts
    SegmentsOverlapTranscripts = Filtered_AS_Segments.intersect(protein_coding_bed, wo=True, sorted=True, f=1, s=True)

    ### loop over intersections and determine NMD status of each isoform with our without the AS_segment
    AS_segments_coords_set = set([(l[0], int(l[1]), int(l[2]), l[5]) for l in adjusted_orth_segs])

    logging.info(f"Iterating through {len(SegmentsOverlapTranscripts)} intersections between AS_segments and protein_coding transcripts and determining effects with or without AS segment included...")
    SuccessfullyAssembledIsoformPairs = []
    rows = []
    UniqJuncsSet = set()
    for i, l in enumerate(SegmentsOverlapTranscripts):
        try:
            if i % 1000 == 0: logging.debug(f'processed {i} intersections...')
            transcript = protein_coding_bedline_dict[l[12]]
            AS_segment = bedparse.bedline(l[0:9] + ["1", str(int(l[2])-int(l[1])), "0"])
            row = {"transcript": transcript.name, "AS_segment": AS_segment.name, "color":AS_segment.color, "strand":transcript.strand} 
            # if AS_segment.name == "hum.55945.s6" and transcript.name: break
            # if i ==39: break
            if AS_segment.color == "#1b9e77":  # CE
                long_isoform, short_isoform = process_ce_segment(AS_segment, transcript)
            elif AS_segment.color in ["#7570b3", "#e7298a"]:  # AA or AD
                long_isoform, short_isoform = process_aa_ad_segment(AS_segment, transcript)
            else:
                long_isoform = None
                short_isoform = None
            if long_isoform and short_isoform and not bedlines_are_equal(long_isoform, short_isoform) and short_isoform.introns():
                SuccessfullyAssembledIsoformPairs.append(AS_segment.color)
                # print(bedlines_are_equal(long_isoform, transcript), bedlines_are_equal(short_isoform, transcript))
                WhichIsoformIsAnnotated = "LongIsoform" if bedlines_are_equal(long_isoform, transcript) else "ShortIsoform"
                long_junctions = set([Get_bedline_coords(junc) for junc in long_isoform.introns().bed12tobed6()])
                short_junctions = set([Get_bedline_coords(junc) for junc in short_isoform.introns().bed12tobed6()])
                unique_long_junctions = ','.join([f'{chr}_{start}_{stop}_{strand}' for chr, start, stop, strand in (long_junctions - short_junctions)])
                unique_short_junctions = ','.join([f'{chr}_{start}_{stop}_{strand}' for chr, start, stop, strand in (short_junctions - long_junctions)])
                UniqJuncsSet.update(long_junctions ^ short_junctions)
                long_isoform_exons = long_isoform.bed12tobed6()
                Notes = ""
                NMDFinderB_long = get_NMD_detective_B_classification(translate_from_cds_startCodon(long_isoform, transcript.cdsStart, transcript.cdsEnd, fasta_obj))
                NMDFinderB_short = get_NMD_detective_B_classification(translate_from_cds_startCodon(short_isoform, transcript.cdsStart, transcript.cdsEnd, fasta_obj))
                ExIndexForAS = GetLongIsoformBlockForAS_Event(long_isoform, short_isoform)
                # print(ExIndexForAS, len(long_isoform_exons))
                if ExIndexForAS == 0 or ExIndexForAS == len(long_isoform_exons) - 1:
                    Notes = "AS at terminal exon"
                elif Get_bedline_coords(long_isoform_exons[ExIndexForAS-1]) in AS_segments_coords_set:
                    Notes = "Left ex is also AS segment"
                elif Get_bedline_coords(long_isoform_exons[ExIndexForAS+1]) in AS_segments_coords_set:
                    Notes = "Right ex is also AS segment"
                row.update({"LongIsoform_UniqueJuncs":unique_long_junctions, "ShortIsoform_UniqueJuncs":unique_short_junctions, "LongIsoform_NMDFinderB":NMDFinderB_long, "ShortIsoform_NMDFinderB":NMDFinderB_short, "Notes":Notes, "WhichIsoformIsAnnotated":WhichIsoformIsAnnotated})
                # if WhichIsoformIsAnnotated=="LongIsoform" and NMDFinderB_long=="No CDS": break
        except bedparse.BEDexception: pass
        except:
            breakpoint()
        rows.append(row)
    logging.debug(Counter(SuccessfullyAssembledIsoformPairs))
    Results = pd.DataFrame(rows)

    logging.info("Writing output...")
    Results.to_csv(args.tsv_out, sep='\t', index=False)
    if args.junclist_out.endswith('.gz'):
        with gzip.open(args.junclist_out, 'wt') as fout:  # 'wt' mode opens the file in text mode for writing
            _ = fout.write("chrom\n")
            for chr,start,stop,strand in UniqJuncsSet:
                _ = fout.write(f'{chr}:{start}:{stop+1}:clu_n_{strand}\n')
    else:
        with open(args.junclist_out, 'w') as fout:
            _ = fout.write("chrom\n")
            for chr,start,stop,strand in UniqJuncsSet:
                _ = fout.write(f'{chr}:{start}:{stop+1}:clu_n_{strand}\n')


if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        main("-TsvOut scratch/test_ASSegmentResults.tsv.gz -JuncListOut scratch/test_ASSegmentResults.juncs.gz -AS_segmentsIn kaessman_AS_dat/AS_segment_lists/Human_ensemblv75.bed -AnnotatedTranscriptsIn GenomeFiles/Human_ensemblv75/Reannotated.B.bed.gz -fa /project2/yangili1/bjf79/ReferenceGenomes/Human_ensemblv75/Reference.fa -v".split(' '))
        main("-TsvOut scratch/test.Chicken_ensemblv84.ASSegmentResults.tsv.gz -JuncListOut scratch/test.Chicken_ensemblv84.ASSegmentResults.juncs.gz -AS_segmentsIn scratch/AS_segments.Chicken_ensemblv84.bed -AnnotatedTranscriptsIn GenomeFiles/Chicken_ensemblv84/Reannotated.B.bed.gz -fa /project2/yangili1/bjf79/ReferenceGenomes/Chicken_ensemblv84/Reference.fa -v".split(' '))
    else:
        main()
