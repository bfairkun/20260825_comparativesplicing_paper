#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : bedparse_translate_transcripts
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
import subprocess
from io import StringIO
from collections import defaultdict
import shutil
import pandas as pd
import csv
import logging

def reorder_gtf(gtf_stringio, output_gtf):
    """Reorder GTF lines by gene and transcript, ensuring correct hierarchical sorting,
    and write to output file."""
    # Read the GTF data into a pandas DataFrame
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gtf_stringio, sep='\t', names=columns, comment='#')

    # Extract gene_id and transcript_id from the attribute column
    df['gene_id'] = df['attribute'].str.extract(r'gene_id "([^"]+)"')
    df['transcript_id'] = df['attribute'].str.extract(r'transcript_id "([^"]+)"')

    # Determine if the line is a gene line
    df['is_gene'] = df['transcript_id'].isnull()

    # Create a column for the start position of the parent gene
    df['start_parent_gene'] = df.groupby('gene_id')['start'].transform('min')

    # Create a column for the start position of the parent transcript
    df['start_parent_transcript'] = df.groupby('transcript_id')['start'].transform('min')

    # Fill NaNs in start_parent_transcript with start_parent_gene for gene lines
    df['start_parent_transcript'] = df['start_parent_transcript'].fillna(df['start_parent_gene'])

    df = df.copy()
    df.loc[:, 'seqname'] = df['seqname'].astype(str)

    # Sort by seqname, start_parent_gene, gene_id, is_gene, start_parent_transcript, transcript_id, start
    df_sorted = df.sort_values(by=['seqname', 'start_parent_gene', 'gene_id', 'is_gene', 'start_parent_transcript', 'transcript_id', 'start'], ascending=[True, True, True, False, True, True, True]).reset_index(drop=True)

    # Drop the helper columns before writing to the output
    df_sorted = df_sorted.drop(columns=['gene_id', 'transcript_id', 'is_gene', 'start_parent_gene', 'start_parent_transcript'])

    # Write the sorted DataFrame to the output GTF file
    df_sorted.to_csv(output_gtf, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
# reorder_gtf("scratch/Chicken_ensemblv84.Reannotated.gtf", "scratch/Chicken_ensemblv84.Reannotated.sorted.gtf")
# reorder_gtf("GenomeFiles/Chicken_ensemblv84/Reannotated.B.gtf", "scratch/Chicken_ensemblv84.Reannotated.sorted.gtf")

def run_bedparse_gtf2bed(gtf_file, *args):
    """
    wrapper around command line `bedparse gtf2bed`.
    """
    # Construct the bedparse gtf2bed command
    command = ['bedparse', 'gtf2bed', gtf_file]
    # Add any additional optional arguments
    command.extend(args)
    # Run the command and capture the output
    result = subprocess.run(command, capture_output=True, text=True)
    # Check for errors
    if result.returncode != 0:
        print(f"Error running bedparse gtf2bed: {result.stderr}")
        return None
    return result.stdout

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


def insert_marks_for_longset_ORF(sequence, require_ATG=True, require_STOP=True):
    """
    return sequence with "^" to mark start codon, "*" to mark stop for longest ORF. If not require_ATG, can start translation from beginning, which could be useful if transcript starts are not well defined. If not require_STOP, translation does not need stop codon at end. "|" characters (which I use mark splice junctions) are ignored, allowing codons to cross exon junctions.
    """
    if require_ATG and require_STOP:
        regex = r"(?=(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)))"
    elif require_ATG and not require_STOP:
        regex = r"(?=(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)*))"
    elif (not require_ATG) and require_STOP:
        regex = r"(?=((?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)))"
    else:
        regex = r"(?=(\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)*))"
    try:
        longest_orf_match = max(re.findall(regex, sequence, flags=re.IGNORECASE), key = len)
    except ValueError:
        longest_orf_match = None
    if longest_orf_match:
        start_codon_pos = sequence.find(longest_orf_match)
        stop_codon_pos = start_codon_pos + len(longest_orf_match)
        return sequence[0:start_codon_pos] + "^" + longest_orf_match + "*" + sequence[stop_codon_pos:]
    else:
        return sequence
    
def insert_marks_for_first_ORF(sequence, require_STOP=True, min_ORF_len=0):
    """
    return sequence with "^" to mark start codon, "*" to mark stop for longest ORF. If not require_STOP, translation does not need stop codon at end. "|" characters (which I use mark splice junctions) are ignored, allowing codons to cross exon junctions. min_ORF_len in codons after start codon.
    """
    if require_STOP:
        regex = r"^.*?(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN]){" + str(min_ORF_len) + r",}((?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A))).*$"
    else:
        regex = r"^.*?(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN]){" + str(min_ORF_len) + r",}((?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A))?).*$"
    first_ORF_match = re.match(regex, sequence, flags=re.IGNORECASE)
    # import pdb; pdb.set_trace()
    if first_ORF_match:
        start_codon_pos = first_ORF_match.start(1)
        if first_ORF_match.group(2):
            stop_codon_pos = first_ORF_match.end(1)
            return sequence[0:start_codon_pos] + "^" + first_ORF_match.group(1) + "*" + sequence[stop_codon_pos:]
        else:
            return sequence[0:start_codon_pos] + "^" + sequence[start_codon_pos:]
    else:
        return sequence
    
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
# insert_marks_for_defined_ORF("GGGATGAAAGGGAAA|GGG|T|AA|GGGAAA", 3)
# insert_marks_for_defined_ORF("GGGATGAAAGGGAAA|GGG|TT|AA|GGGAAA", 3)

import re
def get_NMD_detective_B_classification(sequence):
    """
    sequence should be marked with '^' for start, '*' for stop, and '|' for splice juncs
    """
    CDS = re.search(r"\^(\w+)\**", sequence.replace("|", ""))
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
# get_NMD_detective_B_classification("ACGTACG|CACGT")
# get_NMD_detective_B_classification("A^ATGACG|CACGT")
# get_NMD_detective_B_classification("A^CGTACG|CAC*GT")
# get_NMD_detective_B_classification("A^CGTA*CG|CACGT")
# get_NMD_detective_B_classification("ACG^TACG|" + "A" * 407 +"A*A|CACGT")
# get_NMD_detective_B_classification("ACG^TA" + "A" * 100 + "CG|" + "A" * 40 +"A*A|CACGT")
# get_NMD_detective_B_classification("ACG^TA" + "A" * 150 + "CG|A*" + "A"*60 + "|CACGT")

def calculate_frames(bedline):
    """Calculate the frame for each CDS block."""
    CDS_bed12 = bedline.cds()
    if CDS_bed12:
        frames = []
        current_frame = 0
        bed6_blocks = list(CDS_bed12.bed12tobed6())
        if bedline.strand == '+':
            for cds in bed6_blocks:
                frames.append(current_frame)
                block_length = cds.end - cds.start
                current_frame = (current_frame + block_length) % 3
        elif bedline.strand == '-':
            for cds in reversed(bed6_blocks):
                frames.insert(0, current_frame)
                block_length = cds.end - cds.start
                current_frame = (current_frame + block_length) % 3
        # Adjusting the frames gets the right answer below. haven't figured out why. but i'm sure it makes the function work.
        frames = [2 if frame == 1 else 1 if frame == 2 else frame for frame in frames]
        return frames
    
def extract_codon(bedline, codon_type='start'):
    """
    Extract the start or stop codon from a BED12 object.
    
    Args:
        bedline (bedparse.bedline): A BED12 object.
        codon_type (str): Either 'start' or 'stop' to specify which codon to extract.
    
    Returns:
        bedparse.bedline: A BED12 object representing the codon.
    """
    if bedline.cds() is None:
        raise ValueError("The BED line does not contain any CDS region.")

    CDS_bed12 = bedline.cds()
    bed6_blocks = list(CDS_bed12.bed12tobed6())
    
    if codon_type == 'start':
        codon_blocks = []
        codon_length = 3
        if bedline.strand == '+':
            for block in bed6_blocks:
                block_length = block.end - block.start
                if block_length >= codon_length:
                    codon_blocks.append((block.start, block.start + codon_length))
                    break
                else:
                    codon_blocks.append((block.start, block.end))
                    codon_length -= block_length
        elif bedline.strand == '-':
            for block in reversed(bed6_blocks):
                block_length = block.end - block.start
                if block_length >= codon_length:
                    codon_blocks.append((block.end - codon_length, block.end))
                    break
                else:
                    codon_blocks.append((block.start, block.end))
                    codon_length -= block_length
            codon_blocks = list(reversed(codon_blocks))

    elif codon_type == 'stop':
        codon_blocks = []
        codon_length = 3
        if bedline.strand == '+':
            for block in reversed(bed6_blocks):
                block_length = block.end - block.start
                if block_length >= codon_length:
                    codon_blocks.append((block.end - codon_length, block.end))
                    break
                else:
                    codon_blocks.append((block.start, block.end))
                    codon_length -= block_length
            codon_blocks = list(reversed(codon_blocks))
        elif bedline.strand == '-':
            for block in bed6_blocks:
                block_length = block.end - block.start
                if block_length >= codon_length:
                    codon_blocks.append((block.start, block.start + codon_length))
                    break
                else:
                    codon_blocks.append((block.start, block.end))
                    codon_length -= block_length

    else:
        raise ValueError("codon_type must be either 'start' or 'stop'.")

    block_sizes = ",".join([str(e - s) for s, e in codon_blocks])
    block_starts = ",".join([str(s - codon_blocks[0][0]) for s, _ in codon_blocks])

    return bedparse.bedline([
        bedline.chr,
        codon_blocks[0][0],
        codon_blocks[-1][1],
        bedline.name,
        bedline.score,
        bedline.strand,
        codon_blocks[0][0],
        codon_blocks[-1][1],
        bedline.color,
        len(codon_blocks),
        block_sizes,
        block_starts
    ])

def gtf_formatted_bedline_tx(bedline, source='.', attributes_str=''):
    """
    input bedparse.bedline transcript object, return 9 field gtf line string, where the source and attributes fields are provided as optional argument, since those information cannot be inferred from the bedline object 
    """
    gtf_fields = [bedline.chr, source, 'transcript', str(bedline.start + 1), str(bedline.end), str(bedline.score), bedline.strand, '.', attributes_str]
    return '\t'.join(gtf_fields) + '\n'

def gtf_formatted_bedline_exons(bedline, source='.', attributes_str=''):
    """
    input bedparse.bedline transcript object, return 9 field gtf line strings, where the source and attributes fields are provided as optional argument, since those information cannot be inferred from the bedline object 
    """
    string_to_return = ""
    for exon in bedline.bed12tobed6(appendExN=True):
        gtf_fields = [exon.chr, source, 'exon', str(exon.start + 1), str(exon.end), str(exon.score), exon.strand, '.', attributes_str + f' exon_id "{exon.name}";']
        string_to_return += '\t'.join(gtf_fields) + '\n'
    return string_to_return

def gtf_formatted_bedline_cds(bedline, source='.', attributes_str=''):
    """
    input bedparse.bedline transcript object, return 9 field gtf line strings, where the source and attributes fields are provided as optional argument, since those information cannot be inferred from the bedline object 
    """
    string_to_return = ""
    if bedline.cds():
        for cds, frame in zip(bedline.cds().bed12tobed6(), calculate_frames(bedline)):
            gtf_fields = [cds.chr, source, 'CDS', str(cds.start + 1), str(cds.end), str(cds.score), cds.strand, str(frame), attributes_str]
            string_to_return += '\t'.join(gtf_fields) + attributes_str + '\n'
    return string_to_return



def gtf_formatted_bedline_utr_start_stop(bedline, source='.', attributes_str='', filter_NF_in_attributes_str=True):
    """
    input bedparse.bedline transcript object, return 9 field gtf line strings, where the source and attributes fields are provided as optional argument, since those information cannot be inferred from the bedline object. Note that start and stop have frame attribute that may be non-zero if ATG cross exon boundary. In Gencode (and maybe ensembl... have to check), when a start codon traverses an exon junction, exists as two lines, one for each piece of the codon. If the attributes string contains 'NF' (as in 'mRNA_end_NF'), don't write out start or stop codon.
    """
    string_to_return = ""
    if filter_NF_in_attributes_str and '_NF' in attributes_str:
        return string_to_return
    bedline_cds = bedline.cds()
    if bedline_cds:
        start = extract_codon(bedline_cds, codon_type='start')
        for cds, frame in zip(start.bed12tobed6(), calculate_frames(start)):
            gtf_fields = [cds.chr, source, 'start_codon', str(cds.start + 1), str(cds.end), str(cds.score), cds.strand, str(frame), attributes_str]
            string_to_return += '\t'.join(gtf_fields) + '\n'
        stop = extract_codon(bedline_cds, codon_type='stop')
        for cds, frame in zip(stop.bed12tobed6(), calculate_frames(stop)):
            gtf_fields = [cds.chr, source, 'stop_codon', str(cds.start + 1), str(cds.end), str(cds.score), cds.strand, str(frame), attributes_str]
            string_to_return += '\t'.join(gtf_fields) + '\n'
        UTR_fiveprime = bedline.utr(which=5)
        if UTR_fiveprime:
            for exon in UTR_fiveprime.bed12tobed6():
                gtf_fields = [exon.chr, source, 'UTR', str(exon.start + 1), str(exon.end), str(exon.score), exon.strand, '.', attributes_str]
                string_to_return += '\t'.join(gtf_fields) + '\n'
        UTR_threeprime = bedline.utr(which=3)
        if UTR_threeprime:
            for exon in UTR_threeprime.bed12tobed6():
                gtf_fields = [exon.chr, source, 'UTR', str(exon.start + 1), str(exon.end), str(exon.score), exon.strand, '.', attributes_str]
                string_to_return += '\t'.join(gtf_fields) + '\n'
    return string_to_return

def bed12_formatted_bedline(bedline, attributes_str=''):
    bed12_fields = [bedline.chr, bedline.start, bedline.end, bedline.name, bedline.score, bedline.strand, bedline.cdsStart, bedline.cdsEnd, bedline.color, bedline.nEx, bedline.exLengths, bedline.exStarts]
    return '\t'.join([str(i) for i in bed12_fields]) + attributes_str + '\n'
    
def get_transcript_length(bedline):
    return sum([int(i) for i in bedline.exLengths.split(',')])


def get_thickStart_thickStop_from_marked_seq(bedline, ORF_marked_sequence):
    CDS = re.search(r"\^\w+\*", ORF_marked_sequence.replace("|", ""))
    if CDS:
        TranscriptLength = sum([int(i) for i in bedline.exLengths.split(',')])
        relativeStart = CDS.span()[0]
        relativeStop = CDS.span()[1] - 2
        FivePrimeEdge = get_absolute_pos(bedline, relativeStart)
        ThreePrimeEdge = get_absolute_pos(bedline, relativeStop)
        return tuple(sorted([FivePrimeEdge, ThreePrimeEdge]))
    return bedline.cdsStart, bedline.cdsEnd

def Fill_A_and_B_to_same(A, B):
    if A == "." and B != ".":
        A = B
    elif B == "." and A != ".":
        B = A
    elif B != A:
        B = A
    return(A, B)

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

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Process GTF and FASTA files with optional transcript attributes. By default, if gene_type is protein_coding but no CDS is annotated (or if coding sequence ends at end of transcript) will manually search for longest ORF and add CDS to output gtf along with NMDFinderB transcript classification. The output gtf will have start_codon and stop_codon features, allowing use with Yangs script for classifying splice junctions as productive or unproductive.")
    parser.add_argument('-i', dest='gtf_in', required=True, help='Input GTF file')
    parser.add_argument('-o', dest='gtf_out', required=True, help='Output GTF file')
    parser.add_argument('-fa', dest='fasta_in', required=True, help='Input FASTA file')
    parser.add_argument('-bed12_out', dest='bed12_out', required=False, help='Optional bed12 out. One line per transcript. Attributes as extra columns')
    parser.add_argument('-extra_attributes', dest='extra_attributes', default='transcript_support_level,tag,ccds_id',
                        help='Extra transcript attributes (comma delimited quoted string). default: "%(default)s"')
    parser.add_argument('-translation_approach', dest='translation_approach', choices=['A', 'B', 'C', 'D', 'E'], default='a', help='Approach to annotate CDS in output for genes where gene_biotype/gene_type == "protein_coding", some of which may not have annotated CDS in input (eg transcript_type=="processed_transcript"). Possible approaches: (A) using annotated ORF if ORF is present in input with 5UTR and 3UTR. If 3UTR is present and 5UTR is absent (suggesting stop codon is annotated but start codon may be outside of the transcript bounds), search for longest ORF within transcript bounds where a start codon is not required at the beginning of ORF. Similarly, if 5UTR is present but 3UTR is absent, search for longest ORF without requiring stop codon. If neither UTR is present, or if no CDS is annotated in input, search for longest ORF, not requiring start or stop within transcript bounds. I think this approach might be useful to correctly identify the ORF, even if transcript bounds are not accurate, but it has the downside that true "processed_transcripts" with a internal TSS that eliminates the correct start codon, may be erronesously be classified as "Last exon" (ie productive) transcripts by NMDFinderB. (B) Use only annotated CDS. In effect, output gtf is the same except start_codon and stop_codon features are added even if not present in input. This would be useful for dealing with "processed_transcripts" properly by NMDFinder, but I havent checked whether CDS annotations in poorly annotated species (eg lamprey, chicken, etc) are reasonable, which could be a problem for Yangs script. (C) use annotated CDS if present, and use first ATG if no CDS present (minimum ORF length of 30 codons, not including start or stop). (D) Same as (C) but no minimum ORF length. (E) Use first ORF with minimum length > 42, regardless of annotation')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase verbosity')
    return parser.parse_args(args)


def setup_logging(verbose):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

def main(args=None):
    args = parse_args(args)
    setup_logging(args.verbose)

    # Example usage of the parsed arguments
    logging.info("Input GTF file: %s", args.gtf_in)
    logging.info("Output GTF file: %s", args.gtf_out)
    logging.info("Input FASTA file: %s", args.fasta_in)
    if args.extra_attributes:
        extra_attributes = args.extra_attributes.split(',')
        logging.info("Extra attributes: %s", extra_attributes)
    else:
        logging.info("No extra attributes provided.")
    
    # Open the FASTA file using pysam
    fasta_obj = pysam.FastaFile(args.fasta_in)

    required_attributes = 'gene_id,gene_name,transcript_id,transcript_name,gene_biotype,gene_type,transcript_biotype,transcript_type'
    attributes_to_extract = ','.join([required_attributes, args.extra_attributes])

    # gtf2bed
    logging.info(f'Reading in gtf and converting transcripts to bedparse.bedline objects...')
    bed12 = run_bedparse_gtf2bed(args.gtf_in, '--extraFields', attributes_to_extract)
    NumTrancsripts = len(bed12.splitlines())
    logging.info(f'Read in {NumTrancsripts} from gtf.')

    gtf_stringio = StringIO()
    # Initialize an empty dictionary to keep coordinates of each transcript for each gene. Will need to write out gene coordinates as most extensive span of child transcripts.
    gene_dict = defaultdict(list)
    bed_out_fh = None
    if args.bed12_out:
        bed_out_fh = open(args.bed12_out, 'w')
    for i,l in enumerate(bed12.splitlines()):
        if i >= 0:
            # print(i)
            if i % 1000 == 0: logging.debug(f'processed {i} trancsripts for output...')
            lsplit = l.split('\t')
            # required_attributes = 'gene_id,gene_name,transcript_id,transcript_name,gene_biotype,gene_type,transcript_biotype,transcript_type'

            gene_id, gene_name, transcript_id, transcript_name, gene_biotype, gene_type, transcript_biotype, transcript_type = lsplit[12:20]
            gene_biotype, gene_type = Fill_A_and_B_to_same(gene_biotype, gene_type)
            transcript_biotype, transcript_type = Fill_A_and_B_to_same(transcript_biotype, transcript_type)
            gene_id, gene_name = Fill_A_and_B_to_same(gene_id, gene_name)
            transcript_id, transcript_name = Fill_A_and_B_to_same(transcript_id, transcript_name)
            extra_attribute_values = lsplit[20:]
            transcript_attributes = ' '.join([f'{tag} "{value}";' for (tag, value) in zip(attributes_to_extract.split(','), ([gene_id, gene_name, transcript_id, transcript_name, gene_biotype, gene_type, transcript_biotype, transcript_type] + extra_attribute_values)) if value!='.'])
            transcript = bedparse.bedline(lsplit[0:12])    
            # import pdb; pdb.set_trace()
            gene_dict[gene_id].append((transcript.chr, transcript.start, transcript.end, transcript.strand))
            source = "input_gtf"
            transcript_out = transcript
            NMDFinderB = "NA"
            if args.translation_approach == 'A':
                if gene_biotype == "protein_coding":
                    if transcript.cds() and transcript.utr(which=5) and transcript.utr(which=3):
                        source = "input_gtf"
                        sequence = transcript.extract_sequence(fasta_obj, AddMarksForORF=True)
                    elif transcript.cds() and transcript.utr(which=5) and not transcript.utr(which=3):
                        source = "LongestORF_NoStopRequired"
                        sequence = insert_marks_for_longset_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    elif transcript.cds() and not transcript.utr(which=5) and transcript.utr(which=3):
                        source = "LongestORF_NoStartRequired"
                        sequence = insert_marks_for_longset_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    else:
                        source = "LongestORF_NeitherRequired"
                        sequence = insert_marks_for_longset_ORF(transcript.extract_sequence(fasta_obj), require_ATG=False, require_STOP = False)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    NMDFinderB = get_NMD_detective_B_classification(sequence).replace(' ', '_')
            elif args.translation_approach == 'B':
                if gene_biotype == "protein_coding":
                    sequence = transcript.extract_sequence(fasta_obj, AddMarksForORF=True)
                    NMDFinderB = get_NMD_detective_B_classification(sequence).replace(' ', '_')
            elif args.translation_approach == 'C':
                if gene_biotype == "protein_coding":
                    if transcript.cds():
                        sequence = transcript.extract_sequence(fasta_obj, AddMarksForORF=True)
                        source = "input_gtf"
                    else:
                        source = "FirstORF_NoStopRequired"
                        sequence = insert_marks_for_first_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False, min_ORF_len = 30)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    NMDFinderB = get_NMD_detective_B_classification(sequence).replace(' ', '_')
            elif args.translation_approach == 'D':
                if gene_biotype == "protein_coding":
                    if transcript.cds():
                        sequence = transcript.extract_sequence(fasta_obj, AddMarksForORF=True)
                        source = "input_gtf"
                    else:
                        source = "FirstORF_NoStopRequired"
                        sequence = insert_marks_for_first_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False, min_ORF_len = 0)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    NMDFinderB = get_NMD_detective_B_classification(sequence).replace(' ', '_')
            elif args.translation_approach == 'E':
                source = "FirstORF_NoStopRequired"
                sequence = insert_marks_for_first_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False, min_ORF_len = 42)
                thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                NMDFinderB = get_NMD_detective_B_classification(sequence).replace(' ', '_')
            if transcript_out.cds():
                CDSLen = get_transcript_length(transcript_out.cds())
            else:
                CDSLen = 0
            transcript_attributes += f' tag "NMDFinderB:{NMDFinderB}"; tag "CDSLen:{CDSLen}";'
            _ = gtf_stringio.write(gtf_formatted_bedline_tx(transcript_out, source=source, attributes_str=transcript_attributes))
            _ = gtf_stringio.write(gtf_formatted_bedline_exons(transcript_out, source=source, attributes_str=transcript_attributes))
            _ = gtf_stringio.write(gtf_formatted_bedline_cds(transcript_out, source=source, attributes_str=transcript_attributes))
            _ = gtf_stringio.write(gtf_formatted_bedline_utr_start_stop(transcript_out, source=source, attributes_str=transcript_attributes))
            if bed_out_fh is not None:
                _ = bed_out_fh.write(bed12_formatted_bedline(transcript, attributes_str='\t' + '\t'.join([gene_id, transcript_id, gene_biotype, gene_type, transcript_biotype, transcript_type] + extra_attribute_values + [NMDFinderB, str(CDSLen)])))
    if bed_out_fh is not None: bed_out_fh.close()

    # Identify the gene coordinates as the widest coordinates among all the transcripts belonging to that gene
    gene_coordinates = {}
    for gene, coords in gene_dict.items():
        chroms = set(c[0] for c in coords)
        if len(chroms) > 1:
            raise ValueError(f"Transcripts for gene {gene} are on different chromosomes.")
        chrom = chroms.pop()
        strands = set(c[3] for c in coords)
        if len(strands) > 1:
            raise ValueError(f"Transcripts for gene {gene} are on different strands.")
        strand = strands.pop()
        start = min(c[1] for c in coords)
        stop = max(c[2] for c in coords)
        gene_coordinates[gene] = (chrom, start, stop, strand)

    # Write out gtf lines for gene features. Only gene_id attribute included
    for gene, coords in gene_coordinates.items():
        chrom, start, stop, strand = coords
        _ = gtf_stringio.write(f'{chrom}\tinput_gtf\tgene\t{start+1}\t{stop}\t.\t{strand}\t.\tgene_id "{gene}"; gene_name "{gene}";\n')

    logging.info('sorting and writing out gtf')
    # output filehandle
    gtf_stringio.seek(0)
    with open(args.gtf_out, 'w') as output_fh:
        # transfer commented headers from original gtf
        with open(args.gtf_in, 'r') as input_fh:
            for l in input_fh:
                if l.startswith('#'):
                    _ = output_fh.write(l)
                else:
                    break
            _ = output_fh.write(f"#! args: {args}\n")
        reorder_gtf(gtf_stringio, output_fh)


if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        # main("-i /project2/yangili1/bjf79/2024_comparativesplicing/code/GenomeFiles/Human_ensemblv75/Reference.gtf -o scratch/Human_ensemblv75.Reannotated.gtf -fa /project2/yangili1/bjf79/ReferenceGenomes/Human_ensemblv75/Reference.fa -v -translation_approach A".split(' '))
        main("-i /project2/yangili1/bjf79/2024_comparativesplicing/code/GenomeFiles/Chicken_ensemblv84/Reference.gtf -o scratch/Chicken_ensemblv84.Reannotated.gtf -fa /project2/yangili1/bjf79/ReferenceGenomes/Chicken_ensemblv84/Reference.fa -v -translation_approach A".split(' '))
        # main("-i /project2/yangili1/bjf79/2024_comparativesplicing/code/GenomeFiles/Human_ensemblv75/Reference.gtf -o scratch/Human_ensemblv75.Reannotated.gtf -fa /project2/yangili1/bjf79/ReferenceGenomes/Human_ensemblv75/Reference.fa -v -translation_approach A".split(' '))
    else:
        main()



# args=parse_args("-i /project2/yangili1/bjf79/2024_comparativesplicing/code/GenomeFiles/Human_ensemblv75/Reference.gtf -o scratch/Human_ensemblv75.Reannotated.gtf -fa /project2/yangili1/bjf79/ReferenceGenomes/Human_ensemblv75/Reference.fa -v -translation_approach A".split(' '))

# main("-i /project2/yangili1/bjf79/2024_comparativesplicing/code/GenomeFiles/Human_ensemblv75/Reference.gtf -o scratch/Human_ensemblv75.Reannotated.gtf -fa /project2/yangili1/bjf79/ReferenceGenomes/Human_ensemblv75/Reference.fa -v -translation_approach A".split(' '))
