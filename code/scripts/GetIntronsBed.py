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
import gzip
import pysam
import subprocess
import argparse


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


def extract_sequence(self, fasta_obj):
    """
    ...method to be used on bedparse.bedline object. will monkey patch it into bedparse.bedline class
    will also add '|' between blocks to mark exon junctions
    """
    # Extract information from the BED entry
    chrom = self.chr
    start = self.start
    end = self.end
    strand = self.strand if hasattr(self, 'strand') else '+'
    # block_sizes = list(map(int, self.exLengths.split(','))) if self.exLengths else [self.end - self.start]
    # block_starts = list(map(int, self.exStarts.split(','))) if self.exStarts else [0]
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
    return sequence

def reverse_complement(seq):
    # Helper function to reverse complement a DNA sequence
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]

# Monkey-patch the BEDLine class
bedparse.bedline.extract_sequence = extract_sequence


def main(gtf_file, fasta_file, output_bed):
    """Main function to parse GTF, extract introns, and write to BED."""
    fasta_obj = pysam.FastaFile(fasta_file)
    bed12 = run_bedparse_gtf2bed(gtf_file, '--extraFields', 'gene_id,transcript_id,gene_type,gene_name,transcript_type,transcript_support_level,tag,transcript_name,gene_biotype,transcript_biotype')
    with open(output_bed, 'w') as out_bed:
        for i,l in enumerate(bed12.splitlines()):
            lsplit = l.split('\t')
            extras = lsplit[12:]
            bl = bedparse.bedline(lsplit[0:12])
            if bl.introns():
                for intron in bl.introns().bed12tobed6(appendExN=True):
                    if intron.strand == '-':
                        donor = bedparse.bedline([intron.chr, intron.end - 6, intron.end +3, intron.name + '_donor', intron.score, intron.strand])
                        acceptor = bedparse.bedline([intron.chr, intron.start - 3, intron.start +20, intron.name + '_acceptor', intron.score, intron.strand])
                    else:
                        donor = bedparse.bedline([intron.chr, intron.start - 3, intron.start +6, intron.name + '_donor', intron.score, intron.strand])
                        acceptor = bedparse.bedline([intron.chr, intron.end - 20, intron.end +3, intron.name + '_acceptor', intron.score, intron.strand])
                    out_bed.write('\t'.join([str(field) for field in [intron.chr, intron.start, intron.end, intron.name, intron.score, intron.strand, donor.extract_sequence(fasta_obj), acceptor.extract_sequence(fasta_obj)]]) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract introns from a GTF file and write to a BED file with splice site sequences.")
    parser.add_argument("gtf_file", help="Input GTF file")
    parser.add_argument("fasta_file", help="Reference genome FASTA file")
    parser.add_argument("output_bed", help="Output BED file with introns and splice site sequences")
    args = parser.parse_args()
    main(args.gtf_file, args.fasta_file, args.output_bed)