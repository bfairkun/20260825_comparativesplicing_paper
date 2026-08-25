import argparse
import sys
import logging
from Bio.AlignIO.MafIO import MafIndex
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from collections import Counter

def setup_logging(verbose):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Calculate dN/dS from an indexed MAF file.")
    parser.add_argument('--index', dest="maf_index", help="Path to the MAF index file", required=True)
    parser.add_argument('--maf', dest="maf_file", help="Path to the MAF file", required=True)
    parser.add_argument('--ref', dest="reference", help="Reference genome name in MAF file", required=True)
    parser.add_argument('--start', dest="start", type=int, help="Start position of the interval (1-based, fully closed, as in IGV display)", required=True)
    parser.add_argument('--stop', dest="stop", type=int, help="Stop position of the interval (1-based, fully closed, as in IGV display)", required=True)
    parser.add_argument('--strand', dest="strand", choices=["+", "-"], help="Strand of the sequence", required=True)
    parser.add_argument('--frame', dest="frame", type=int, choices=[0, 1, 2], help="Reading frame (0, 1, or 2)", required=True)
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase verbosity')
    return parser.parse_args(args)

def extract_alignment(maf_index, start, stop):
    results = maf_index.search([start], [stop])
    alignments = []
    
    for aln in results:
        ref_record = next((rec for rec in aln if rec.id.startswith(maf_index._target_seqname)), None)
        if ref_record is None:
            continue  # Skip if no reference sequence found

        # Extract original start position from annotations
        seq_start = ref_record.annotations.get("start", None)
        if seq_start is None:
            raise KeyError(f"Missing 'start' annotation in MAF for reference sequence {ref_record.id}")

        # Determine column indices for exact start:stop range
        col_start = max(0, start - seq_start)
        col_end = max(0, stop - seq_start)

        # Slice the alignment along its columns
        trimmed_alignment = aln[:, col_start:col_end]

        # Update the reference record's start annotation (correct after trimming)
        for rec in trimmed_alignment:
            if rec.id == ref_record.id:
                rec.annotations["start"] = start  # Set to requested `start`, now correct after slicing

        alignments.append(trimmed_alignment)

    return alignments

def calculate_dn_ds(alignments, reference, strand, frame=0):
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    dn, ds = 0, 0
    current_frame = frame  # Start with user-provided frame shift
    
    # Convert reversed alignments to a list so we can index it
    alignments = list(alignments) if strand == "+" else list(reversed(alignments))
    
    for i, aln in enumerate(alignments):
        ref_seq = None
        species_seqs = {}
        ref_start = None  # Store the reference start position for frame tracking

        for record in aln:
            seq = str(record.seq).replace("-", "")  # Remove gaps
            if record.id.startswith(reference):
                ref_seq = seq
                ref_start = record.annotations["start"]  # Extract reference start
            else:
                species_seqs[record.id] = seq
        
        if not ref_seq:
            continue  # Skip if no reference sequence

        # Reverse complement for minus strand cases (AFTER extracting start)
        if strand == "-":
            ref_seq = str(Seq(ref_seq).reverse_complement())
            species_seqs = {sp: str(Seq(seq).reverse_complement()) for sp, seq in species_seqs.items()}

        # Apply the current frame shift
        ref_seq = ref_seq[current_frame:]
        species_seqs = {sp: seq[current_frame:] for sp, seq in species_seqs.items()}
        
        codons = [ref_seq[i:i+3] for i in range(0, len(ref_seq), 3) if len(ref_seq[i:i+3]) == 3]
        
        for species, seq in species_seqs.items():
            species_codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
            for ref_codon, species_codon in zip(codons, species_codons):
                if ref_codon in standard_table.stop_codons or species_codon in standard_table.stop_codons:
                    continue
                if ref_codon != species_codon:
                    if Seq(ref_codon).translate() == Seq(species_codon).translate():
                        ds += 1
                    else:
                        dn += 1
        
        # If there's another alignment, update the frame shift
        if i < len(alignments) - 1:  # Check if this isn't the last alignment
            next_ref_start = alignments[i + 1][0].annotations["start"]
            shift = (next_ref_start - ref_start) % 3
            current_frame = shift  # Update frame shift for next alignment
    return dn / ds if ds > 0 else float("inf")

def main(args=None):
    args = parse_args(args)
    setup_logging(args.verbose)
    #Convert to 0-based start coord
    args.start = args.start - 1
    maf_index = MafIndex(args.maf_index, args.maf_file, args.reference)
    alignments = extract_alignment(maf_index, args.start, args.stop)
    # breakpoint()
    dn_ds_ratio = calculate_dn_ds(alignments, args.reference, args.strand, frame=args.frame)
    logging.debug(f"dN/dS ratio: {dn_ds_ratio}")
    return(dn_ds_ratio)

if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        # some tests
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 10510027 --stop 10510127 --strand - --frame 1 -v".split(' '))
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 10510027 --stop 10510127 --strand + --frame 1".split(' '))
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 10510027 --stop 10510127 --strand + --frame 0".split(' '))
        #THOC5 coding exons
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 29531831 --stop 29531963 --strand - --frame 0".split(' '))
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 29529162 --stop 29529239 --strand - --frame 2".split(' '))
        #THOC5 coding exons, wrong frame
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 29531831 --stop 29531963 --strand - --frame 1".split(' '))
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 29529162 --stop 29529239 --strand - --frame 1".split(' '))
        #RBFOX2 coding exons
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 29531831 --stop 29531963 --strand - --frame 0".split(' '))
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 35768257 --stop 35768349 --strand - --frame 0".split(' '))
        #An RBFOX2 coding exon split amongst 2 alignments
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 35746473 --stop 35746561 --strand - --frame 0".split(' '))
        #poison exon1
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 35772931 --stop 35772971 --strand - --frame 0".split(' '))
        #poison exon2
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 35767368 --stop 35767441 --strand - --frame 0".split(' '))
        #KCTD17 coding exon
        main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 37053100 --stop 37053208 --strand + --frame 0".split(' '))
    main()
