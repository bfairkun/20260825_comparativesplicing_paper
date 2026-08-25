import argparse
import pyBigWig
import pandas as pd
import sys

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Append mean BigWig signal to BED intervals.")
    parser.add_argument('-b', '--bedfile', required=True, help='Input BED file')
    parser.add_argument('-w', '--bigwigfile', required=True, help='Input BigWig file')
    parser.add_argument('-o', '--outfile', required=True, help='Output BED file with mean signal appended')
    return parser.parse_args(args)

def calculate_mean_signal(row, bw):
    try:
        signal = bw.values(row[0], row[1], row[2], numpy=True)
        if len(signal) == 0:
            return 0.0
        return signal.mean()
    except RuntimeError as e:
        # Handle invalid interval bounds or other runtime errors
        print(f"Error processing interval {row[0]}:{row[1]}-{row[2]}: {e}")
        return 0.0

def process_bed_file(bed_file, bigwig_file, output_file):
    """
    Process a BED file and append the mean BigWig signal to each interval.
    """
    bw = pyBigWig.open(bigwig_file)
    chunk_size = 10000  # Adjust chunk size based on memory constraints and performance
    # Open the output file for writing
    with open(output_file, 'w') as out_file:
        # Read the BED file in chunks
        reader = pd.read_csv(bed_file, sep='\t', header=None, chunksize=chunk_size)
        for chunk in reader:
            # Add a new column with the mean signal
            chunk[3] = chunk.apply(lambda row: calculate_mean_signal(row, bw), axis=1)
            # Write the chunk to the output file
            chunk.to_csv(out_file, sep='\t', header=False, index=False, mode='a')
    bw.close()

def main(args=None):
    """
    Main function to append mean BigWig signal to BED intervals.
    """
    # Parse command-line arguments
    parsed_args = parse_args(args)
    # Extract the arguments
    bed_file = parsed_args.bedfile
    bigwig_file = parsed_args.bigwigfile
    output_file = parsed_args.outfile
    # Process the BED file and calculate mean signal
    print(f"Processing BED file: {bed_file}")
    print(f"Using BigWig file: {bigwig_file}")
    print(f"Output will be saved to: {output_file}")
    # Call the processing function
    process_bed_file(bed_file, bigwig_file, output_file)
    print("Processing complete.")

if __name__ == "__main__":
    if hasattr(sys, 'ps1'):
        main("-b conservation/human.AS_segments.bed -w conservation/PhyloP.hg19.bw -o scratch/test.bed".split(' '))
    main()


