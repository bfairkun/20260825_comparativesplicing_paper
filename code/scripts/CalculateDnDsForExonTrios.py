import pandas as pd
import sys
sys.path.append('./scripts')
from MAF_to_dnds import main

# Read the file with pandas
df = pd.read_csv('./scratch/dnds_exon_trios.tsv.gz', sep='\t')
df['ce_frame'] = df.apply(lambda row: row['downstream_frame'] if row['strand'] == '+' else row['upstream_frame'], axis=1)


pd.set_option('display.max_columns', None)

# Print the first few rows of the DataFrame for context
print(df.head())



# Define a function to apply the main function to each row
def apply_main_function(row, pos_col, frame_col):
    chrom, pos_range = row[pos_col].split(":")
    start, stop = pos_range.split("-")
    start_int, stop_int = int(start), int(stop)
    
    # Check if the region has positive length
    if stop_int - start_int <= 0:
        return pd.Series([float('nan')])
    
    args = [
        '--index', f"multiple_sequence_alignment_analysis/maf/{chrom}.mafindex",
        '--maf', f"multiple_sequence_alignment_analysis/maf/{chrom}.maf",
        '--ref', f"hg38.{chrom}",
        '--start', start,
        '--stop', stop,
        '--strand', row['strand'],
        '--frame', str(row[frame_col])
    ]
    # print(args)
    result = main(args)
    return pd.Series([result])

main("--index scratch/chr22.mafindex --maf scratch/chr22.maf --ref hg38.chr22 --start 37053100 --stop 37053208 --strand + --frame 0".split(' '))
main(['--index', 'multiple_sequence_alignment_analysis/maf/chrX.mafindex', '--maf', 'multiple_sequence_alignment_analysis/maf/chrX.maf', '--ref', 'hg38.chrX', '--start', '100630798', '--stop', '100630866', '--strand', '-', '--frame', 0])
main(['--index', 'multiple_sequence_alignment_analysis/maf/chr20.mafindex', '--maf', 'multiple_sequence_alignment_analysis/maf/chr20.maf', '--ref', 'hg38.chr20', '--start', '50936148', '--stop', '50936262', '--strand', '-', '--frame', '1'])

# Apply the function to each row and add the new column to the dataframe
df_test = df.head(10)
df_test['dn_ds_ratio'] = df_test.apply(apply_main_function, axis=1, args=('upstream_exon', 'upstream_frame'))

df_test2 = df.loc[df['gene_name'] == 'ENSG00000196876.19']
df_test2['dn_ds_ratio'] = df_test2.apply(apply_main_function, axis=1, args=('upstream_exon', 'upstream_frame'))
df_test2

df.head(10).apply(apply_main_function, axis=1, args=('upstream_exon', 'upstream_frame'))

df['dn_ds_ratio_us'] = df.apply(apply_main_function, axis=1, args=('upstream_exon', 'upstream_frame'))
df['dn_ds_ratio_ds'] = df.apply(apply_main_function, axis=1, args=('downstream_exon', 'downstream_frame'))
df['dn_ds_ratio_ce'] = df.apply(apply_main_function, axis=1, args=('cassette_exon', 'ce_frame'))


# Save the updated dataframe to a new file
df.to_csv('./scratch/dnds_exon_trios.dnds_added.tsv.gz', sep='\t', index=False)