rule CopyGtf_AndMakeSameFormat:
    """
    Different ensembl versions are formatted a bit differently
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        "GenomeFiles/{GenomeName}/Reference.gtf"
    wildcard_constraints:
        GenomeName = "|".join([g for g in STAR_genomes.index if g != "Human_ensemblv75"])
    shell:
        """
        cp {input} {output}
        """

rule CopyGtf_AndMakeSameFormat_human:
    """
    Different ensembl versions are formatted a bit differently
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        "GenomeFiles/{GenomeName}/Reference.gtf"
    wildcard_constraints:
        GenomeName = "Human_ensemblv75"
    shell:
        """
        awk '$1 ~ /^#/ {{print}} $3 == "gene" && $1 !~ /^#/ {{print}} $3 != "gene" && $1 !~ /^#/ {{print $0" transcript_biotype \\""$2"\\";"}}' {input} > {output}
        """

rule gtf_to_bed12_copiedgtf:
    """
    For ensembl gtfs, which may need some manual checking for formatting... Different ensembl versions are formatted a bit differently
    """
    input:
        gtf = "GenomeFiles/{GenomeName}/Reference.gtf"
    output:
        bed = "GenomeFiles/{GenomeName}/Reference.bed.gz",
        index = touch("GenomeFiles/{GenomeName}/Reference.bed.gz.indexing_done")
    params:
        tabix_params = GetIndexingParams
    log:
        "logs/gtf_to_bed12_copiedgtf/{GenomeName}.log"
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        (bedparse gtf2bed {input.gtf} --extraFields gene_id,transcript_id,gene_biotype,gene_name,transcript_biotype,transcript_support_level,tag,transcript_name | awk -F'\\t' -v OFS='\\t' '{{$4=$NF; print $0}}' | bedtools sort -i - | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        tabix {params.tabix_params} -f -p bed {output.bed} && touch {output.index}
        """

rule GetSpliceSitesBedSeq:
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    output:
        "GenomeFiles/{GenomeName}/SpliceSites.seq.bed.gz"
    log:
        "logs/GetSpliceSitesBedSeq/{GenomeName}.log"
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        (python scripts/GetIntronsBed.py {input.gtf} {input.fa} /dev/stdout | gzip - > {output}) &> {log}
        """

rule Reannotate_gtf:
    input:
        gtf = "GenomeFiles/{GenomeName}/Reference.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    output:
        bed = "GenomeFiles/{GenomeName}/Reannotated.{TranslationApproach}.bed",
        gtf = "GenomeFiles/{GenomeName}/Reannotated.{TranslationApproach}.gtf"
    conda:
        "../envs/bedparse.yml"
    resources:
        mem_mb = 16000
    log:
        "logs/Reannotate_gtf/{GenomeName}.{TranslationApproach}.log"
    shell:
        """
        python scripts/bedparse_translate_transcripts.py -i {input.gtf} -fa {input.fa} -o {output.gtf} -bed12_out {output.bed} -translation_approach {wildcards.TranslationApproach} -v &> {log}
        """

rule SortAndTabixReannotated:
    input:
        bed = "GenomeFiles/{GenomeName}/Reannotated.{TranslationApproach}.bed",
        gtf = "GenomeFiles/{GenomeName}/Reannotated.{TranslationApproach}.gtf"
    output:
        bed = "GenomeFiles/{GenomeName}/Reannotated.{TranslationApproach}.bed.gz",
        gtf = "GenomeFiles/{GenomeName}/Reannotated.{TranslationApproach}.gtf.gz",
        bed_index = touch("GenomeFiles/{GenomeName}/Reannotated.{TranslationApproach}.bed.gz.indexing_done"),
        gtf_index = touch("GenomeFiles/{GenomeName}/Reannotated.{TranslationApproach}.gtf.gz.indexing_done")
    log:
        "logs/SortAndTabixReannotated/{GenomeName}.{TranslationApproach}.log"
    params:
        tabix_params = GetIndexingParams
    resources:
        mem_mb = 16000
    shell:
        """
        bedtools sort -i {input.bed} | bgzip -c /dev/stdin > {output.bed}
        awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k4,4n -k5,5n"}}' {input.gtf} | bgzip -c /dev/stdin > {output.gtf}
        tabix {params.tabix_params} -f -p bed {output.bed} && touch {output.bed_index}
        tabix {params.tabix_params} -f -p gff {output.gtf} && touch {output.gtf_index}
        """

rule MakeMazin_AS_SegmentBeds:
    input:
        "kaessman_AS_dat/Supplementary_Data/Supplementary_Data_9.csv"
    output:
        expand("kaessman_AS_dat/AS_segment_lists/{GenomeName}.bed", GenomeName = MazinGenomes)
    params:
        Prefix = "kaessman_AS_dat/AS_segment_lists/"
    log:
        "logs/MakeMazin_AS_SegmentBeds"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/Mazin_AS_Segments_ToBed9.R {input} {params.Prefix} &> {log}
        """

rule NMDFinderB_EffectsOfMazinSegments:
    input:
        AS_segments_bed = "kaessman_AS_dat/AS_segment_lists/{GenomeName}.bed",
        transcripts_bed = "GenomeFiles/{GenomeName}/Reannotated.B.bed.gz",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    log:
        "logs/NMDFinderB_EffectsOfMazinSegments/{GenomeName}.log"
    output:
        junclist = "kaessmanAnalysis/JuncListsForLeafcutter2/{GenomeName}.juncs.gz",
        results = "kaessmanAnalysis/NMDFinderB_ASSegments/{GenomeName}.tsv.gz"
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        python scripts/AS_ExonSegments_to_AS_Junctions.py -TsvOut {output.results} -JuncListOut {output.junclist} -AS_segmentsIn {input.AS_segments_bed} -AnnotatedTranscriptsIn {input.transcripts_bed} -fa {input.fa} -v &> {log}
        """

rule leafcutter2_ClassifyJuncs:
    input:
        junclist = "kaessmanAnalysis/JuncListsForLeafcutter2/{GenomeName}.juncs.gz",
        gtf = "GenomeFiles/{GenomeName}/Reannotated.B.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    output:
        "kaessmanAnalysis/JuncListsForLeafcutter2/{GenomeName}._junction_classifications.txt"
    log:
        "logs/leafcutter2_ClassifyJuncs/{GenomeName}.log"
    params:
        rundir = "kaessmanAnalysis/JuncListsForLeafcutter2"
    resources:
        mem_mb = 16000
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        python scripts/leafcutter2_bfairkun/scripts/SpliceJunctionClassifier.py -c {input.junclist} -G {input.fa} -A {input.gtf} -v -r {params.rundir} -o {wildcards.GenomeName}. &> {log}
        """

rule leafcutter_ClassifyJuncs_to_AS_segments:
    input:
        transcripts_bed = "GenomeFiles/{GenomeName}/Reannotated.B.bed.gz",
        junc_classifications = "kaessmanAnalysis/JuncListsForLeafcutter2/{GenomeName}._junction_classifications.txt",
        AS_segments = "kaessmanAnalysis/NMDFinderB_ASSegments/{GenomeName}.tsv.gz"
    output:
        JuncClassifications_Full = "kaessmanAnalysis/leaf2_to_AS_segments/{GenomeName}.full.tsv.gz",
        JuncClassifications_Collapsed = "kaessmanAnalysis/leaf2_to_AS_segments/{GenomeName}.collapsed.tsv.gz",
    log:
        "logs/leafcutter_ClassifyJuncs_to_AS_segments/{GenomeName}.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/LeafcutterJuncClassificationsToAS_segments.R {input.junc_classifications} {input.AS_segments} {input.transcripts_bed} {output.JuncClassifications_Full} {output.JuncClassifications_Collapsed} &> {log}
        """

rule Download_hg19_PhyloP:
    output:
        "conservation/PhyloP.hg19.bw"
    shell:
        """
        wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
        """

rule Get_human_Mazin_AS_Segments:
    input:
        "kaessman_AS_dat/Supplementary_Data/Supplementary_Data_9.csv"
    output:
        "conservation/human.AS_segments.bed"
    run:
        # Read the CSV file with pandas
        df = pd.read_csv(input[0], quotechar='"', skipinitialspace=True)

        # Filter for rows where the first column starts with "hum"
        df_filtered = df[df.iloc[:, 0].str.startswith('hum')]

        # Select and rename the relevant columns
        df_filtered = df_filtered.iloc[:, [3, 4, 5, 0]]
        df_filtered.columns = ['chrom', 'start', 'end', 'name']

        # prepend chr to chromosomes
        df_filtered['chrom'] = 'chr' + df_filtered['chrom'].astype(str)

        # Remove quotes from the 'name' column
        df_filtered['name'] = df_filtered['name'].str.replace('"', '')

        # Filter out invalid intervals where end <= start
        df_filtered = df_filtered[df_filtered['end'] > df_filtered['start']]

        # Write the output to a BED file
        df_filtered.to_csv(output[0], sep='\t', header=False, index=False)

# rule Get_human_Mazin_AS_Segment_PhyloPScores:
#     input:
#         PhyloP = "conservation/PhyloP.hg19.bw"
