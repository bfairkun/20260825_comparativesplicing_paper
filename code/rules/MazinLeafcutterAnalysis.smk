rule ConvertLeacfutterJuncNames_To_True_Bed_Coords:
    input:
        Clustered_juncs = "rna-seq/SplicingAnalysis/leafcutter/{GenomeName}/clustering/leafcutter_perind.counts.gz",
        All_Observed_juncs = "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz"
    output:
        Clustered = "MazinLeafcutterAnalysis/ClassifyJuncs/{GenomeName}.Clustered.juncList.tsv.gz",
        AllObserved = "MazinLeafcutterAnalysis/ClassifyJuncs/{GenomeName}.AllObserved.juncList.tsv.gz"
    log:
        "logs/ConvertLeacfutterJuncNames_To_True_Bed_Coords/{GenomeName}.log"
    conda:
        "../envs/r_essentials.yml"
    resources:
        mem_mb = 24000
    shell:
        """
        Rscript scripts/PrepAllJuncsFor_JunctionClassifier.R {input.Clustered_juncs} {input.All_Observed_juncs} {output.Clustered} {output.AllObserved} 2> {log}
        """

def GetParamsForReformatGtf(wildcards):
    if wildcards.GenomeName in ["Chicken_UCSC.galGal6_ensv101", "Opossum_UCSC.monDom5_ensvUnknown", "Macaque_UCSC.rheMac10_ensv101", "Rabbit_UCSC.oryCun2_ensv101"]:
        return "-infer_gene_type_approach B -infer_transcript_type_approach B"
    elif wildcards.GenomeName in ["Rat_UCSC.rn7_RefSeqv108"]:
        return "-infer_gene_type_approach B -infer_transcript_type_approach B -transcript_name_attribute_name transcript_id"
    elif wildcards.GenomeName in ["Human_UCSC.hg38_GencodeComprehensive46", "Mouse_UCSC.mm39_GencodeComprehensive46"]:
        # infer_transcript_type_approach A more appropriate since nonsense_mediated_decay transcripts do have CDS child features
        return "-infer_gene_type_approach A -infer_transcript_type_approach A -transcript_name_attribute_name transcript_id -gene_name_attribute_name gene_id"
    else:
        return ""


def GetFastaLinkForGencode(wildcards):
    if wildcards.GenomeName == "Human_UCSC.hg38_GencodeComprehensive46":
        return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
    elif wildcards.GenomeName == "Mouse_UCSC.mm39_GencodeComprehensive46":
        return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.primary_assembly.genome.fa.gz"

rule DownloadHumanAndMouse_GencodeFa:
    output:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.GencodePrimary.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.GencodePrimary.fa.fai",
    params:
        link = GetFastaLinkForGencode
    wildcard_constraints:
        GenomeName = "|".join(["Human_UCSC.hg38_GencodeComprehensive46", "Mouse_UCSC.mm39_GencodeComprehensive46"])
    shell:
        """
        wget -O- {params.link} | gunzip -c > {output.fa}
        samtools faidx {output.fa}
        """

rule GatherGencodeFa:
    input:
        expand(config['GenomesPrefix'] + "{GenomeName}/Reference.GencodePrimary.fa", GenomeName = ["Human_UCSC.hg38_GencodeComprehensive46", "Mouse_UCSC.mm39_GencodeComprehensive46"])


def GetFastaFor_ReformatGTF(wildcards):
    if wildcards.GenomeName in ["Human_UCSC.hg38_GencodeComprehensive46", "Mouse_UCSC.mm39_GencodeComprehensive46"]:
        return config['GenomesPrefix'] + "{GenomeName}/Reference.GencodePrimary.fa"
    else:
        return config['GenomesPrefix'] + "{GenomeName}/Reference.fa",

rule ReformatGTF:
    """
    UCSC mouse and human fasta is not compatible (slightly different chrom
    names for patches/scaffolds with the Gencode gtf for which I am using
    because it nicely has the most tags for transcript_type and gene_type).
    Therefore, input function uses the Gencode fasta for human and mouse
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        fa = GetFastaFor_ReformatGTF
    output:
        gtf = "MazinLeafcutterAnalysis/ReformatedGTFs/{GenomeName}.gtf",
        bed = "MazinLeafcutterAnalysis/ReformatedGTFs/{GenomeName}.bed"
    log:
        "logs/ReformatGTF/{GenomeName}.log"
    params:
        SpeciesSpecificParams = GetParamsForReformatGtf,
        extra = ""
    resources:
        mem_mb = 24000
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        python scripts/leafcutter2_bfairkun/scripts/Reformat_gtf.py -i {input.gtf} -o {output.gtf} -fa {input.fa} -v {params.SpeciesSpecificParams} {params.extra} -bed12_out {output.bed} 2> {log}
        """

rule Gather_Reformatted_GTF:
    input:
        expand("MazinLeafcutterAnalysis/ReformatedGTFs/{GenomeName}.gtf", GenomeName=CordosMoreiraGenomes)

rule leafcutter2_ClassifyJuncs_ClusterPerInd:
    """
    Need to reformat gtfs to be more gencode like. Attributes for genes is gene_name
    """
    input:
        junclist = "MazinLeafcutterAnalysis/ClassifyJuncs/{GenomeName}.{JunctionSet}.juncList.tsv.gz",
        gtf = "MazinLeafcutterAnalysis/ReformatedGTFs/{GenomeName}.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    output:
        "MazinLeafcutterAnalysis/ClassifyJuncs/{GenomeName}.{JunctionSet}._junction_classifications.txt"
    log:
        "logs/leafcutter2_ClassifyJuncs_ClusterPerInd/{GenomeName}.{JunctionSet}.log"
    params:
        rundir = "MazinLeafcutterAnalysis/ClassifyJuncs"
    wildcard_constraints:
        JunctionSet = "Clustered|AllObserved"
    resources:
        mem_mb = 52000
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        python scripts/leafcutter2_daiuc/scripts/ForwardSpliceJunctionClassifier.py -c {input.junclist} -G {input.fa} -A {input.gtf} -v -r {params.rundir} -o {wildcards.GenomeName}.{wildcards.JunctionSet}. --max_juncs 1000 &> {log}
        """

rule Add_protein_coding_to_all_lines:
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
    output:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.ProteinCodingTagAddedToAll.gtf",
    log:
        "logs/Add_protein_coding_to_all_lines/{GenomeName}.log"
    shadow:
        "shallow"
    conda:
        "../envs/bedparse.yml"
    resources:
        mem_mb = 42000
    shell:
        """
        exec 2> {log}
        agat_convert_sp_gxf2gxf.pl -g {input} -o {wildcards.GenomeName}.gtf
        sed -E 's/(transcript_type|gene_type)=[^;]*;//g; s/$/;transcript_type=protein_coding;gene_type=protein_coding;/' {wildcards.GenomeName}.gtf > {output}
        """

rule SpearmanCoef_Splicing:
    input:
        PSI = "rna-seq/SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/{Metric}.sorted.bed.gz",
        samples = "config/Cordoso_Moreira_SampleList.tsv"
    output:
        "MazinLeafcutterAnalysis/SplicingSpearmanCoefs/{GenomeName}.{Metric}.tsv.gz"
    resources:
        mem_mb = 52000
    log:
        "logs/SpearmanCoef_Splicing/{GenomeName}.{Metric}.log"
    conda:
        "../module_workflows/rna_seq/envs/r_2.yml"
    wildcard_constraints:
        Metric = "PSI|PSI_ByMax"
    shell:
        """
        Rscript scripts/SpearmanCor_Mazin_LeafcutterPSI.R {input.PSI} {input.samples} {output} &> {log}
        """


use rule rna_seq_featurecounts_allUnstranded as featurecounts_allUnstranded_ReformattedGtf with:
    input:
        bam = rna_seq.ExpandAllSamplesInFormatStringFromGenomeNameWildcard("rna-seq/Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam"),
        index = rna_seq.ExpandAllSamplesInFormatStringFromGenomeNameWildcard("rna-seq/Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.indexing_done"),
        gtf = "MazinLeafcutterAnalysis/ReformatedGTFs/{GenomeName}.gtf",
    output:
        counts = "featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt",
        summary = "featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt.summary",
    conda:
        "../module_workflows/rna_seq/envs/subread_featureCounts.yml"
    log:
        "logs/featureCounts/{GenomeName}.AllUnstranded.log"
    params:
        strand = '-s 0',
        extra = "-g gene_name"

rule featureCounts_to_Matrix:
    input:
        counts = "featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt",
    output:
        counts = "MazinLeafcutterAnalysis/Expression/{GenomeName}.counts.tsv.gz",
        logrpkm = "MazinLeafcutterAnalysis/Expression/{GenomeName}.log2rpkm.tsv.gz",
        logtpm = "MazinLeafcutterAnalysis/Expression/{GenomeName}.log2tpm.tsv.gz",
        logcpm = "MazinLeafcutterAnalysis/Expression/{GenomeName}.log2cpm.tsv.gz",
    log:
        "logs/featureCounts_to_Matrix/{GenomeName}.log"
    conda:
        "../envs/r_essentials.yml"
    resources:
        mem_mb = 16000
    shell:
        """
        Rscript scripts/FeatureCounts_to_Mat.R {input} {output.counts} {output.logrpkm} {output.logtpm} {output.logcpm} 2> {log}
        """


rule SpearmanCoef_Expression:
    input:
        logrpkm = "MazinLeafcutterAnalysis/Expression/{GenomeName}.{Metric}.tsv.gz",
        samples = "config/Cordoso_Moreira_SampleList.tsv"
    output:
        Basic = "MazinLeafcutterAnalysis/ExpressionSpearmanCoefs/{GenomeName}.{Metric}.tsv.gz",
        Extended = "MazinLeafcutterAnalysis/ExpressionSpearmanCoefs/{GenomeName}.{Metric}.Extended.tsv.gz",
    resources:
        mem_mb = 52000
    log:
        "logs/SpearmanCoef_Splicing/{GenomeName}.{Metric}.log"
    conda:
        "../module_workflows/rna_seq/envs/r_2.yml"
    wildcard_constraints:
        Metric = "log2rpkm|log2tpm|log2cpm"
    shell:
        """
        Rscript scripts/SpearmanCor_Mazin_log2RPKM.R {input.logrpkm} {input.samples} {output.Basic} {output.Extended} &> {log}
        """

rule GetExtra_GTF_From_link:
    output:
        "MazinLeafcutterAnalysis/GTFs/Ensembl/{GenomeName}.gtf"
    params:
        link = lambda wildcards: Extra_Gtfs.loc[wildcards.GenomeName, 'EnsemblGtf']
    wildcard_constraints:
        GenomeName = "|".join(GenomesWithEnsembl)
    shell:
        """
        wget -O- {params.link} | zcat - > {output}
        """

use rule GetExtra_GTF_From_link as GetExtra_GTF_From_link_RefSeq with:
    output:
        "MazinLeafcutterAnalysis/GTFs/RefSeq/{GenomeName}.gtf"
    params:
        link = lambda wildcards: Extra_Gtfs.loc[wildcards.GenomeName, 'RefSeqGtf']
    wildcard_constraints:
        GenomeName = "|".join(GenomesWithRefSeq)

# def GetParamsForReformatGtf_ByStyle(wildcards):
#     if wildcards.GtfStyle == "RefSeq":
#         return "-infer_gene_type_approach B -infer_transcript_type_approach B -transcript_name_attribute_name transcript_id"
#     elif wildcards.GtfStyle == "Ensembl":
#         return "-infer_gene_type_approach B -infer_transcript_type_approach B -transcript_name_attribute_name transcript_id"

use rule ReformatGTF as ReformatGTF_ExtraGTF_RefSeq with:
    input:
        gtf = "MazinLeafcutterAnalysis/GTFs/{GtfStyle}/{GenomeName}.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
    output:
        gtf = "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/{GtfStyle}/{GenomeName}.gtf",
        bed = "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/{GtfStyle}/{GenomeName}.bed",
    log:
        "logs/ReformatGTF_ExtraGTF/{GtfStyle}/{GenomeName}.log"
    params:
        SpeciesSpecificParams = "",
        extra = "-infer_gene_type_approach B -infer_transcript_type_approach B -transcript_name_attribute_name transcript_id -gene_name_attribute_name gene_id" 
    wildcard_constraints:
        GtfStyle = "RefSeq|Ensembl"


def Get_Gtf_For_ClassifyJuncs_ExtraGtfs(wildcards):
    if wildcards.GtfStyle == "Gencode":
        return "MazinLeafcutterAnalysis/ReformatedGTFs/{GenomeName}.gtf"
    elif wildcards.GtfStyle == "Ensembl" and wildcards.GenomeName == "Mouse_UCSC.mm39_GencodeComprehensive46":
        return "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Ensembl/Mouse_UCSC.mm39_GencodeComprehensive46.FixedOut.gtf"
    else:
        return "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/{GtfStyle}/{GenomeName}.gtf"

use rule leafcutter2_ClassifyJuncs_ClusterPerInd as leafcutter2_ClassifyJuncs_ClusterPerInd_NewGtfs with:
    input:
        junclist = "MazinLeafcutterAnalysis/ClassifyJuncs/{GenomeName}.{JunctionSet}.juncList.tsv.gz",
        gtf = Get_Gtf_For_ClassifyJuncs_ExtraGtfs,
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    output:
        "MazinLeafcutterAnalysis/ClassifyJuncs/{GtfStyle}/{GenomeName}.{JunctionSet}._junction_classifications.txt"
    log:
        "logs/leafcutter2_ClassifyJuncs_ClusterPerInd_NewGtfs/{GtfStyle}.{GenomeName}.{JunctionSet}.log"
    params:
        rundir = "MazinLeafcutterAnalysis/ClassifyJuncs/{GtfStyle}"
    wildcard_constraints:
        JunctionSet = "Clustered|AllObserved"
    # resources:
    #     mem_mb = 48000
    # conda:
    #     "../envs/bedparse.yml"
    # shell:
    #     """
    #     python scripts/leafcutter2_daiuc/scripts/ForwardSpliceJunctionClassifier.py -c {input.junclist} -G {input.fa} -A {input.gtf} -v -r {params.rundir} -o {wildcards.GenomeName}.{wildcards.JunctionSet}. &> {log}
    #     """

rule Gather_New_Gtfs:
    input:
        expand("MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Ensembl/{GenomeName}.TranscriptTypes.tsv.gz", GenomeName = GenomesWithEnsembl.unique()),
        expand("MazinLeafcutterAnalysis/Reformated_ExtraGTFs/RefSeq/{GenomeName}.TranscriptTypes.tsv.gz", GenomeName = GenomesWithRefSeq.unique()),
        expand("MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Gencode/{GenomeName}.TranscriptTypes.tsv.gz", GenomeName = ['Human_UCSC.hg38_GencodeComprehensive46','Mouse_UCSC.mm39_GencodeComprehensive46']),
        expand("MazinLeafcutterAnalysis/ClassifyJuncs/Ensembl/{GenomeName}.AllObserved._junction_classifications.txt", GenomeName = GenomesWithEnsembl.unique()),
        expand("MazinLeafcutterAnalysis/ClassifyJuncs/RefSeq/{GenomeName}.AllObserved._junction_classifications.txt", GenomeName = GenomesWithRefSeq.unique()),
        expand("MazinLeafcutterAnalysis/ClassifyJuncs/Gencode/{GenomeName}.AllObserved._junction_classifications.txt", GenomeName = ['Human_UCSC.hg38_GencodeComprehensive46','Mouse_UCSC.mm39_GencodeComprehensive46']),

rule TranscriptTypes_ForGtf:
    input:
        gtf = "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/{GtfStyle}/{GenomeName}.gtf",
    output:
        "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/{GtfStyle}/{GenomeName}.TranscriptTypes.tsv.gz",
    shell:
        """
        awk -v OFS='\\t' 'BEGIN {{print "gname", "gtype", "tname", "ttype"}} $3 == "transcript" {{match($0, /gene_name "([^"]+)"/, gname); match($0, /gene_type "([^"]+)"/, gtype); match($0, /transcript_type "([^"]+)"/, ttype); match($0, /transcript_name "([^"]+)"/, tname); print gname[1], gtype[1], tname[1], ttype[1]}}' {input.gtf} | gzip - > {output}
        """

use rule TranscriptTypes_ForGtf as TranscriptTypes_ForGtf_Gencode with:
    input:
        gtf = "MazinLeafcutterAnalysis/ReformatedGTFs/{GenomeName}.gtf",
    output:
        "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Gencode/{GenomeName}.TranscriptTypes.tsv.gz",
    wildcard_constraints:
        GenomeName = "Human_UCSC.hg38_GencodeComprehensive46|Mouse_UCSC.mm39_GencodeComprehensive46"


rule Gather_Annotated_beds:
    input:
        expand("MazinLeafcutterAnalysis/ReformatedGTFs/{GenomeName}.bed", GenomeName = CordosMoreiraGenomes)
rule Fix_Mouse_UCSC_Ensembl_gtf:
    input:
        gtf = "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Ensembl/Mouse_UCSC.mm39_GencodeComprehensive46.Fixed.bed",
        fa = config['GenomesPrefix'] + "Mouse_UCSC.mm39_GencodeComprehensive46/Reference.fa",
    output:
        gtf = "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Ensembl/Mouse_UCSC.mm39_GencodeComprehensive46.FixedOut.gtf",
        bed = "MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Ensembl/Mouse_UCSC.mm39_GencodeComprehensive46.FixedOut.bed",
    log:
        "logs/Fix_Mouse_UCSC_Ensembl_gtf.log"
    params:
        SpeciesSpecificParams = "",
        extra = "-infer_gene_type_approach B -infer_transcript_type_approach B -transcript_name_attribute_name transcript_id -gene_name_attribute_name gene_id" 
    params:
        SpeciesSpecificParams = GetParamsForReformatGtf,
        extra = ""
    resources:
        mem_mb = 24000
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        python scripts/leafcutter2_bfairkun/scripts/Reformat_gtf.py -input_type bed12 -i {input.gtf} -o {output.gtf} -fa {input.fa} -v {params.SpeciesSpecificParams} {params.extra} -bed12_out {output.bed} 2> {log}
        """

rule QQ_norm_splicing:
    input:
        PSI_bed = "rna-seq/SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/PSI.sorted.bed.gz",
        sample_table = "config/Cordoso_Moreira_SampleList.tsv"
    output:
        stats = "MazinLeafcutterAnalysis/Splicing_qqNormed/{GenomeName}/DevAS.tsv.gz",
        PSI_tables = directory("MazinLeafcutterAnalysis/Splicing_qqNormed/{GenomeName}/PSI_tables")
    resources:
        mem_mb = 52000
    log:
        "logs/QQ_norm_splicing/{GenomeName}.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        mkdir {output.PSI_tables}
        Rscript scripts/QQNorm_PSITable_ByTissue.R {input.PSI_bed} {input.sample_table} {output.PSI_tables}/ {output.stats} &> {log}
        """

rule leafcutter_ds_contrasts:
    input:
        groupfile = "config/ContrastGroupFiles/{treatment}.txt",
        numers = lambda wildcards: f"rna-seq/SplicingAnalysis/leafcutter/{Cordoso_contrasts.loc[wildcards.treatment]['GenomeName']}/clustering/leafcutter_perind_numers.counts.gz",
        wrapper = "scripts/leafcutter_ds_py.py"
    output:
        "MazinLeafcutterAnalysis/Contrasts_ds/{treatment}_effect_sizes.txt",
        "MazinLeafcutterAnalysis/Contrasts_ds/{treatment}_cluster_significance.txt"
    threads: 4
    wildcard_constraints:
        treatment = "|".join(PAPER_DEV_CONTRASTS)
    resources:
        ntasks = 5,
        mem_mb = 24000
    params:
        Prefix = "MazinLeafcutterAnalysis/Contrasts_ds/",
        ExtraParams = "-i 2 -g 2"
    conda:
        "../envs/leafcutter_ds_py.yml"
    log:
        "logs/leafcutter_ds/{treatment}.log"
    shell:
        """
        python {input.wrapper} -p {threads} -o {params.Prefix}{wildcards.treatment} {params.ExtraParams} {input.numers} {input.groupfile} &> {log}
        """

use rule leafcutter_ds_contrasts as leafcutter_ds_contrasts_AdultTissues with:
    input:
        groupfile = "config/ContrastAdultTissueGroupFiles/{treatment}.txt",
        numers = lambda wildcards: f"rna-seq/SplicingAnalysis/leafcutter/{Cordoso_TissueContrasts.loc[wildcards.treatment]['GenomeName']}/clustering/leafcutter_perind_numers.counts.gz",
        wrapper = "scripts/leafcutter_ds_py.py"
    wildcard_constraints:
        treatment = "|".join(PAPER_ADULT_CONTRASTS)
    params:
        Prefix = "MazinLeafcutterAnalysis/ContrastAdultTissues_ds/",
        ExtraParams = "-i 2 -g 2"
    output:
        "MazinLeafcutterAnalysis/ContrastAdultTissues_ds/{treatment}_effect_sizes.txt",
        "MazinLeafcutterAnalysis/ContrastAdultTissues_ds/{treatment}_cluster_significance.txt"
    log:
        "logs/leafcutter_ds_AdultTissues/{treatment}.log"

rule Gather_AdultTissueContrasts:
    input:
        expand("MazinLeafcutterAnalysis/ContrastAdultTissues_ds_tidy/{treatment}.joined.tsv.gz", treatment = PAPER_ADULT_CONTRASTS)

rule tidy_leafcutter_ds_and_de:
    input:
        groupfile = "config/ContrastGroupFiles/{treatment}.txt",
        classifications = lambda wildcards: f"MazinLeafcutterAnalysis/ClassifyJuncs/{Cordoso_contrasts.loc[wildcards.treatment]['GenomeName']}.AllObserved._junction_classifications.txt",
        ClusterSig = "MazinLeafcutterAnalysis/Contrasts_ds/{treatment}_cluster_significance.txt",
        EffectSizes = "MazinLeafcutterAnalysis/Contrasts_ds/{treatment}_effect_sizes.txt",
        Counts = lambda wildcards: f"MazinLeafcutterAnalysis/Expression/{Cordoso_contrasts.loc[wildcards.treatment]['GenomeName']}.counts.tsv.gz",
    output:
        de = "MazinLeafcutterAnalysis/Contrasts_ds_tidy/{treatment}.de.tsv.gz",
        tidy = "MazinLeafcutterAnalysis/Contrasts_ds_tidy/{treatment}.joined.tsv.gz",
    resources:
        mem_mb = 24000
    wildcard_constraints:
        treatment = "|".join(PAPER_DEV_CONTRASTS)
    log:
        "logs/tidy_leafcutter_ds_and_de/{treatment}.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/tidy_leafcutter_ds_results.R {input.groupfile} {input.classifications} {input.ClusterSig} {input.EffectSizes} {input.Counts} {output.de} {output.tidy} &> {log}
        """

use rule tidy_leafcutter_ds_and_de as tidy_leafcutter_ds_and_de_AdultTissues with:
    input:
        groupfile = "config/ContrastAdultTissueGroupFiles/{treatment}.txt",
        classifications = lambda wildcards: f"MazinLeafcutterAnalysis/ClassifyJuncs/{Cordoso_TissueContrasts.loc[wildcards.treatment]['GenomeName']}.AllObserved._junction_classifications.txt",
        ClusterSig = "MazinLeafcutterAnalysis/ContrastAdultTissues_ds/{treatment}_cluster_significance.txt",
        EffectSizes = "MazinLeafcutterAnalysis/ContrastAdultTissues_ds/{treatment}_effect_sizes.txt",
        Counts = lambda wildcards: f"MazinLeafcutterAnalysis/Expression/{Cordoso_TissueContrasts.loc[wildcards.treatment]['GenomeName']}.counts.tsv.gz",
    output:
        de = "MazinLeafcutterAnalysis/ContrastAdultTissues_ds_tidy/{treatment}.de.tsv.gz",
        tidy = "MazinLeafcutterAnalysis/ContrastAdultTissues_ds_tidy/{treatment}.joined.tsv.gz",
    wildcard_constraints:
        treatment = "|".join(PAPER_ADULT_CONTRASTS)
    log:
        "logs/tidy_leafcutter_ds_and_de_AdultTissues/{treatment}.log"

rule SummariseJuncClassificationsBySample:
    input:
        juncfiles = lambda wildcards: expand("rna-seq/SplicingAnalysis/juncfiles/{sample}.junccounts.tsv.gz", sample=CordosoMoreira_df.loc[(CordosoMoreira_df['ID_Species']== wildcards.Species) & (CordosoMoreira_df['JuncFileGroupForParralelization']== int(wildcards.Chunk))].index),
        Annotations = lambda wildcards: f"MazinLeafcutterAnalysis/ClassifyJuncs/{CordosMoreiraGenomes_dict[wildcards.Species]}.AllObserved._junction_classifications.txt"
    output:
        "MazinLeafcutterAnalysis/SummarisedClassificationsBySample/{Species}_chunks/{Chunk}.tsv.gz"
    log:
        "logs/SummariseJuncClassificationsBySample/{Species}/{Chunk}.log"
    resources:
        mem_mb = 24000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/SummariseJuncFilesByProductivityCounts.R {input.Annotations} {output} {input.juncfiles} &> {log}
        """

rule Concat_SummarisedJuncClassificationChunks:
    input:
        expand("MazinLeafcutterAnalysis/SummarisedClassificationsBySample/{{Species}}_chunks/{Chunk}.tsv.gz", Chunk=[i for i in range(1,11)])
    output:
        "MazinLeafcutterAnalysis/SummarisedClassificationsBySample/{Species}.tsv.gz"
    shell:
        """
        cat {input} > {output}
        """


rule Gather_leafcutter_ds_contrasts:
    input:
        expand("MazinLeafcutterAnalysis/Contrasts_ds_tidy/{treatment}.joined.tsv.gz", treatment = PAPER_DEV_CONTRASTS),
        # expand("MazinLeafcutterAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc.gz", GenomeName = CordosMoreiraGenomes)

rule Gather_PSI_beds:
    input:
        expand("rna-seq/SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/PSI_ByMax.sorted.bed.gz.indexing_done", GenomeName=CordosMoreiraGenomes),


use rule rna_seq_ConcatJuncFilesAndKeepUniq as ConcatJuncFilesAndKeepUniq_NoEmbryo with:
    input:
        ExpandAllSamplesInFormatStringFromGenomeNameWildcard_NoEmbryo("rna-seq/SplicingAnalysis/juncfiles/{sample}.junc"),
    output:
        "MazinLeafcutterAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc.gz"
    log:
        "logs/ConcatJuncFilesAndKeepUniq_NoEmbryo/{GenomeName}.log"
    resources:
        mem_mb = rna_seq.GetMemForSuccessiveAttempts(24000, 48000)

# rule ARHGAP17_FlankingConservation:
#     input:
#         bw = "/project2/yangili1/bjf79/20211209_JingxinRNAseq/code/PhyloP/SourceTrack/PhyloP.hg38.bw",
#         bed = "MazinLeafcutterAnalysis/ReformatedGTFs/Human_UCSC.hg38_GencodeComprehensive46.bed"
#     output:
#         "scratch/ARHGAP17.CassetteExon.flank.MeanPhyloP.bed"
#     shell:
#         """
#         python scripts/GetMeanBigwigOverBedIntervals.py -b {input.bed}
#         """

