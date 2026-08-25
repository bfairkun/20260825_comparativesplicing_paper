# An example collection of Snakemake rules imported in the main Snakefile.
rule RemoveTrailingWhitespaceFromGtf:
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
    output:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.NoTrailingWhitespace.gtf",
    shell:
        """
        sed  -e 's/[[:space:]]*$//' {input} > {output}
        """

use rule rna_seq_annotate_juncfiles as annotate_juncfiles_NoWhitespaceGtf with:
    input:
        fa = rna_seq.FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa"),
        fai = rna_seq.FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai"),
        gtf = rna_seq.FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.NoTrailingWhitespace.gtf"),
        juncs = "rna-seq/SplicingAnalysis/juncfiles/{sample}.junc",
    wildcard_constraints:
        sample = "|".join(samples.loc[samples['STARGenomeName'] == "StarletSeaAnemone_RefSeq_GCF_932526225.1"].index)
    output:
        counts = "rna-seq/SplicingAnalysis/juncfiles/{sample}.junccounts.tsv.gz"
    conda:
        "../module_workflows/rna_seq/envs/regtools.yml"

use rule rna_seq_AnnotateConcatedUniqJuncFile_basic as AnnotateConcatedUniqJuncFile_basic_NoWhitespaceGtf with:
    input:
        junc = "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.NoTrailingWhitespace.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa"
    wildcard_constraints:
        GenomeName = "|".join(["StarletSeaAnemone_RefSeq_GCF_932526225.1"])
    output:
        "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz"
    conda:
        "../module_workflows/rna_seq/envs/regtools.yml"
