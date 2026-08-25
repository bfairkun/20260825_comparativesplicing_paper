# Terminal rules that reduce pipeline intermediates to the small plot-ready tables in
# output/. These are the only inputs the figure notebook reads, so a user who does not want
# to re-run alignment/liftover/clustering can regenerate every figure from the shipped
# output/ directory alone.
#
# All of these are cheap (seconds to a few minutes) given the upstream intermediates.

EXPORT = "scripts/Export_FigureData.py"


def species_of(genome):
    return genome.split("_")[0]

CORDOSO_GENOMES = [
    "Human_UCSC.hg38_GencodeComprehensive46",
    "Macaque_UCSC.rheMac10_ensv101",
    "Mouse_UCSC.mm39_GencodeComprehensive46",
    "Rat_UCSC.rn7_RefSeqv108",
    "Rabbit_UCSC.oryCun2_ensv101",
    "Opossum_UCSC.monDom5_ensvUnknown",
    "Chicken_UCSC.galGal6_ensv101",
]
# species -> the annotation source used for that species in the main figures
NON_HUMAN_GENOMES = [g for g in CORDOSO_GENOMES if not g.startswith("Human_")]

# annotation source x genome combinations that actually exist
EXTRA_GTF_COMBOS = [
    ("Gencode", "Human_UCSC.hg38_GencodeComprehensive46"),
    ("Gencode", "Mouse_UCSC.mm39_GencodeComprehensive46"),
    ("Ensembl", "Chicken_UCSC.galGal6_ensv101"),
    ("Ensembl", "Human_UCSC.hg38_GencodeComprehensive46"),
    ("Ensembl", "Macaque_UCSC.rheMac10_ensv101"),
    ("Ensembl", "Mouse_UCSC.mm39_GencodeComprehensive46"),
    ("Ensembl", "Opossum_UCSC.monDom5_ensvUnknown"),
    ("Ensembl", "Rabbit_UCSC.oryCun2_ensv101"),
    ("RefSeq", "Chicken_UCSC.galGal6_ensv101"),
    ("RefSeq", "Human_UCSC.hg38_GencodeComprehensive46"),
    ("RefSeq", "Macaque_UCSC.rheMac10_ensv101"),
    ("RefSeq", "Mouse_UCSC.mm39_GencodeComprehensive46"),
    ("RefSeq", "Rabbit_UCSC.oryCun2_ensv101"),
    ("RefSeq", "Rat_UCSC.rn7_RefSeqv108"),
]
EXTRA_CLASSIFICATIONS = [
    f"MazinLeafcutterAnalysis/ClassifyJuncs/{src}/{g}.AllObserved._junction_classifications.txt"
    for src, g in EXTRA_GTF_COMBOS
]
EXTRA_TRANSCRIPT_TYPES = [
    f"MazinLeafcutterAnalysis/Reformated_ExtraGTFs/{src}/{g}.TranscriptTypes.tsv.gz"
    for src, g in EXTRA_GTF_COMBOS
]

# SCN8A exons plotted in the conservation panel: two constitutive, one poison exon
SCN8A_EXONS = {"exon17": "chr12:51774189-51774362",
               "exon18n": "chr12:51780202-51780271",
               "exon19": "chr12:51786542-51786826"}


rule Export_NoEmbryoSampleFlags:
    """Which libraries are used for which developmental analyses (read at parse time)."""
    input:
        samples = "config/Cordoso_Moreira_SampleList.tsv",
        stages = "../data/Stages_AsIn_CordosoMoreira_Recoded.txt",
        star_samples = "config/samples.tsv",
    output:
        "../data/Cordoso_SampleList.WhichAreNonHumanEmbryo.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} noembryo-samples --samples {input.samples} --stages {input.stages} "
        "--star-samples {input.star_samples} --out {output}"


rule Export_JuncClassTallies:
    """Fig 3b and the alternate-annotation panels. 1.8 GB + 230 MB -> ~2 KB."""
    input:
        classifications = EXTRA_CLASSIFICATIONS,
        counts = expand("MazinLeafcutterAnalysis/ObservedJuncsAnnotations/{g}.uniq.junc.gz",
                        g=CORDOSO_GENOMES),
    output:
        "../output/DataToPlot.JuncClassTallies.tsv.gz",
    params:
        count_dir = "MazinLeafcutterAnalysis/ObservedJuncsAnnotations",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} junc-class-tallies --classifications {input.classifications} "
        "--count-dir {params.count_dir} --out {output}"


rule Export_TranscriptsPerGene:
    """Fig 3a and the alternate-annotation panel a. 12 MB -> ~2 KB."""
    input:
        transcript_types = EXTRA_TRANSCRIPT_TYPES,
        mouse_fix = ("MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Ensembl/"
                     "Mouse_UCSC.mm39_GencodeComprehensive46.FixedOut.bed"),
    output:
        hist = "../output/DataToPlot.TranscriptsPerGene.tsv.gz",
        totals = "../output/DataToPlot.TranscriptsPerGene.Totals.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} transcripts-per-gene "
        "--transcript-types {input.transcript_types} "
        "--mouse-ensembl-fix {input.mouse_fix} "
        "--out {output.hist} --out-totals {output.totals}"


rule Export_ContrastSummaries:
    """Splicing-vs-expression coupling for every contrast. 2.7 GB / 404 files -> ~13 KB."""
    input:
        dev = expand("MazinLeafcutterAnalysis/Contrasts_ds_tidy/{t}.joined.tsv.gz",
                     t=PAPER_DEV_CONTRASTS),
        adult = expand("MazinLeafcutterAnalysis/ContrastAdultTissues_ds_tidy/{t}.joined.tsv.gz",
                       t=PAPER_ADULT_CONTRASTS),
    output:
        "../output/DataToPlot.ContrastSummaries.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} contrast-summaries --joined {input.dev} {input.adult} --out {output}"


rule Export_BetaBetaAdultBrainKidney:
    """Main Fig 3d scatter points, one adult tissue pair."""
    input:
        expand("MazinLeafcutterAnalysis/ContrastAdultTissues_ds_tidy/"
               "{sp}_Brain-Kidney.joined.tsv.gz", sp=["Chicken", "Human"]),
    output:
        "../output/DataForReplotting_BetaBeta_AdultBrainKidney.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} betabeta-points --joined {input} --out {output}"


rule Export_JuncClassCountsBySample:
    """Per-library junction-productivity read counts (Fig 3c input)."""
    input:
        expand("MazinLeafcutterAnalysis/SummarisedClassificationsBySample/{sp}.tsv.gz",
               sp=[species_of(g) for g in CORDOSO_GENOMES]),
    output:
        "../output/DataToPlot.JuncClassCountsBySample.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} junc-class-by-sample --summaries {input} --out {output}"


rule Export_TestisPercentUnproductive:
    """Percent unproductive reads per library."""
    input:
        "../output/DataToPlot.JuncClassCountsBySample.tsv.gz",
    output:
        "../output/DataToPlot.Testis.Percent.UP.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} testis-percent-up --counts {input} --out {output}"


rule Export_SCN8A_PhyloP:
    """Per-base PhyloP over the three plotted SCN8A exons."""
    input:
        bedgraph = "../data/conservation/SCN8A.PhyloP.bg",
        dnds = "../data/conservation/dnds_exon_trios.dnds_added.tsv.gz",
        mean_phylop = "../data/conservation/dnds_exon_trios.regions.MeanPhyloP.bed",
    output:
        "../output/DataForReplotting_SCN8A_PhyloP.tsv.gz",
    params:
        labels = " ".join(SCN8A_EXONS.keys()),
        regions = " ".join(SCN8A_EXONS.values()),
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} scn8a-phylop --bedgraph {input.bedgraph} --dnds {input.dnds} "
        "--mean-phylop {input.mean_phylop} --exon-labels {params.labels} "
        "--exon-regions {params.regions} --out {output}"


rule Export_dNdS_PhyloP_Trios:
    """dN/dS and mean PhyloP per cassette-exon trio."""
    input:
        dnds = "../data/conservation/dnds_exon_trios.dnds_added.tsv.gz",
        mean_phylop = "../data/conservation/dnds_exon_trios.regions.MeanPhyloP.bed",
        gene_names = "../data/Ensembl.TranscriptInfo.tsv.gz",
    output:
        "../output/DataToPlot.dNdS_PhyloP_ExonTrios.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} dnds-phylop-trios --dnds {input.dnds} "
        "--mean-phylop {input.mean_phylop} --gene-names {input.gene_names} --out {output}"


rule Export_ARHGAP17_FlankingPhyloP:
    """Mean PhyloP over the 50 bp flanking each ARHGAP17 exon."""
    input:
        "../data/conservation/ARHGAP17.exonflanks.MeanPhyloP.bed",
    output:
        "../output/DataForReplotting_ARHGAP17_FlankingPhyloP.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} arhgap17-flanks --bed {input} --out {output}"


rule Export_ConservedEventCascade:
    """Conserved-event counts after each successive filter. ~2.3 GB -> 70 numbers."""
    input:
        lifted = expand("LiftoverJuncs/AsFlanks/{g}.Lifted.bed.gz", g=NON_HUMAN_GENOMES),
        classifications = expand(
            "MazinLeafcutterAnalysis/ClassifyJuncs/{g}.AllObserved._junction_classifications.txt",
            g=CORDOSO_GENOMES),
        clustered = expand(
            "MazinLeafcutterAnalysis/ClassifyJuncs/{g}.Clustered.juncList.tsv.gz",
            g=CORDOSO_GENOMES),
        observed = expand(
            "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/{g}.uniq.annotated.tsv.gz",
            g=CORDOSO_GENOMES),
        joined = "scratch/Dat.DifferentialSplicing.joinedByHumanCoords.tsv.gz",
    output:
        "../output/DataToPlot.ConservedEventCascade.tsv.gz",
    params:
        observed_dir = "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} conserved-event-cascade --lifted {input.lifted} "
        "--classifications {input.classifications} --clustered {input.clustered} "
        "--observed-dir {params.observed_dir} --joined-by-human {input.joined} "
        "--out {output}"


rule Export_NMDFactorExpression:
    """Core NMD-factor expression deltas. 119 MB + 12 MB -> ~155 KB."""
    input:
        expression = expand("MazinLeafcutterAnalysis/Expression/{g}.log2rpkm.tsv.gz",
                            g=CORDOSO_GENOMES),
        orthologs = "../data/Ensembl.GeneHumanHomologs.tsv.gz",
        transcript_info = "../data/Ensembl.TranscriptInfo.tsv.gz",
        sample_stages = "../data/Cordoso_SampleList.WhichAreNonHumanEmbryo.tsv.gz",
    output:
        "../output/DataForReplotting_NMDFactorExpression.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} nmd-expression --expression {input.expression} "
        "--orthologs {input.orthologs} --transcript-info {input.transcript_info} "
        "--sample-stages {input.sample_stages} --out {output}"


rule Export_DevelopmentalMarkers:
    """Percent-rank position of birth / onset-of-meiosis landmarks, per species."""
    input:
        sample_stages = "../data/Cordoso_SampleList.WhichAreNonHumanEmbryo.tsv.gz",
        stages = "../data/Stages_AsIn_CordosoMoreira_Recoded.txt",
    output:
        "../output/DataToPlot.DevelopmentalMarkers.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} developmental-markers --sample-stages {input.sample_stages} "
        "--stages {input.stages} --out {output}"


rule Export_ConservedAS_Examples:
    """PSI, cluster depth and host-gene expression at the plotted example loci.

    Slices ~1.5 GB of PSI tables down to a few hundred rows.
    """
    input:
        loci = "config/ExampleJuncsToPlot.tsv",
        sample_groups = "config/ExampleJuncSampleGroups.tsv",
        psi = expand("rna-seq/SplicingAnalysis/leafcutter/{g}/juncTableBeds/{m}.sorted.bed.gz",
                     g=CORDOSO_GENOMES, m=["PSI", "PSIDenom"]),
        expression = expand("MazinLeafcutterAnalysis/Expression/{g}.log2rpkm.tsv.gz",
                            g=CORDOSO_GENOMES),
    output:
        "../output/DataForReplotting_ConservedAS_Examples.tsv.gz",
    params:
        psi_dir = "rna-seq/SplicingAnalysis/leafcutter",
        expression_dir = "MazinLeafcutterAnalysis/Expression",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} conserved-as-examples --loci {input.loci} "
        "--sample-groups {input.sample_groups} --psi-dir {params.psi_dir} "
        "--expression-dir {params.expression_dir} --out {output}"


rule Export_ConservedJuncsPastMacaque:
    """Supplementary table of conserved junctions shared with a non-macaque species."""
    input:
        lifted = expand("LiftoverJuncs/AsFlanks/{g}.Lifted.bed.gz", g=NON_HUMAN_GENOMES),
        classifications = expand(
            "MazinLeafcutterAnalysis/ClassifyJuncs/{g}.AllObserved._junction_classifications.txt",
            g=CORDOSO_GENOMES),
        clustered = expand(
            "MazinLeafcutterAnalysis/ClassifyJuncs/{g}.Clustered.juncList.tsv.gz",
            g=CORDOSO_GENOMES),
        observed = expand(
            "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/{g}.uniq.annotated.tsv.gz",
            g=CORDOSO_GENOMES),
        joined = "scratch/Dat.DifferentialSplicing.joinedByHumanCoords.tsv.gz",
    output:
        "../output/ConservedJuncs_PastMacaque.tsv.gz",
    params:
        observed_dir = "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} conserved-juncs-past-macaque --lifted {input.lifted} "
        "--classifications {input.classifications} --clustered {input.clustered} "
        "--observed-dir {params.observed_dir} --joined-by-human {input.joined} --out {output}"


rule Export_CassetteExonSupplement:
    """Supplementary table of conserved-but-unconstrained unproductive cassette exons."""
    input:
        dnds = "../data/conservation/dnds_exon_trios.dnds_added.tsv.gz",
        mean_phylop = "../data/conservation/dnds_exon_trios.regions.MeanPhyloP.bed",
        gene_names = "../data/Ensembl.TranscriptInfo.tsv.gz",
    output:
        "../output/ConservedJuncs_PastMacaque.CassetteExons.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} cassette-exon-supplement --dnds {input.dnds} "
        "--mean-phylop {input.mean_phylop} --gene-names {input.gene_names} --out {output}"


rule all_supplement_tables:
    """Published supplementary data files."""
    input:
        "../output/ConservedJuncs_PastMacaque.tsv.gz",
        "../output/ConservedJuncs_PastMacaque.CassetteExons.tsv.gz",


rule all_plot_data:
    """Everything the figure notebook reads."""
    input:
        "../output/DataToPlot.JuncClassTallies.tsv.gz",
        "../output/DataToPlot.TranscriptsPerGene.tsv.gz",
        "../output/DataToPlot.TranscriptsPerGene.Totals.tsv.gz",
        "../output/DataToPlot.ContrastSummaries.tsv.gz",
        "../output/DataForReplotting_BetaBeta_AdultBrainKidney.tsv.gz",
        "../output/DataToPlot.JuncClassCountsBySample.tsv.gz",
        "../output/DataToPlot.Testis.Percent.UP.tsv.gz",
        "../output/DataForReplotting_SCN8A_PhyloP.tsv.gz",
        "../output/DataToPlot.dNdS_PhyloP_ExonTrios.tsv.gz",
        "../output/DataForReplotting_ARHGAP17_FlankingPhyloP.tsv.gz",
        "../output/DataToPlot.ConservedEventCascade.tsv.gz",
        "../output/DataForReplotting_NMDFactorExpression.tsv.gz",
        "../output/DataToPlot.DevelopmentalMarkers.tsv.gz",
        "../output/DataForReplotting_ConservedAS_Examples.tsv.gz",


rule RefreshEnsemblAnnotations:
    """Re-fetch gene/transcript biotypes and human orthologs from Ensembl via biomaRt.

    NOT part of `all`. The shipped copies in data/ are what the workflow reads, because
    biomaRt talks to a live service and the UCSC tracks these must match are frozen
    snapshots (Ensembl v110 for mm39). Run this to refresh them, diff against data/, and
    copy over deliberately:

        snakemake all_refresh_annotations
        diff <(zcat scratch/Ensembl.TranscriptInfo.tsv.gz) \
             <(zcat ../data/Ensembl.TranscriptInfo.tsv.gz)
    """
    output:
        transcripts = "scratch/Ensembl.TranscriptInfo.tsv.gz",
        homologs = "scratch/Ensembl.GeneHumanHomologs.tsv.gz",
    conda: "../envs/r_essentials.yml"
    shell:
        "Rscript scripts/FetchEnsemblAnnotations.R {output.transcripts} {output.homologs}"


rule all_refresh_annotations:
    """Optional: re-fetch the Ensembl annotation tables shipped in data/."""
    input:
        rules.RefreshEnsemblAnnotations.output


rule FixMouseEnsemblBed:
    """Repair UCSC's mouse Ensembl BED, which stores transcript ids in the gene_id field.

    Closes a gap in the original pipeline: this correction was made once in a notebook
    (analysis/2024-09-23_FixMouseEnsemblBed12ToGtf.Rmd) and its output was then read by
    downstream rules without ever being produced by one.
    """
    input:
        bed = ("MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Ensembl/"
               "Mouse_UCSC.mm39_GencodeComprehensive46.bed"),
        info = "../data/Ensembl.TranscriptInfo.tsv.gz",
    output:
        ("MazinLeafcutterAnalysis/Reformated_ExtraGTFs/Ensembl/"
         "Mouse_UCSC.mm39_GencodeComprehensive46.Fixed.bed"),
    conda: "../envs/py_general.yml"
    shell:
        "python scripts/FixMouseEnsemblBed.py --bed {input.bed} "
        "--transcript-info {input.info} --out {output}"


rule JoinDifferentialSplicingByHumanCoords:
    """Place each species' differential-splicing results beside human, by lifted coordinate.

    Closes another gap in the original pipeline: this join was built in a notebook and
    exported by hand. Only the columns the cascade needs are kept, so the intermediate stays
    a few tens of MB.
    """
    input:
        joined = (expand("MazinLeafcutterAnalysis/Contrasts_ds_tidy/{t}.joined.tsv.gz",
                         t=PAPER_DEV_CONTRASTS)
                  + expand("MazinLeafcutterAnalysis/ContrastAdultTissues_ds_tidy/{t}.joined.tsv.gz",
                           t=PAPER_ADULT_CONTRASTS)),
        lifted = expand("LiftoverJuncs/AsFlanks/{g}.Lifted.bed.gz", g=NON_HUMAN_GENOMES),
    output:
        "scratch/Dat.DifferentialSplicing.joinedByHumanCoords.tsv.gz",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} join-by-human-coords --joined {input.joined} "
        "--lifted {input.lifted} --out {output}"


rule MakeContrastGroupFiles:
    """Derive the leafcutter group files from the supplementary contrast table.

    The table lists every contrast the paper used and the libraries in each arm, so it is the
    single authority; previously these 223 files were written once by a notebook and read by
    downstream rules without any rule producing them.
    """
    input:
        table = CONTRAST_TABLE,
    output:
        dev = expand("config/ContrastGroupFiles/{t}.txt", t=PAPER_DEV_CONTRASTS),
        adult = expand("config/ContrastAdultTissueGroupFiles/{t}.txt",
                       t=PAPER_ADULT_CONTRASTS),
    params:
        dev_dir = "config/ContrastGroupFiles",
        adult_dir = "config/ContrastAdultTissueGroupFiles",
    conda: "../envs/py_general.yml"
    shell:
        "python {EXPORT} contrast-group-files --table {input.table} "
        "--dev-dir {params.dev_dir} --adult-dir {params.adult_dir}"
