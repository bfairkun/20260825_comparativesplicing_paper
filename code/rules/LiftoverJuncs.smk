
rule Download_Chain_to_hg38:
    output:
        "ChainFiles/ToHg38/{GenomeName}.chain.gz"
    params:
        link = lambda wildcards: ChainFiles.loc[wildcards.GenomeName]['Genome_to_Hg38_Chain_Link']
    shell:
        """
        wget -O {output} {params.link}
        """

rule Download_Chain_from_hg38:
    output:
        "ChainFiles/FromHg38/{GenomeName}.chain.gz"
    params:
        link = lambda wildcards: ChainFiles.loc[wildcards.GenomeName]['Hg38_to_Genome_Chain_Link']
    shell:
        """
        wget -O {output} {params.link}
        """

rule Gather_Chains:
    input:
        expand("ChainFiles/FromHg38/{GenomeName}.chain.gz", GenomeName = ChainFiles.index),
        expand("ChainFiles/ToHg38/{GenomeName}.chain.gz", GenomeName = ChainFiles.index),

rule LiftoverJuncs:
    """
    Liftover to hg38 and also back to original species... Dylan found that liftover and crossmap are not equivalent... Liftover seems to do better, more results.
    """
    input:
        ChainToHg38 = "ChainFiles/ToHg38/{GenomeName}.chain.gz",
        Juncs = "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz"
    output:
        LiftedOver = "LiftoverJuncs/AsJuncs/{GenomeName}.Lifted.bed.gz",
        index = touch("LiftoverJuncs/AsJuncs/{GenomeName}.Lifted.bed.gz.indexing_done"),
    log:
        "logs/LiftoverJuncsToHuman/{GenomeName}.log"
    params: GetIndexingParams
    conda:
        "../envs/crossmap.yml"
    resources:
        mem_mb = rna_seq.GetMemForSuccessiveAttempts(24000, 54000)
    shell:
        """
        exec 2> {log}
        set +o pipefail;
        zcat {input.Juncs} | awk -v OFS='\\t' 'NR>1 {{ print $1,$2,$3-1,$4,$5,$6 }}' | CrossMap.py bed --unmap-file /dev/null {input.ChainToHg38} /dev/stdin /dev/stdout | sort -k 1,1 -k2,2n | bgzip /dev/stdin -c > {output.LiftedOver}
        tabix -p bed {output.LiftedOver} && touch {output.index}
        """

rule LiftoverFlankingBlocks:
    """
    Liftover juncs as flanking blocks bed12 to hg38 and also back to original species
    """
    input:
        ChainToHg38 = "ChainFiles/ToHg38/{GenomeName}.chain.gz",
        Juncs = "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz"
    output:
        LiftedOver = "LiftoverJuncs/AsFlanks/{GenomeName}.Crossmapped.bed.gz",
        index = touch("LiftoverJuncs/AsFlanks/{GenomeName}.Crossmapped.bed.gz.indexing_done"),
    log:
        "logs/CrossmappedJuncFlankingSites/{GenomeName}.log"
    params: GetIndexingParams
    conda:
        "../envs/crossmap.yml"
    resources:
        mem_mb = rna_seq.GetMemForSuccessiveAttempts(24000, 54000)
    shell:
        """
        exec 2> {log}
        set +o pipefail;
        zcat {input.Juncs} | awk -v OFS='\\t' 'NR>1 {{ $3=$3-1; print $1,$2-1,$3+1,$4,$5,$6,$2-1,$3+1,"255,0,0",2,"1,1", "0,"$3-$2+1 }}' | CrossMap.py bed --unmap-file /dev/null {input.ChainToHg38} /dev/stdin /dev/stdout | sort -k 1,1 -k2,2n | bgzip /dev/stdin -c > {output.LiftedOver}
        tabix -p bed {output.LiftedOver} && touch {output.index}
        """

rule LiftoverFlankingBlocks_LiftOver:
    """
    Liftover juncs as flanking blocks bed12 to hg38 and also back to original species
    """
    input:
        ChainToHg38 = "ChainFiles/ToHg38/{GenomeName}.chain.gz",
        Juncs = "rna-seq/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz"
    output:
        LiftedOver = "LiftoverJuncs/AsFlanks/{GenomeName}.Lifted.bed.gz",
        index = touch("LiftoverJuncs/AsFlanks/{GenomeName}.Lifted.bed.gz.indexing_done"),
    log:
        "logs/LiftoverJuncFlankingSites/{GenomeName}.log"
    params: GetIndexingParams
    conda:
        "../envs/liftover.yml"
    resources:
        mem_mb = rna_seq.GetMemForSuccessiveAttempts(24000, 54000)
    shell:
        """
        exec 2> {log}
        set +o pipefail;
        zcat {input.Juncs} | awk -v OFS='\\t' 'NR>1 {{ $3=$3-1; print $1,$2-1,$3+1,$4,$5,$6,$2-1,$3+1,"255,0,0",2,"1,1", "0,"$3-$2+1 }}' | liftOver /dev/stdin {input.ChainToHg38} /dev/stdout /dev/null | sort -k 1,1 -k2,2n | bgzip /dev/stdin -c > {output.LiftedOver}
        tabix -p bed {output.LiftedOver} && touch {output.index}
        """

rule Liftover_Back_To_Original_Species_Flanks:
    input:
        ChainFromHg38 = "ChainFiles/FromHg38/{GenomeName}.chain.gz",
        LiftedOver = "LiftoverJuncs/AsFlanks/{GenomeName}.Crossmapped.bed.gz",
    output:
        LiftedBackToOriginal = "LiftoverJuncs/AsFlanks/{GenomeName}.CrossmappedBack.bed.gz",
        index2 = touch("LiftoverJuncs/AsFlanks/{GenomeName}.CrossmappedBack.bed.gz.indexing_done"),
    log:
        "logs/Crossmapped_Back_To_Original_Species_Flanks/{GenomeName}.log"
    params: GetIndexingParams
    conda:
        "../envs/crossmap.yml"
    resources:
        mem_mb = rna_seq.GetMemForSuccessiveAttempts(24000, 54000)
    shell:
        """
        exec 2> {log}
        set +o pipefail;
        CrossMap.py bed --unmap-file /dev/null {input.ChainFromHg38} {input.LiftedOver} /dev/stdout | sort -k 1,1 -k2,2n | bgzip /dev/stdin -c > {output.LiftedBackToOriginal}
        tabix -p bed {params} {output.LiftedBackToOriginal} && touch {output.index2}
        """

rule Liftover_Back_To_Original_Species_Flanks_liftover:
    input:
        ChainFromHg38 = "ChainFiles/FromHg38/{GenomeName}.chain.gz",
        LiftedOver = "LiftoverJuncs/AsFlanks/{GenomeName}.Lifted.bed.gz",
    output:
        LiftedBackToOriginal = "LiftoverJuncs/AsFlanks/{GenomeName}.LiftedBack.bed.gz",
        index2 = touch("LiftoverJuncs/AsFlanks/{GenomeName}.LiftedBack.bed.gz.indexing_done"),
    log:
        "logs/Liftover_Back_To_Original_Species_Flanks/{GenomeName}.log"
    params: GetIndexingParams
    conda:
        "../envs/liftover.yml"
    resources:
        mem_mb = rna_seq.GetMemForSuccessiveAttempts(24000, 54000)
    shell:
        """
        exec 2> {log}
        set +o pipefail;
        liftOver {input.LiftedOver} {input.ChainFromHg38} /dev/stdout /dev/null | sort -k 1,1 -k2,2n | bgzip /dev/stdin -c > {output.LiftedBackToOriginal}
        tabix -p bed {params} {output.LiftedBackToOriginal} && touch {output.index2}
        """

use rule Liftover_Back_To_Original_Species_Flanks as Liftover_Back_To_Original_Species_Juncs with:
    input:
        ChainFromHg38 = "ChainFiles/FromHg38/{GenomeName}.chain.gz",
        LiftedOver = "LiftoverJuncs/AsJuncs/{GenomeName}.Crossmapped.bed.gz",
    output:
        LiftedBackToOriginal = "LiftoverJuncs/AsJuncs/{GenomeName}.CrossmappedBack.bed.gz",
        index2 = touch("LiftoverJuncs/AsJuncs/{GenomeName}.CrossmappedBack.bed.gz.indexing_done"),
    log:
        "logs/Liftover_Back_To_Original_Species_Juncs/{GenomeName}.log"

use rule Liftover_Back_To_Original_Species_Flanks_liftover as Liftover_Back_To_Original_Species_Juncs_liftover with:
    input:
        ChainFromHg38 = "ChainFiles/FromHg38/{GenomeName}.chain.gz",
        LiftedOver = "LiftoverJuncs/AsJuncs/{GenomeName}.Lifted.bed.gz",
    output:
        LiftedBackToOriginal = "LiftoverJuncs/AsJuncs/{GenomeName}.LiftedBack.bed.gz",
        index2 = touch("LiftoverJuncs/AsJuncs/{GenomeName}.LiftedBack.bed.gz.indexing_done"),
    log:
        "logs/Liftover_Back_To_Original_Species_Juncs/{GenomeName}.log"

rule GatherLiftedJuncs:
    input:
        # expand("LiftoverJuncs/AsFlanks/{GenomeName}.LiftedBack.bed.gz", GenomeName = "Macaque_UCSC.rheMac10_ensv101"),
        # expand("LiftoverJuncs/AsJuncs/{GenomeName}.LiftedBack.bed.gz", GenomeName = "Macaque_UCSC.rheMac10_ensv101"),
        expand("LiftoverJuncs/AsFlanks/{GenomeName}.LiftedBack.bed.gz", GenomeName = [g for g in CordosMoreiraGenomes if g != "Human_UCSC.hg38_GencodeComprehensive46"]),
        expand("LiftoverJuncs/AsFlanks/{GenomeName}.Lifted.bed.gz", GenomeName = [g for g in CordosMoreiraGenomes if g != "Human_UCSC.hg38_GencodeComprehensive46"]),
        # expand("LiftoverJuncs/AsFlanks/{GenomeName}.Crossmapped.bed.gz", GenomeName = [g for g in CordosMoreiraGenomes if g != "Human_UCSC.hg38_GencodeComprehensive46"]),
        # expand("LiftoverJuncs/AsFlanks/{GenomeName}.CrossmappedBack.bed.gz", GenomeName = [g for g in CordosMoreiraGenomes if g != "Human_UCSC.hg38_GencodeComprehensive46"]),


        # expand("LiftoverJuncs/AsJuncs/{GenomeName}.LiftedBack.bed.gz", GenomeName = [g for g in CordosMoreiraGenomes if g != "Human_UCSC.hg38_GencodeComprehensive46"]),

# seaborn pandas numpy bioconda::crossmap anaconda::biopython ipykernel conda-forge::jupyterlab bioconda::pybigwig bioconda::pybedtools

# non_stop_decay can have CDS
# gene_biotype can be inferred by presence of transcripts with CDS
# ...but would have to make exceptions by NMDFinderB or by tags in that case
# flag for --infer_transcript_type_approach [UseExisting|NMDetectiveB|ByCDS]
# flag for --infer_gene_type_approach [UseExisting|ContainsCodingChildTranscripts]
# flag --transcript_type_attribute_name
# flag --gene_type_attribute_name

# gene_id and transcript_id must be defined (or the equivalents)
# transcript_biotype can be inferred by NMDFinderB (maybe future)
