#rule DownloadFastqFromLink:
#    """
#    Download w aspera if possible, then try ftp link, then try fasterq-dump
#    """
#    input:
#        aspera_key = config['aspera_key']
#    output:
#        fastq = "CordosoMoreira_Fastq/{sample}.fastq.gz",
#    log:
#        "logs/DownloadFastqFromAsperaLink/{sample}.log"
#    wildcard_constraints:
#    shadow: "shallow"
#    params:
#        aspera_link  = lambda wildcards: CordosoMoreira_df.loc[wildcards.sample]['fastq_aspera']
#    shell:
#        """
#        #Use aspera if aspera key file and aspera links parameters are defined (not empty strings)
#        if [[ ! -z "{input.aspera_key}" && ! -z "{params.aspera_link}" ]]; then
#            for link in {params.aspera_link}
#            do
#                tmpfile=$(mktemp -p . tmp.download.XXXXXXXX.fastq.gz)
#                ascp -v -QT -l 300m -P33001 -i {input.aspera_key} era-fasp@${{link}} $tmpfile &>> {log}
#                cat $tmpfile >> {output.fastq}
#                rm $tmpfile
#            done
#        fi
#        """

# rule Download_CordosoMoreira_globus_batch_MakeBatches:
#     output:
#         expand("Downloads/Cordos_Moreira_DownloadBatches/{n}.txt", n=range(5))
#     run:
#         def process_group(group, mod_value):
#             output_file = f'Downloads/Cordos_Moreira_DownloadBatches/{mod_value}.txt'
#             with open(output_file, 'w') as f:
#                 for idx, row in group.iterrows():
#                     col1 = row['fastq_aspera'].replace('fasp.sra.ebi.ac.uk:', '')
#                     col2 = f"{os.getcwd()}/CordosoMoreira_Fastq/{idx}.fastq.gz"
#                     f.write(f"{col1} {col2}\n")
#         # Group by the RowNumber_mod_5 column and process each group
#         for mod_value, group in CordosoMoreira_df.groupby('RowNumber_mod_5'):
#             process_group(group, mod_value)


rule DownloadFastqFromENA:
    """Fetch one Cardoso-Moreira library from ENA.

    Uses aspera when `aspera_key` is configured, otherwise plain FTP, so the workflow does
    not depend on any one user's aspera install. Both URL columns are already present in
    config/Cordoso_Moreira_SampleList.tsv.
    """
    output:
        FASTQ_DIR + "/{SampleTitle}.fastq.gz",
    params:
        ftp = lambda wildcards: CordosoMoreira_df.loc[wildcards.SampleTitle, "fastq_ftp"],
        aspera = lambda wildcards: str(
            CordosoMoreira_df.loc[wildcards.SampleTitle, "fastq_aspera"]
        ).replace("fasp.sra.ebi.ac.uk:", ""),
        key = config.get("aspera_key", ""),
    wildcard_constraints:
        SampleTitle = "|".join(CordosoMoreira_df.index)
    log:
        "logs/DownloadFastqFromENA/{SampleTitle}.log"
    shell:
        """
        set -euo pipefail
        if [ -n "{params.key}" ] && [ -f "{params.key}" ]; then
            ascp -QT -l 300m -P33001 -i {params.key} \
                era-fasp@fasp.sra.ebi.ac.uk:{params.aspera} {output} &> {log}
        else
            wget -O {output} "ftp://{params.ftp}" &> {log}
        fi
        """


# Kept for the record: the original bulk transfer used globus with hardcoded endpoint UUIDs
# and an empty input, so it could not be driven by the DAG.
#
# rule Download_CordosoMoreira_globus_batch_Transfer:
#     output:
#         expand("CordosoMoreira_Fastq/{SampleTitle}.fastq.gz",
#                SampleTitle=CordosoMoreira_df.index),
#     shell:
#         "globus transfer --batch {input} <endpoint-uuid> <endpoint-uuid>"


rule DownloadFastaAndGtf:
    output:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    params:
        fa_link = lambda wildcards: STAR_genomes.loc[wildcards.GenomeName]['FastaLink'],
        gtf_link = lambda wildcards: STAR_genomes.loc[wildcards.GenomeName]['GtfLink'],
    shell:
        """
        wget -O- {params.fa_link} | zcat > {output.fa}
        wget -O- {params.gtf_link} | zcat > {output.gtf}
        """

rule faidxGenome:
    input:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
    output:
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    shell:
        """
        samtools faidx {input.fa}
        """


rule gtf_to_bed12:
    """
    This shell command works with Gencode formatted gtf, haven't looked at
    Ensembl gtfs
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        bed = config['GenomesPrefix'] + "{GenomeName}/Reference.bed.gz",
        index = touch(config['GenomesPrefix'] + "{GenomeName}/Reference.bed.gz.indexing_done")
    params:
        tabix_params = GetIndexingParams
    log:
        "logs/gtf_to_bed12/{GenomeName}.log"
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        (bedparse gtf2bed {input.gtf} --extraFields gene_id,transcript_id,gene_type,gene_name,transcript_type,transcript_support_level,tag,transcript_name | awk -F'\\t' -v OFS='\\t' '{{$4=$NF; print $0}}' | bedtools sort -i - | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        tabix {params.tabix_params} -f -p bed {output.bed} && touch {output.index}
        """


rule gtftools:
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        exons = temp(config['GenomesPrefix'] + "{GenomeName}/GTFTools/Exons.bed"),
        # genes = temp(config['GenomesPrefix'] + "{GenomeName}/GTFTools/Genes.bed"),
        introns = temp(config['GenomesPrefix'] + "{GenomeName}/GTFTools/Introns.bed"),
        splicesites = temp(config['GenomesPrefix'] + "{GenomeName}/GTFTools/SpliceSites.bed"),
    log:
        "logs/gtftools/{GenomeName}.log"
    conda:
        "../envs/gtftools.yml"
    shell:
        """
        gtftools -e {output.exons} -i {output.introns} -q {output.splicesites} {input} &> {log}
        """

rule sort_and_index_gtftools_output:
    input:
        config['GenomesPrefix'] + "{GenomeName}/GTFTools/{RegionType}.bed",
    output:
        bed = config['GenomesPrefix'] + "{GenomeName}/GTFTools/{RegionType}.bed.gz",
        index = touch(config['GenomesPrefix'] + "{GenomeName}/GTFTools/{RegionType}.bed.gz.indexing_done")
    params:
        tabix_params = GetIndexingParams
    shell:
        """
        bedtools sort -i {input} | bgzip -c /dev/stdin > {output.bed}
        tabix {params.tabix_params} -f -p bed {output.bed} && touch {output.index}
        """

rule MakeIGV_GenomeFile:
    input:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
        bed = config['GenomesPrefix'] + "{GenomeName}/Reference.bed.gz",
    output:
        config['GenomesPrefix'] + "{GenomeName}/Reference.igv.genome.json"
    params:
        name = lambda wildcards: STAR_genomes.loc[wildcards.GenomeName]['Notes'],
        index_suffix = GetIndexSuffix
    run:
        import jinja2
        # Read the template file
        template_content = """{
          "id": "{{name}}",
          "name": "{{name}}",
          "fastaURL": "Reference.fa",
          "indexURL": "Reference.fa.fai",
          "tracks": [
            {
              "name": "Genes",
              "format": "bed",
              "url": "Reference.bed.gz",
              "indexURL": "Reference.bed.gz.{{index_suffix}}",
              "indexed": true,
              "hidden" : false,
              "searchable": true,
              "removable": false
            }
          ]
        }
        """
        # Create a Jinja2 template from the content
        template = jinja2.Template(template_content)
        # Render the template with the wildcards and params
        rendered_content = template.render(id=wildcards.GenomeName, name=params.name, index_suffix = params.index_suffix)
        # Write the rendered content to the output file
        with open(output[0], 'w') as f:
            f.write(rendered_content)
