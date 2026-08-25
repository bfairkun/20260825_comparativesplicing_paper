
rule Download_100WayVertebrate_MAF_hg38:
    output:
        "multiple_sequence_alignment_analysis/maf/{fn_base}.maf.gz"
    params:
        expected_md5=lambda wildcards: maf_filelist.loc[maf_filelist['fn_base']==wildcards.fn_base, "md5sum"].iloc[0],
        link = lambda wildcards: maf_filelist.loc[maf_filelist['fn_base']==wildcards.fn_base, 'link'].iloc[0], 
    log:
        "logs/Download_100WayVertebrate_MAF_hg38/{fn_base}.log"
    shell:
        """
        wget -O {output} {params.link} 2> {log}
        actual_md5=$(md5sum {output} | awk '{{print $1}}')
        if [ "$actual_md5" != "{params.expected_md5}" ]; then
            echo "ERROR: md5sum mismatch for {wildcards.fn_base}" >&2
            rm -f {output}
            exit 1
        fi
        """

rule unzip_MAF_and_index:
    input:
        "multiple_sequence_alignment_analysis/maf/{fn_base}.maf.gz"
    output:
        maf = "multiple_sequence_alignment_analysis/maf/{fn_base}.maf",
        index = "multiple_sequence_alignment_analysis/maf/{fn_base}.mafindex"
    conda:
        "../envs/bedparse.yml"
    resources:
        mem_mb = 48000
    shell:
        """
        zcat {input} > {output.maf}
        python scripts/IndexMAF.py {output.maf} {output.index} hg38.{wildcards.fn_base}
        """

rule Collect_MAF:
    input:
        expand("multiple_sequence_alignment_analysis/maf/{fn_base}.mafindex", fn_base=maf_filelist['fn_base'])
