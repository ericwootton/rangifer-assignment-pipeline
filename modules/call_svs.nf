/*
========================================================================================
    CALL SVs (WGS mode)
========================================================================================
    Per-sample structural variant calling using Delly.
    Three steps: per-sample discovery → merge sites → per-sample genotyping at merged sites.
*/

process CALL_SVS_DISCOVER {
    tag "${sample}"
    label 'process_medium'

    input:
    tuple val(sample), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(sample), path("${sample}.delly.bcf"), path("${sample}.delly.bcf.csi"), emit: bcf

    script:
    """
    delly call -g ${fasta} ${bam} -o ${sample}.delly.bcf
    """
}

process MERGE_SV_SITES {
    tag "merge_sv_sites"
    label 'process_medium'

    input:
    path bcfs
    path csis

    output:
    path "merged_sites.bcf", emit: sites_bcf
    path "merged_sites.bcf.csi", emit: sites_csi

    script:
    """
    ls *.delly.bcf > bcf_list.txt
    delly merge -o merged_sites.bcf \$(cat bcf_list.txt | tr '\\n' ' ')
    """
}

process GENOTYPE_SVS {
    tag "${sample}"
    label 'process_medium'

    input:
    tuple val(sample), path(bam), path(bai)
    path fasta
    path fai
    path sites_bcf
    path sites_csi

    output:
    tuple val(sample), path("${sample}.geno.bcf"), path("${sample}.geno.bcf.csi"), emit: bcf

    script:
    """
    delly call -g ${fasta} -v ${sites_bcf} ${bam} -o ${sample}.geno.bcf
    """
}
