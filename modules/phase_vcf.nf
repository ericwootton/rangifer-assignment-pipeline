/*
========================================================================================
    PHASE VCF (WGS mode)
========================================================================================
    Reference-free phasing using Eagle.
*/

process PHASE_VCF {
    tag "phase"
    label 'process_high'

    publishDir "${params.outdir}/wgs_analyses/phasing", mode: 'copy'

    input:
    path vcf
    path tbi

    output:
    tuple path("phased.vcf.gz"), path("phased.vcf.gz.tbi"), emit: phased_vcf

    script:
    """
    eagle \\
        --vcf ${vcf} \\
        --outPrefix phased \\
        --numThreads ${task.cpus}

    bcftools index -t phased.vcf.gz
    """
}
