/*
========================================================================================
    EXTRACT VCF GENOTYPES
========================================================================================
    Uses bcftools to extract genotypes from the WGS-derived VCF into a TSV format
    that can be merged with chip data in the R-based merge step.
*/

process EXTRACT_VCF_GENOTYPES {
    tag "extract_vcf"
    label 'process_medium'

    publishDir "${params.outdir}/data_integration", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)

    output:
    path "wgs_genotypes.tsv", emit: genotype_tsv
    path "wgs_samples.txt",   emit: sample_list
    path "wgs_snp_info.tsv",  emit: snp_info

    script:
    """
    # Extract sample list
    bcftools query -l ${vcf} > wgs_samples.txt

    # Extract SNP info (CHROM POS REF ALT)
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${vcf} > wgs_snp_info.tsv

    # Extract genotypes as 0/0, 0/1, 1/1, ./.
    # Header line with sample names
    echo -ne "CHROM\\tPOS\\tREF\\tALT" > wgs_genotypes.tsv
    while read sample; do
        echo -ne "\\t\${sample}" >> wgs_genotypes.tsv
    done < wgs_samples.txt
    echo "" >> wgs_genotypes.tsv

    # Data lines
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ${vcf} >> wgs_genotypes.tsv

    echo "Extracted \$(wc -l < wgs_samples.txt) samples"
    echo "Extracted \$(tail -n+2 wgs_genotypes.tsv | wc -l) SNPs"
    """
}
