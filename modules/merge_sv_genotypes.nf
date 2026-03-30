/*
========================================================================================
    MERGE SV GENOTYPES (WGS mode)
========================================================================================
    Merge per-sample genotyped SV BCFs into a multi-sample VCF.
    Filter to inversions (PASS, >= 1kb), code genotypes as 0/1/2 matrix.
*/

process MERGE_SV_GENOTYPES {
    tag "merge_sv_geno"
    label 'process_medium'

    publishDir "${params.outdir}/wgs_analyses/sv", mode: 'copy'

    input:
    path bcfs
    path csis

    output:
    path "sv_inversions.vcf.gz", emit: sv_vcf
    path "sv_genotype_matrix.tsv", emit: sv_matrix
    path "sv_stats.txt", emit: stats

    script:
    """
    set -euo pipefail

    # Merge all genotyped BCFs
    ls *.geno.bcf > bcf_list.txt
    bcftools merge --file-list bcf_list.txt -Oz -o all_sv.vcf.gz
    bcftools index -t all_sv.vcf.gz

    # Filter: inversions only, PASS, min size 1kb
    bcftools view -i 'INFO/SVTYPE="INV" && FILTER="PASS" && INFO/END-POS>=1000' \\
        all_sv.vcf.gz -Oz -o sv_inversions.vcf.gz
    bcftools index -t sv_inversions.vcf.gz

    N_INV=\$(bcftools view -H sv_inversions.vcf.gz | wc -l)
    N_SAMPLES=\$(bcftools query -l sv_inversions.vcf.gz | wc -l)

    # Code genotypes: 0/0→0, 0/1→1, 1/1→2, ./.→NA
    echo -ne "sv_id" > sv_genotype_matrix.tsv
    bcftools query -l sv_inversions.vcf.gz | while read s; do
        echo -ne "\\t\${s}" >> sv_genotype_matrix.tsv
    done
    echo "" >> sv_genotype_matrix.tsv

    bcftools query -f '%CHROM:%POS-%INFO/END[\\t%GT]\\n' sv_inversions.vcf.gz \\
        | sed 's|0/0|0|g; s|0/1|1|g; s|1/0|1|g; s|1/1|2|g; s|\\./\\.|NA|g' \\
        >> sv_genotype_matrix.tsv

    cat > sv_stats.txt <<EOF
SV Genotyping Statistics
========================
Total inversions (PASS, >= 1kb): \${N_INV}
Samples: \${N_SAMPLES}
EOF

    echo "SV genotyping complete: \${N_INV} inversions across \${N_SAMPLES} samples"
    """
}
