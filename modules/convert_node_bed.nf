/*
========================================================================================
    CONVERT NODE GENOTYPES TO PLINK BED
========================================================================================
    Converts a node's 012-coded reference genotype TSV to PLINK .bed/.bim/.fam
    for ADMIXTURE input at each hierarchy level.
*/

process CONVERT_NODE_BED {
    tag "${node_dir.name}"
    label 'process_low'

    input:
    tuple path(node_dir), path(node_meta)

    output:
    tuple val("${node_dir.name}"), path("${node_dir.name}.bed"), path("${node_dir.name}.bim"), path("${node_dir.name}.fam"), emit: plink_bed

    script:
    def node_name = node_dir.name
    """
    #!/usr/bin/env bash
    set -euo pipefail

    GENO="${node_dir}/reference_genotypes.tsv"

    # Convert 012 matrix to PLINK .ped format
    # Header has: sample_id SNP1 SNP2 ...
    # Data has: sampleID 0 1 2 NA ...

    # Build .map file from header (SNP names are CHROM:POS)
    # Use chromosome 1 for all SNPs (ADMIXTURE doesn't use chromosome structure)
    head -1 "\$GENO" | tr '\\t' '\\n' | tail -n +2 | awk -F: '{
        print 1, \$0, 0, NR
    }' > "${node_name}.map"

    # Build .ped file from genotype data
    tail -n +2 "\$GENO" | awk -F'\\t' '{
        # FID IID PID MID SEX PHENO genotypes...
        printf "%s %s 0 0 0 -9", \$1, \$1
        for (i = 2; i <= NF; i++) {
            if (\$i == "0") printf " A A"
            else if (\$i == "1") printf " A B"
            else if (\$i == "2") printf " B B"
            else printf " 0 0"
        }
        printf "\\n"
    }' > "${node_name}.ped"

    # Convert to binary PLINK format
    plink --file "${node_name}" --make-bed --out "${node_name}" --allow-no-sex --allow-extra-chr
    """
}
