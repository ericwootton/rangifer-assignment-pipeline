/*
========================================================================================
    CONVERT GENOTYPE MATRIX TO PLINK BED FORMAT
========================================================================================
    Converts the 012-coded genotype matrix to PLINK .bed/.bim/.fam format
    for use by fastSTRUCTURE.
*/

process CONVERT_TO_BED {
    tag "convert_bed"
    label 'process_low'

    input:
    path genotype_matrix

    output:
    tuple path("fs_input.bed"), path("fs_input.bim"), path("fs_input.fam"), emit: plink_bed

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    echo "Preparing PLINK input from genotype matrix..."

    # Create .map file (dummy chromosome since we just need structure analysis)
    head -1 ${genotype_matrix} | tr '\\t' '\\n' | tail -n+2 | awk '{print 1, \$0, 0, NR}' OFS='\\t' > input.map

    # Create .ped file from genotype matrix
    # PED format: FID IID PID MID Sex Pheno Geno1a Geno1b ...
    # Convert 012 to allele pairs: 0 -> A A, 1 -> A B, 2 -> B B, NA -> 0 0
    tail -n+2 ${genotype_matrix} | awk -F'\\t' '{
        sample = \$1
        printf "%s\\t%s\\t0\\t0\\t0\\t0", sample, sample
        for (i = 2; i <= NF; i++) {
            if (\$i == "" || \$i == "NA") printf "\\t0\\t0"
            else if (\$i == 0) printf "\\tA\\tA"
            else if (\$i == 1) printf "\\tA\\tB"
            else if (\$i == 2) printf "\\tB\\tB"
            else printf "\\t0\\t0"
        }
        printf "\\n"
    }' > input.ped

    N_SAMPLES=\$(wc -l < input.ped)
    N_SNPS=\$(wc -l < input.map)
    echo "Converting \${N_SAMPLES} samples x \${N_SNPS} SNPs to bed format..."

    plink --file input --make-bed --out fs_input --allow-no-sex --allow-extra-chr 2>&1
    """
}
