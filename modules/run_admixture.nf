/*
========================================================================================
    ADMIXTURE
========================================================================================
    Run ADMIXTURE across K=2 to K=15 to estimate ancestry proportions.
    Uses cross-validation error to select optimal K.

    Expects pre-converted PLINK .bed/.bim/.fam from CONVERT_TO_BED process.
*/

process RUN_ADMIXTURE {
    tag "admixture"
    label 'process_high'

    publishDir "${params.outdir}/structure_discovery/admixture", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path "faststructure_results", emit: structure_results, type: 'dir'
    path "choosek_output.txt",    emit: choosek

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    N_SAMP=\$(wc -l < ${fam})
    echo "Running ADMIXTURE on \${N_SAMP} samples"

    mkdir -p faststructure_results

    MAX_K=${params.max_k_structure}
    if [ \$MAX_K -ge \$N_SAMP ]; then
        MAX_K=\$((N_SAMP - 1))
    fi

    BEST_K=2
    BEST_CV=999999

    for K in \$(seq 2 \$MAX_K); do
        echo "Running K = \$K"
        admixture --cv -j\${SLURM_CPUS_PER_TASK:-4} --seed=42 ${bed} \$K 2>&1 | tee faststructure_results/admixture_K\${K}.log

        # Move output files to results dir
        mv ${bed.baseName}.\${K}.Q faststructure_results/ 2>/dev/null || true
        mv ${bed.baseName}.\${K}.P faststructure_results/ 2>/dev/null || true

        # Parse CV error
        CV=\$(grep "CV error" faststructure_results/admixture_K\${K}.log | awk '{print \$NF}')
        if [ -n "\$CV" ]; then
            echo "K=\$K CV error: \$CV"
            LESS=\$(awk "BEGIN {print (\$CV < \$BEST_CV)}")
            if [ "\$LESS" = "1" ]; then
                BEST_CV=\$CV
                BEST_K=\$K
            fi
        fi
    done

    echo "Optimal K (lowest CV error) = \$BEST_K (CV = \$BEST_CV)" > choosek_output.txt
    echo "Model complexity that maximizes marginal likelihood = \$BEST_K" >> choosek_output.txt
    cat choosek_output.txt
    """
}
