/*
========================================================================================
    COMPUTE FST + DXY (WGS mode, per node)
========================================================================================
    Uses pixy to compute windowed Fst and Dxy between groups at each hierarchy node.
*/

process COMPUTE_FST_DXY {
    tag "${node_name}"
    label 'process_medium'

    publishDir "${params.outdir}/wgs_analyses/fst_dxy/${node_name}", mode: 'copy'

    input:
    tuple val(node_name), path(node_dir), path(node_meta)
    tuple path(vcf), path(tbi)

    output:
    tuple val(node_name), path("${node_name}_pixy_fst.txt"), path("${node_name}_pixy_dxy.txt"), emit: results

    script:
    """
    set -euo pipefail

    # Create populations file from node labels
    awk -F'\\t' 'NR>1 {print \$1"\\t"\$2}' ${node_dir}/reference_labels.tsv > populations.txt

    # Extract samples for this node
    cut -f1 populations.txt > node_samples.txt

    # Subset VCF to node samples
    bcftools view -S node_samples.txt ${vcf} -Oz -o node.vcf.gz
    bcftools index -t node.vcf.gz

    # Run pixy
    pixy \\
        --stats fst dxy \\
        --vcf node.vcf.gz \\
        --populations populations.txt \\
        --window_size 10000 \\
        --n_cores ${task.cpus} \\
        --output_prefix ${node_name}_pixy

    # Rename outputs
    if [ -f ${node_name}_pixy_fst.txt ]; then
        echo "Fst computed for node ${node_name}"
    fi
    if [ -f ${node_name}_pixy_dxy.txt ]; then
        echo "Dxy computed for node ${node_name}"
    fi
    """
}
