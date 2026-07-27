process ASSIGN_GENE_CONSERVATION_RANK {
    tag "$status"
    scratch false

    input:
    tuple val(status), path(orthogroups), path(gene_count), path(unassigned_genes), path(tree), path(taxonomy), path(colors_yaml), val(clades)

    output:
    path "$status"

    script:
    clades_arg = status == "with_user_clades" && clades ? "--clades '${clades}'" : ""

    """
    python_path="/shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python"
    mkdir -p "$status"

    \$python_path /shared/projects/metainvert/orthoexplorer/bin/assign_gene_conservation_rank.py \
        --gene-counts "${gene_count}" \
        --orthogroups "${orthogroups}" \
        --unassigned "${unassigned_genes}" \
        --tree "${tree}" \
        --taxonomy "${taxonomy}" \
        --outdir "$status" \
        --colors "${colors_yaml}" \
        ${clades_arg}
    """

    stub:
    """
    conda activate ${projectDir}/envs/orthoexplorer
    assign_gene_conservation_rank.py --help >/dev/null

    mkdir -p "$status"
    touch "$status/.stub"
    """
}