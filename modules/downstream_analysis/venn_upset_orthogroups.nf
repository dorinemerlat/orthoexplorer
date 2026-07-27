process VENN_UPSET_ORTHOGROUPS {
    tag ""
    scratch false

    input:
    tuple path(orthogroups), path(gene_count), path(taxonomy), path(colors_yaml), val(clades)

    output:
    path("eligible_orthogroups.txt"), emit: eligible_orthogroups
    path("intersection_counts.tsv"), emit: intersection_counts
    path("intersection_gene_lists"), emit: intersection_gene_lists
    path("intersection_genes_all.tsv"), emit: intersection_genes
    path("orthogroup_presence_absence_by_clade.tsv"), emit: presence_absence
    path("upset_orthogroups.{png,pdf,svg}"), emit: upset
    path("venn_orthogroups.{png,pdf,svg}"), emit: venn

    script:
    if (!clades) {
        error "VENN_UPSET_ORTHOGROUPS requires at least two comma-separated clades"
    }

    """
    python_path="/shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python"

    PYTHONNOUSERSITE=1 \$python_path /shared/projects/metainvert/orthoexplorer/bin/venn_upset_orthogroups.py \
        --gene-counts "${gene_count}" \
        --orthogroups "${orthogroups}" \
        --taxonomy "${taxonomy}" \
        --clades "${clades}" \
        --outdir . \
        --colors "${colors_yaml}"
    """

    stub:
    """
    python_path="/shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python"
    \$python_path /shared/projects/metainvert/orthoexplorer/bin/venn_upset_orthogroups.py --help >/dev/null

    mkdir -p pangenome
    touch .stub
    """
}
