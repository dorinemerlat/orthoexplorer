process VENN_UPSET_ORTHOGROUPS {
    tag ""
    scratch false

    input:
    tuple path(orthogroups), path(gene_count), path(taxonomy), path(colors_yaml), val(clades)

    output:
    path "pangenome", emit: results

    script:
    if (!clades) {
        error "VENN_UPSET_ORTHOGROUPS requires at least two comma-separated clades"
    }

    """
    python_path="/shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python"
    mkdir -p pangenome

    PYTHONNOUSERSITE=1 \$python_path /shared/projects/metainvert/orthoexplorer/bin/venn_upset_orthogroups.py \
        --gene-counts "${gene_count}" \
        --orthogroups "${orthogroups}" \
        --taxonomy "${taxonomy}" \
        --clades "${clades}" \
        --outdir pangenome \
        --colors "${colors_yaml}"
    """

    stub:
    """
    python_path="/shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python"
    \$python_path /shared/projects/metainvert/orthoexplorer/bin/venn_upset_orthogroups.py --help >/dev/null

    mkdir -p pangenome
    touch pangenome/.stub
    """
}
