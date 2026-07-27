process EXTRACT_USER_GFF_GO {
    tag "${id}"

    input:
    tuple val(id), val(meta), path(gff)

    output:
    tuple val(id), val(meta), path("${id}_go.tsv")

    script:
    """
    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python \
        /shared/projects/metainvert/orthoexplorer/bin/extract_user_gff_go.py \
        --gff ${gff} \
        --species-name "${meta.name}" \
        --output ${id}_go.tsv
    """

    stub:
    """
    touch ${id}_go.tsv
    """
}