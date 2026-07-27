process ADD_ORTHOGROUPS_TO_ANNOTATIONS {
    tag "${name}"

    input:
    tuple val(name), val(meta), path(annotations)
    each path(orthogroups)

    output:
    tuple val(name), val(meta), path("${name}_go_with_orthogroups.tsv")

    script:
    """
    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python \
        /shared/projects/metainvert/orthoexplorer/bin/add_orthogroups_to_annotations.py \
        --annotations ${annotations} \
        --orthogroups ${orthogroups} \
        --output ${name}_go_with_orthogroups.tsv
    """

    stub:
    """
    touch ${name}_go_with_orthogroups.tsv
    """
}