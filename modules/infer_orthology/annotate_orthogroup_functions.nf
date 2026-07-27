process ANNOTATE_ORTHOGROUP_FUNCTIONS {
    tag "orthogroup-functional-annotations"
    memory "8 GB"
    
    input:
    path(annotation_tables)

    output:
    path("orthogroup_functional_annotations.tsv")

    script:
    """
    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python \
        /shared/projects/metainvert/orthoexplorer/bin/annotate_orthogroup_functions.py \
        --annotations ${annotation_tables.join(' ')} \
        --output orthogroup_functional_annotations.tsv
    """

    stub:
    """
    touch orthogroup_functional_annotations.tsv
    """
}