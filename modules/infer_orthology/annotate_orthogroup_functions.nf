process ANNOTATE_ORTHOGROUP_FUNCTIONS {
    tag ""
    memory "8 GB"

    input:
    path(annotation_tables)

    output:
    path("orthogroup_functional_annotations.tsv"), emit: orthogroup
    path("orthogroup_functional_annotations_without_all_functions.tsv"), emit: orthogroup_without_all_functions

    script:
    """
    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python \
        /shared/projects/metainvert/orthoexplorer/bin/annotate_orthogroup_functions.py \
        --annotations ${annotation_tables.join(' ')} \
        --output orthogroup_functional_annotations.tsv

    cut -f 1,2,3,4,5,6,7,8,9,10 orthogroup_functional_annotations.tsv > orthogroup_functional_annotations_without_all_functions.tsv
    """

    stub:
    """
    touch orthogroup_functional_annotations.tsv
    """
}