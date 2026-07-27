process BUILD_ORTHOGROUP2GO {
    tag "orthogroup2go"

    input:
    path annotation_tables

    output:
    path "orthogroup2go_union.tsv"

    script:
    """
    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python \
        /shared/projects/metainvert/orthoexplorer/bin/build_orthogroup2go.py \
        --annotations ${annotation_tables.join(' ')} \
        --output orthogroup2go_union.tsv
    """

    stub:
    """
    touch orthogroup2go_union.tsv
    """
}