process ADD_CAFE_FILTER_STATUS {
    tag "${ratio_threshold}"

    input:
    path annotations
    path blacklist
    val ratio_threshold

    output:
    path "orthogroups_annotation_cafe.tsv"

    script:
    """
    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python \
        /shared/projects/metainvert/orthoexplorer/bin/add_cafe_filter_status.py \
        --annotations ${annotations} \
        --blacklist ${blacklist} \
        --ratio-threshold ${ratio_threshold} \
        --output orthogroups_annotation_cafe.tsv
    """

    stub:
    """
    touch orthogroups_annotation_cafe.tsv
    """
}