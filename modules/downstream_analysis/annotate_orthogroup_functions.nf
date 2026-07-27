process ANNOTATE_ORTHOGROUP_FUNCTIONS {
    tag ""

    input:
    tuple path(orthogroups), path(gffs, stageAs: 'gff/*'), path(taxonomy)

    output:
    path "gene_function_annotations.tsv", emit: gene_annotations
    path "orthogroup_function_summary.tsv", emit: function_summary
    path "orthogroup_product_counts.tsv", emit: product_counts
    path "orthogroup_go_support.tsv", emit: go_support
    path "orthogroup2go_union.tsv", emit: orthogroup2go
    path "annotation_summary.tsv", emit: annotation_summary

    script:
    """
    annotate_orthogroup_functions.py \
        --orthogroups ${orthogroups} \
        --gffs ${gffs} \
        --species-mapping ${taxonomy} \
        --output-prefix orthogroup
    """
}