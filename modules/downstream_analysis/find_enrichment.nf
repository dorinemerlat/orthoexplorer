process FIND_ENRICHMENT {
    tag "${study_set_name}"

    input:
    tuple val(study_set_name), path(study_set)
    each path(population_set)
    each path(orthogroup2go)
    each path(obo)

    output:
    tuple val(study_set_name), path("${study_set_name}_go_enrichment.tsv"), emit: enrichment

    script:
    """
    # Extract unique orthogroups from the intersection gene table.
    tail -n +2 ${study_set} |
        cut -f 2 |
        sort -u \
        > study_orthogroups.txt

    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/find_enrichment.py\
        study_orthogroups.txt \
        ${population_set} \
        ${orthogroup2go} \
        --obo ${obo} \
        --alpha 0.05 \
        --method fdr_bh \
        --pval_field fdr_bh \
        --outfile ${study_set_name}_go_enrichment.tsv
    """

    stub:
    """
    touch ${study_set_name}_go_enrichment.tsv
    """
}