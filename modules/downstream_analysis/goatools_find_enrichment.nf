process GOATOOLS_FIND_ENRICHMENT {
    tag "${study_set_name}"

    input:
    tuple val(study_set_name), path(study_set)
    each path(population_set)
    each path(association)
    each path(obo)

    output:
    tuple val(study_set_name), path("${study_set_name}_go_enrichment.tsv")

    script:
    """
    # Extract unique orthogroups from the intersection gene table.
    tail -n +2 "${study_set}" |
        cut -f2 |
        sort -u \
        > study_orthogroups.txt

    goatools find_enrichment \
        --obo "${obo}" \
        --indent \
        --outfile "${study_set_name}_go_enrichment.tsv"
        study_orthogroups.txt \
        "${population_set}" \
        "${association}" \

    # GOATOOLS does not create an output file when no significant result is found.
    if [[ ! -s "${study_set_name}_go_enrichment.tsv" ]]; then
        touch "${study_set_name}_go_enrichment.tsv"
    fi
    """

    stub:
    """
    touch "${study_set_name}_go_enrichment.tsv"
    """
}