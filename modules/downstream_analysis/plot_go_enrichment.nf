process PLOT_GO_ENRICHMENT {
    tag "${study_set_name}/max_terms_${max_terms}/fdr_threshold_${fdr_threshold}"
    scratch false
    
    input:
    tuple val(study_set_name), path(enrichment_table)
    val max_terms
    val fdr_threshold

    output:
    tuple val(study_set_name), path("go_enrichment_dotplot.{png,pdf,svg}")

    script:
    """
    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python \
        /shared/projects/metainvert/orthoexplorer/bin/plot_go_enrichment.py \
        --input "${enrichment_table}" \
        --output-prefix go_enrichment_dotplot \
        --title "${study_set_name}" \
        --max-terms ${max_terms} \
        --fdr-threshold ${fdr_threshold}
    """

    stub:
    """
    plot_go_enrichment.py --help >/dev/null

    touch go_enrichment_dotplot.png
    touch go_enrichment_dotplot.pdf
    touch go_enrichment_dotplot.svg
    """
}