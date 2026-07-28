process DOWNLOAD_GOSLIM {
    tag ""

    output:
    path("goslim_generic.obo")

    script:
    """
    wget https://current.geneontology.org/ontology/subsets/goslim_generic.obo
    """

    stub:
    """
    touch goslim_generic.obo
    """
}