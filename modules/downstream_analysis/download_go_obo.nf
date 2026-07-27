process DOWNLOAD_GO_OBO {
    tag ""

    output:
    path("go-basic.obo")

    script:
    """
    wget https://current.geneontology.org/ontology/go-basic.obo
    """

    stub:
    """
    touch go-basic.obo
    """
}