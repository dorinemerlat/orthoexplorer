process DOWNLOAD_GENE2GO {
    tag ""

    output:
    path("gene2go.gz")
    path("gene2go.gz")

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    """

    stub:
    """
    touch gene2go.gz
    """
}
