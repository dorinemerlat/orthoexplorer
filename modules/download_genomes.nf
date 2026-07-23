process DOWNLOAD_GENOME {
    tag "$specie"

    input:
    tuple val(specie), val(meta)

    output:
    tuple val(specie), val(meta), path("${specie}.fna"), path("${specie}.gff3")

    script:
    """
    module load ncbi-datasets-cli
    datasets download genome taxon ${meta.taxid}  \
        --assembly-source RefSeq \
        --include genome,gff3 \

    unzip -q ncbi_dataset.zip

    cp ncbi_dataset/data/*/*_genomic.fna ${specie}.fna
    cp ncbi_dataset/data/*/genomic.gff ${specie}.gff3
    """
}