process DOWNLOAD_GENOME {
    tag "${name}"

    input:
    tuple val(name), val(meta)

    output:
    tuple val(name), val(meta), path("${name}_assembly.fna"), path("${name}_annotations.gff")

    script:
    """
    datasets download genome taxon ${meta.taxid} \
        --assembly-source RefSeq \
        --include genome,gff3 \
        --filename ncbi_dataset.zip

    unzip -q ncbi_dataset.zip

    genome_file=\$(find ncbi_dataset/data -type f -name '*_genomic.fna' | head -n 1)
    gff_file=\$(find ncbi_dataset/data -type f -name 'genomic.gff' | head -n 1)

    if [[ -z "\${genome_file}" || -z "\${gff_file}" ]]; then
        echo "ERROR: no RefSeq genome/GFF3 found for taxid ${meta.taxid}" >&2
        exit 1
    fi

    cp "\${genome_file}" ${name}_assembly.fna
    cp "\${gff_file}" ${name}_annotations.gff
    """

    stub:
    """
    datasets --help >/dev/null
    unzip -v >/dev/null

    touch ${name}_assembly.fna
    touch ${name}_annotations.gff
    """
}
