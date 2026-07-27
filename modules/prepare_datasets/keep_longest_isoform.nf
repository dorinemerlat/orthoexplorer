process KEEP_LONGEST_ISOFORM {
    tag "${name}"

    input:
    tuple val(name), val(meta), path(genome), path(gff)

    output:
    tuple val(name), val(meta), path("${name}.fa"), emit: fasta
    tuple val(name), val(meta), path("${name}.gff"), emit: gff

    script:
    """
    agat_sp_keep_longest_isoform.pl \
        --gff ${gff} \
        -o ${name}.gff

    agat_sp_extract_sequences.pl \
        -g ${name}.gff \
        -f ${genome} \
        --protein \
        -o ${name}.fa
    """

    stub:
    """
    agat_sp_keep_longest_isoform.pl --help >/dev/null
    agat_sp_extract_sequences.pl --help >/dev/null

    touch ${name}.gff
    touch ${name}.fa
    """
}
