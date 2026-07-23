process KEEP_LONGEST_ISOFORM {
    tag "${specie}"
    // label 'agat'
    scratch false
    stageInMode 'copy'
    label 'retry_with_backoff'
    maxRetries 5
    
    input:
    tuple val(specie), val(meta), path(genome), path(gff)

    output:
    tuple val(specie), val(meta), path("${specie}.gff"), emit: gff
    tuple val(specie), val(meta), path("${specie}.fasta"), emit: fasta
    
    script:
    """
    # extract the longest isoform for each gene using AGAT
    module load agat
    agat_sp_keep_longest_isoform.pl --gff $gff -o ${specie}.gff

    # extract proteins sequences
    agat_sp_extract_sequences.pl -g ${specie}.gff -f $genome -o ${specie}.fasta --protein
    """

    stub:
    """
    touch ${specie}.gff ${specie}.fasta
    """
}