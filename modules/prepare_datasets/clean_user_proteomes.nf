process CLEAN_USER_PROTEOMES {
    tag "${name}"
    stageInMode 'copy'

    input:
    tuple val(name), val(meta), path(fasta)

    output:
    tuple val(name), val(meta), path("${name}.fa")

    script:
    """
    # avoid modifying the original files in case they are needed later
    cp ${fasta} ${name}.fa.tmp

    # Remove terminal or internal '*' characters from protein sequences.
    # FASTA header lines are left unchanged.
    awk '
        /^>/ {
            print
            next
        }
        {
            gsub(/\\*/, "")
            print
        }
    ' ${name}.fa.tmp > ${name}.fa
    """
    
    stub:
    """
    touch ${name}.fa
    """
}