process CLEAN_USER_FILES {
    tag "${name}"
    stageInMode 'copy'

    input:
    tuple val(name), val(meta), path(fasta), path(gff)

    output:
    tuple val(name), val(meta), path("${name}.fa"), path("${name}.gff")

    script:
    """
    # avoid modifying the original files in case they are needed later
    cp ${fasta} ${name}.fa.tmp
    cp ${gff} ${name}.gff.tmp

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


    # Remove any GFF3 features that do not have a corresponding sequence in the FASTA file.
    # agat_sp_keep_longest_isoform.pl \
    #     --gff ${name}.gff.tmp \
    #     -o ${name}.gff \

    cp ${name}.gff.tmp ${name}.gff
    """
    
    stub:
    """
    touch ${name}.fa
    touch ${name}.gff
    """
}