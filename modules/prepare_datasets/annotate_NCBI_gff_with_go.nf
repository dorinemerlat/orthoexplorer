process ANNOTATE_NCBI_GFF_WITH_GO {
    tag "${name}"
    stageInMode 'copy'

    input:
    tuple val(name), val(meta), path(gff)
    path(gene2go)

    output:
    tuple val(name), val(meta), path("${name}_go.tsv")

    script:
    """
    /shared/projects/metainvert/orthoexplorer/envs/orthoexplorer/bin/python \
        /shared/projects/metainvert/orthoexplorer/bin/annotate_NCBI_gff_with_go.py \
        --gff ${gff} \
        --gene2go ${gene2go} \
        --output ${name}_go.tsv \
        --taxid $meta.taxid


    """
    
    stub:
    """
    touch ${name}_go.tsv
    """
}