process VENN_UPSET_ORTHOGROUPS {
    tag '' 
    
    input:
    tuple path(orthogroups), path(gene_count), path(taxonomy), path(species_csv), path(outgroups_csv), path(colors_yaml), val(groups)

    output:
    path "pangenome", emit: results

    script:
    def groups_arg = groups ? "--groups '${groups}'" : ''

    """
    mkdir -p pangenome

    reformat_orthofinder_tables.py \
        --orthogroups ${orthogroups} \
        --gene-counts ${gene_count} \
        --species ${species_csv} \
        --outgroups ${outgroups_csv} \
        --output-prefix reformatted

    venn_upset_orthogroups.py \
        --gene-counts reformatted.Orthogroups.GeneCount.tsv \
        --orthogroups reformatted.Orthogroups.tsv \
        ${groups_arg} \
        --taxonomy ${taxonomy} \
        --outdir pangenome \
        --colors ${colors_yaml}
    """

    stub:
    """
    reformat_orthofinder_tables.py --help >/dev/null
    venn_upset_orthogroups.py --help >/dev/null

    mkdir -p pangenome
    touch pangenome/.stub
    """
}
