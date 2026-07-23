process ASSIGN_GENE_CONSERVATION_RANK {
    tag '' 

    input:
    tuple path(orthogroups), path(gene_count), path(unassigned_genes), path(taxonomy), path(species_csv), path(outgroups_csv), path(colors_yaml)

    output:
    path "gene_conservation", emit: results

    script:
    """
    mkdir -p gene_conservation

    reformat_orthofinder_tables.py \
        --orthogroups ${orthogroups} \
        --gene-counts ${gene_count} \
        --unassigned ${unassigned_genes} \
        --species ${species_csv} \
        --outgroups ${outgroups_csv} \
        --output-prefix reformatted

    tail -n +2 ${outgroups_csv} \
        | cut -d',' -f1 \
        | sed '/^[[:space:]]*\$/d' \
        > outgroups_taxid.txt

    assign_gene_conservation_rank.py \
        --gene-counts reformatted.Orthogroups.GeneCount.tsv \
        --orthogroups reformatted.Orthogroups.tsv \
        --unassigned reformatted.Orthogroups_UnassignedGenes.tsv \
        --taxonomy ${taxonomy} \
        --outdir gene_conservation \
        --outgroups outgroups_taxid.txt \
        --colors ${colors_yaml} \
        --reorder-table
    """

    stub:
    """
    reformat_orthofinder_tables.py --help >/dev/null
    assign_gene_conservation_rank.py --help >/dev/null

    mkdir -p gene_conservation
    touch gene_conservation/.stub
    """
}
