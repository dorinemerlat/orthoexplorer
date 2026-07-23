process ASSIGN_GENE_CONSERVATION_RANK {
    tag "${id}"
    scratch false
    conda "/shared/projects/metainvert/orthoexplorer/envs/orthoexplorer"
    
    input:
    tuple val(id), path(orthogroups), path(gene_count), path(unassigned_genes), path(taxonomy), path(species), path(outgroups), path(colors_yaml)

    output:
    tuple val(id), path("outdir/")

    script:
    def orthogroups_reformatted = "Orthogroups.reformatted.tsv"
    def gene_count_reformatted = "Orthogroups.GeneCount.reformatted.tsv"
    def unassigned_genes_reformatted = "Orthogroups_UnassignedGenes.reformatted.tsv"

    """
    cp ${orthogroups} ${orthogroups_reformatted}
    cp ${gene_count} ${gene_count_reformatted}
    cp ${unassigned_genes} ${unassigned_genes_reformatted}

    mkdir -p outdir

    tail -n +2 ${species} | while IFS=',' read -r taxid species_name proteome; do
        proteome_base=\$(basename "\${proteome}")
        proteome_base=\${proteome_base%.*}

        sed -i "s/\${proteome_base}/\${species_name}/g" ${orthogroups_reformatted}
        sed -i "s/\${proteome_base}/\${species_name}/g" ${gene_count_reformatted}
        sed -i "s/\${proteome_base}/\${species_name}/g" ${unassigned_genes_reformatted}
    done

    tail -n +2 ${outgroups} | while IFS=',' read -r taxid species_name; do
        new_name=\$(echo "\${species_name}" | tr '[:upper:]' '[:lower:]' | sed 's/[ .]/-/g')

        sed -i "s/\${new_name}/\${species_name}/g" ${orthogroups_reformatted}
        sed -i "s/\${new_name}/\${species_name}/g" ${gene_count_reformatted}
        sed -i "s/\${new_name}/\${species_name}/g" ${unassigned_genes_reformatted}
    done

    tail -n +2 ${outgroups} | cut -f1 -d, > outgroups_taxid.txt

    assign_gene_conservation_rank.py \
        --gene-counts ${gene_count_reformatted} \
        --orthogroups ${orthogroups_reformatted} \
        --unassigned ${unassigned_genes_reformatted} \
        --taxonomy ${taxonomy} \
        --outdir outdir \
        --outgroups outgroups_taxid.txt \
        --colors ${colors_yaml} \
        --reorder-table
    """
}