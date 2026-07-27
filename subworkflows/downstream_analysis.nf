include { ASSIGN_GENE_CONSERVATION_RANK } from "${projectDir}/modules/downstream_analysis/assign_gene_conservation_rank"
include { VENN_UPSET_ORTHOGROUPS        } from "${projectDir}/modules/downstream_analysis/venn_upset_orthogroups"

workflow DOWNSTREAM_ANALYSIS {

    take:
    orthogroups
    gene_count
    unassigned_genes
    tree
    datasets
    taxonomy
    colors

    main:
    conservation_analysis_inputs = orthogroups
        .combine(gene_count)
        .combine(unassigned_genes)
        .combine(tree)
        .combine(taxonomy)
        .combine(colors)
        .map { orthogroups, gene_count, unassigned_genes, tree, taxonomy, colors ->
            ['with_all_clades', orthogroups, gene_count, unassigned_genes, tree, taxonomy, colors, null ]
        }

    if (params.clades) {
        conservation_analysis_inputs_with_clades = conservation_analysis_inputs
            .map { status, orthogroups, gene_count, unassigned_genes, tree, taxonomy, colors, clades ->
                ['with_user_clades', orthogroups, gene_count, unassigned_genes, tree, taxonomy, colors, params.clades ]
            }
        
        conservation_analysis_inputs = conservation_analysis_inputs.concat(conservation_analysis_inputs_with_clades)

    } 

    ASSIGN_GENE_CONSERVATION_RANK(conservation_analysis_inputs)

    pangenome_analysis_inputs = conservation_analysis_inputs_with_clades
        .map { status, orthogroups, gene_count, unassigned_genes, tree, taxonomy, colors, clades ->
            [orthogroups, gene_count, taxonomy, colors, clades]
        }

    VENN_UPSET_ORTHOGROUPS(pangenome_analysis_inputs)

    
    // emit:
    // conservation = ASSIGN_GENE_CONSERVATION_RANK.out.results
    // pangenome     = VENN_UPSET_ORTHOGROUPS.out.results
}
