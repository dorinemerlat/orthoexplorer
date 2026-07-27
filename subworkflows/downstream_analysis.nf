include { ASSIGN_GENE_CONSERVATION_RANK } from "${projectDir}/modules/downstream_analysis/assign_gene_conservation_rank"
include { VENN_UPSET_ORTHOGROUPS        } from "${projectDir}/modules/downstream_analysis/venn_upset_orthogroups"
include { DOWNLOAD_GO_OBO               } from "${projectDir}/modules/downstream_analysis/download_go_obo"
include { FIND_ENRICHMENT               } from "${projectDir}/modules/downstream_analysis/find_enrichment"

workflow DOWNSTREAM_ANALYSIS {

    take:
    orthogroups
    gene_count
    unassigned_genes
    orthogroup2go
    tree
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

    DOWNLOAD_GO_OBO()

    VENN_UPSET_ORTHOGROUPS.out.intersection_gene_lists
        .flatten()
        .map { path -> [path.baseName, path] }
        .set { intersection_gene_lists }
    
    FIND_ENRICHMENT(intersection_gene_lists, VENN_UPSET_ORTHOGROUPS.out.background_orthogroups, orthogroup2go, DOWNLOAD_GO_OBO.out)
    
    // emit:
    // conservation = ASSIGN_GENE_CONSERVATION_RANK.out.results
    // pangenome     = VENN_UPSET_ORTHOGROUPS.out.results
}
