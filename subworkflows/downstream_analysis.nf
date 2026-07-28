include { ASSIGN_GENE_CONSERVATION_RANK     } from "${projectDir}/modules/downstream_analysis/assign_gene_conservation_rank"
include { VENN_UPSET_ORTHOGROUPS            } from "${projectDir}/modules/downstream_analysis/venn_upset_orthogroups"
include { DOWNLOAD_GO_OBO                   } from "${projectDir}/modules/downstream_analysis/download_go_obo"
include { DOWNLOAD_GOSLIM                   } from "${projectDir}/modules/downstream_analysis/download_goslim"
include { GOATOOLS_FIND_ENRICHMENT          } from "${projectDir}/modules/downstream_analysis/goatools_find_enrichment"
include { GOATOOLS_FIND_ENRICHMENT_GOSLIM   } from "${projectDir}/modules/downstream_analysis/goatools_find_enrichment_go_slim"
include { PLOT_GO_ENRICHMENT                } from "${projectDir}/modules/downstream_analysis/plot_go_enrichment"

workflow DOWNSTREAM_ANALYSIS {

    take:
    orthogroups
    gene_count
    unassigned_genes
    orthogroup2go
    population_orthogroups
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
    DOWNLOAD_GOSLIM()

    VENN_UPSET_ORTHOGROUPS.out.intersection_gene_lists
        .flatten()
        .map { path -> [path.baseName, path] }
        .set { intersection_gene_lists }
    
    GOATOOLS_FIND_ENRICHMENT(intersection_gene_lists, population_orthogroups, orthogroup2go, DOWNLOAD_GO_OBO.out)
    GOATOOLS_FIND_ENRICHMENT_GOSLIM(intersection_gene_lists, population_orthogroups, orthogroup2go, DOWNLOAD_GO_OBO.out, DOWNLOAD_GOSLIM.out)
    
    GOATOOLS_FIND_ENRICHMENT.out
        .filter { label, file -> file.size() > 0 }
        .map { label, file -> ['without_goslim/' + label, file] }
        .set { enrichment_results_without_goslim }

    GOATOOLS_FIND_ENRICHMENT_GOSLIM.out
        .filter { label, file -> file.size() > 0 }
        .map { label, file -> ['with_goslim/' + label, file] }
        .set { enrichment_results_with_goslim }

    enrichment_results = enrichment_results_without_goslim.concat(enrichment_results_with_goslim)

    PLOT_GO_ENRICHMENT(enrichment_results, 20, 0.05)
    
    // emit:
    // conservation = ASSIGN_GENE_CONSERVATION_RANK.out.results
    // pangenome     = VENN_UPSET_ORTHOGROUPS.out.results
}
