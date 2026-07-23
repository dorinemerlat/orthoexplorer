include { ASSIGN_GENE_CONSERVATION_RANK } from "${projectDir}/modules/downstream_analysis/assign_gene_conservation_rank"
include { VENN_UPSET_ORTHOGROUPS        } from "${projectDir}/modules/downstream_analysis/venn_upset_orthogroups"

workflow DOWNSTREAM_ANALYSIS {

    take:
    orthofinder_results
    taxonomy
    tree
    species_csv
    outgroups_csv
    colors_yaml

    main:
    analysis_inputs = orthofinder_results
        .combine(taxonomy)
        .map { orthogroups, gene_count, unassigned_genes, taxonomy_file ->
            [ orthogroups, gene_count, unassigned_genes, taxonomy_file, species_csv, outgroups_csv, colors_yaml ]
        }

    ASSIGN_GENE_CONSERVATION_RANK(analysis_inputs)

    pangenome_inputs = orthofinder_results
        .combine(taxonomy)
        .map {orthogroups, gene_count, unassigned_genes, taxonomy_file ->
            [ orthogroups, gene_count, taxonomy_file, species_csv, outgroups_csv, colors_yaml, params.clades ] 
        }

    VENN_UPSET_ORTHOGROUPS(pangenome_inputs)

    emit:
    conservation = ASSIGN_GENE_CONSERVATION_RANK.out.results
    pangenome     = VENN_UPSET_ORTHOGROUPS.out.results
}
