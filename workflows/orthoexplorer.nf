include { DOWNLOAD_GENOME               } from '../modules/download_genomes'
include { KEEP_LONGEST_ISOFORM          } from '../modules/keep_longest_isoform'
include { ORTHOFINDER                   } from '../modules/orthofinder'
include { DOWNLOAD_TAXONOMY             } from '../modules/download_taxonomy'
include { ASSIGN_GENE_CONSERVATION_RANK } from '../modules/assign_gene_conservation_rank'
include { VENN_UPSET_ORTHOGROUPS        } from '../modules/venn_upset_orthogroups'

workflow ORTHOEXPLORER {

    main:

    Channel
        .fromPath(params.outgroups)
        .splitCsv(header: true, sep: ',')
        .map { row ->

            def taxid = row.taxid
            def specie = row.name

            def meta = [
                original_name: specie,
                taxid: taxid
            ]

            def specie_reformatted = specie
                .replaceAll(/\s+/, '-')
                .replaceAll(/\./, '-')
                .toLowerCase()

            tuple(specie_reformatted, meta)
        }
        .set { outgroup_species }

    DOWNLOAD_GENOME(outgroup_species)

    KEEP_LONGEST_ISOFORM(DOWNLOAD_GENOME.out)

    Channel
        .fromPath(params.species)
        .splitCsv(header: true, sep: ',')
        .map { row ->

            def taxid = row.taxid
            def specie = row.name
            def proteome = file(row.file)

            def meta = [
                original_name: specie,
                taxid: taxid
            ]

            def specie_reformatted = specie
                .replaceAll(/\s+/, '-')
                .replaceAll(/\./, '-')
                .toLowerCase()

            tuple(specie_reformatted, meta, proteome)
        }
        .set { ingroup_species }

    ingroup_species
        .concat(KEEP_LONGEST_ISOFORM.out.fasta)
        .map {name, meta, file -> ['all', file] }
        .groupTuple() 
        .map { name, files -> files } 
        .set { files } 

    ORTHOFINDER(files)

    // Download taxonomy
    Channel
        .from([[ file(params.species), file(params.outgroups) ]])
        .map { species, outgroups -> ['all', species, outgroups] } 
        .set { all_species }

    Channel
        .from([[ file("/shared/projects/metainvert/orthoexplorer/species_best.csv"), file(params.outgroups) ]])
        .map { species, outgroups -> ['best', species, outgroups] } 
        .set { best_species }

    all_species
        .concat(best_species)
        .set { all_inputs }

    DOWNLOAD_TAXONOMY(all_inputs)

    DOWNLOAD_TAXONOMY.out
        .branch {
            all:  it[0] == 'all'
            best: it[0] == 'best'
        }
        .set { taxonomy }

    taxonomy.all
        .combine(ORTHOFINDER.out)
        .map { id, taxonomy, orthogroups, gene_count, unassigned_genes -> [id, orthogroups, gene_count, unassigned_genes, taxonomy, params.species, params.outgroups, params.colors ] } 
        .set { all_inputs }
    
    all_inputs
        .map { id, orthogroups, gene_count, unassigned_genes, taxonomy, species, outgroups, colors ->
            tuple(
                'best',
                file(orthogroups.toString().replace('Jun23', 'Jun23_1')),
                file(gene_count.toString().replace('Jun23', 'Jun23_1')),
                file(unassigned_genes.toString().replace('Jun23', 'Jun23_1')),
                file('/shared/projects/metainvert/orthoexplorer/cache/download_taxonomy/best/taxonomy_best.tsv'),
                file(params.species),
                file(params.outgroups),
                file(params.colors)
            )
        }
        .set { best_inputs }

    ASSIGN_GENE_CONSERVATION_RANK(all_inputs.concat(best_inputs)) 

    VENN_UPSET_ORTHOGROUPS(all_inputs.concat(best_inputs))
}