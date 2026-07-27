include { DOWNLOAD_GENOME           } from "${projectDir}/modules/prepare_datasets/download_genomes"
include { KEEP_LONGEST_ISOFORM      } from "${projectDir}/modules/prepare_datasets/keep_longest_isoform"
include { DOWNLOAD_GENE2GO          } from "${projectDir}/modules/prepare_datasets/download_gene2go"
include { ANNOTATE_NCBI_GFF_WITH_GO } from "${projectDir}/modules/prepare_datasets/annotate_NCBI_gff_with_go"
include { EXTRACT_USER_GFF_GO       } from "${projectDir}/modules/prepare_datasets/extract_user_gff_go"
include { CLEAN_USER_PROTEOMES      } from "${projectDir}/modules/prepare_datasets/clean_user_proteomes"
include { DOWNLOAD_TAXONOMY         } from "${projectDir}/modules/prepare_datasets/download_taxonomy"


workflow PREPARE_DATASETS {

    take:
    supplied_datasets
    datasets_to_download
    ingroups_csv
    outgroups_csv

    main:    
    /*
    * Download the genome and its annotation for datasets for which
    * no protein FASTA and GFF were supplied.
    */
    DOWNLOAD_GENOME(datasets_to_download)

    /*
    * Produce one canonical protein sequence per gene from each
    * downloaded genome annotation.
    */
    KEEP_LONGEST_ISOFORM(DOWNLOAD_GENOME.out)

    /*
    * Download gene2go.gz from NCBI and annotate the GFF files of downloaded genomes with GO terms.
    */
    trigger = datasets_to_download.first().map { true }
    DOWNLOAD_GENE2GO(trigger)

    /*
    * Annotate the GFF files of downloaded genomes with GO terms from NCBI.
    */
    ANNOTATE_NCBI_GFF_WITH_GO(KEEP_LONGEST_ISOFORM.out.gff, DOWNLOAD_GENE2GO.out)

    /*
    * Join the canonical protein sequences with the GO-annotated GFF files for each downloaded genome.
    */
    KEEP_LONGEST_ISOFORM.out.fasta
        .join(ANNOTATE_NCBI_GFF_WITH_GO.out, by: [0,1])
        .set { downloaded_datasets }

    /*
    * EXTRACT_USER_GFF_GO
    */
    user_gff_files = supplied_datasets.map { id, meta, protein_fasta, gff ->[ id, meta, gff ] }

    EXTRACT_USER_GFF_GO(user_gff_files)

    // /*
    // * Format the fasta files and gff provided by the user
    // */
    user_gff_files = supplied_datasets.map { id, meta, protein_fasta, gff ->[ id, meta, protein_fasta ] }
    CLEAN_USER_PROTEOMES(user_gff_files)

    /*
    * Join the canonical protein sequences with the GO-annotated GFF files for each supplied genome.
    */
    CLEAN_USER_PROTEOMES.out
        .join(EXTRACT_USER_GFF_GO.out, by: [0,1])
        .set { supplied_datasets }

    /*
    * Join the canonical protein sequences with the GO-annotated GFF files for each supplied genome.
    */
    downloaded_datasets
        .concat(supplied_datasets)
        .set { all_datasets }

    /*
     * OrthoFinder only needs the protein FASTA files.
     */
    all_datasets
        .map { id, meta, protein_fasta, go_tsv ->
            protein_fasta
        }
        .collect()
        .set { all_proteomes }

    DOWNLOAD_TAXONOMY(ingroups_csv, outgroups_csv)

    emit:
    datasets = all_datasets
    proteomes = all_proteomes
    taxonomy = DOWNLOAD_TAXONOMY.out
}