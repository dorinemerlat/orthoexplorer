include { DOWNLOAD_GENOME      } from "${projectDir}/modules/prepare_datasets/download_genomes"
include { KEEP_LONGEST_ISOFORM } from "${projectDir}/modules/prepare_datasets/keep_longest_isoform"
include { CLEAN_USER_FILES     } from "${projectDir}/modules/prepare_datasets/clean_user_files"
include { DOWNLOAD_TAXONOMY    } from "${projectDir}/modules/prepare_datasets/download_taxonomy"


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
     *
     * Expected output:
     * [ id, meta, protein_fasta, gff ]
     */
    KEEP_LONGEST_ISOFORM(
        DOWNLOAD_GENOME.out
    )

    /*
     * Both branches must now emit the same tuple structure:
     *
     * [ id, meta, protein_fasta, gff ]
     */
    prepared_datasets = supplied_datasets
        .concat(KEEP_LONGEST_ISOFORM.out)

    /*
     * OrthoFinder only needs the protein FASTA files.
     */
    prepared_datasets
        .map { id, meta, protein_fasta, gff ->
            protein_fasta
        }
        .collect()
        .set { proteome_files }
    
    /*
     * Format the fasta files and gff provided by the user
     */
    CLEAN_USER_FILES(supplied_datasets)

    DOWNLOAD_TAXONOMY(ingroups_csv, outgroups_csv)

    emit:
    datasets = prepared_datasets
    proteomes = proteome_files
    taxonomy = DOWNLOAD_TAXONOMY.out
}