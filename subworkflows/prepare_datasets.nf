include { DOWNLOAD_GENOME           } from "${projectDir}/modules/prepare_datasets/download_genomes"
include { KEEP_LONGEST_ISOFORM      } from "${projectDir}/modules/prepare_datasets/keep_longest_isoform"
include { DOWNLOAD_GENE2GO          } from "${projectDir}/modules/prepare_datasets/download_gene2go"
include { ANNOTATE_NCBI_GFF_WITH_GO } from "${projectDir}/modules/prepare_datasets/annotate_NCBI_gff_with_go"
include { CLEAN_USER_FILES          } from "${projectDir}/modules/prepare_datasets/clean_user_files"
include { DOWNLOAD_TAXONOMY         } from "${projectDir}/modules/prepare_datasets/download_taxonomy"


workflow PREPARE_DATASETS {

    take:
    supplied_datasets
    datasets_to_download
    ingroups_csv
    outgroups_csv
    download_count

    main:
    if (download_count > 0) {       
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
        DOWNLOAD_GENE2GO()

        /*
        * Annotate the GFF files of downloaded genomes with GO terms from NCBI.
        */
        ANNOTATE_NCBI_GFF_WITH_GO(KEEP_LONGEST_ISOFORM.out.gff, DOWNLOAD_GENE2GO.out)
    
        /*
        * Both branches must now emit the same tuple structure:
        */
        KEEP_LONGEST_ISOFORM.out.fasta
            .join(ANNOTATE_NCBI_GFF_WITH_GO.out, by: [0,1])
    
    }
    

    // /*
    //  * OrthoFinder only needs the protein FASTA files.
    //  */
    // prepared_datasets
    //     .map { id, meta, protein_fasta, gff ->
    //         protein_fasta
    //     }
    //     .collect()
    //     .set { proteome_files }
    
    // /*
    //  * Format the fasta files and gff provided by the user
    //  */
    // CLEAN_USER_FILES(supplied_datasets)

    // DOWNLOAD_TAXONOMY(ingroups_csv, outgroups_csv)

    // emit:
    // datasets = prepared_datasets
    // proteomes = proteome_files
    // taxonomy = DOWNLOAD_TAXONOMY.out
}