include { PREPARE_DATASETS    } from "${projectDir}/subworkflows/prepare_datasets"
include { INFER_ORTHOLOGY     } from "${projectDir}/subworkflows/infer_orthology"
include { DOWNSTREAM_ANALYSIS } from "${projectDir}/subworkflows/downstream_analysis"


workflow ORTHOEXPLORER {

    take:
    supplied_datasets
    datasets_to_download
    ingroups_csv
    outgroups_csv
    tree
    colors
    cafe_filter_keywords

    main:
    PREPARE_DATASETS(
        supplied_datasets,
        datasets_to_download,
        ingroups_csv,
        outgroups_csv
    )

    INFER_ORTHOLOGY(
        PREPARE_DATASETS.out.proteomes,
        tree,
        PREPARE_DATASETS.out.taxonomy,
        PREPARE_DATASETS.out.annotations,
        cafe_filter_keywords
    )

    DOWNSTREAM_ANALYSIS(
        INFER_ORTHOLOGY.out.orthogroups,
        INFER_ORTHOLOGY.out.gene_count,
        INFER_ORTHOLOGY.out.unassigned_genes,
        INFER_ORTHOLOGY.out.tree,
        PREPARE_DATASETS.out.taxonomy,
        colors
    )

    // emit:
    // prepared_datasets = PREPARE_DATASETS.out.datasets
    // proteomes         = PREPARE_DATASETS.out.proteomes
    // taxonomy          = PREPARE_DATASETS.out.taxonomy
}