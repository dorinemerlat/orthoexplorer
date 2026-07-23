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
    colors_yaml

    main:
    PREPARE_DATASETS(
        supplied_datasets,
        datasets_to_download,
        ingroups_csv,
        outgroups_csv
    )

    INFER_ORTHOLOGY(
        PREPARE_DATASETS.out.proteomes,
        tree
    )

    // DOWNSTREAM_ANALYSIS(
    //     INFER_ORTHOLOGY.out.orthofinder_results,
    //     PREPARE_DATASETS.out.taxonomy,
    //     INFER_ORTHOLOGY.out.tree,
    //     ingroups_csv,
    //     outgroups_csv,
    //     colors_yaml,
    // )

    // emit:
    // prepared_datasets = PREPARE_DATASETS.out.datasets
    // proteomes         = PREPARE_DATASETS.out.proteomes
    // taxonomy          = PREPARE_DATASETS.out.taxonomy
}