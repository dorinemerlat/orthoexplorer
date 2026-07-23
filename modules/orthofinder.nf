process ORTHOFINDER {
    tag ""
    scratch false
    cpus 50
    time '10d'
    memory '100GB'

    input:
    path(proteomes)

    output:
    tuple path("Results_Jun23/Orthogroups/Orthogroups.tsv"), path("Results_Jun23/Orthogroups/Orthogroups.GeneCount.tsv"), path("Results_Jun23/Orthogroups/Orthogroups_UnassignedGenes.tsv")

    script:
    """
    mkdir proteomes

    cp ${proteomes} proteomes/

    module load orthofinder
    
    orthofinder \
        -f proteomes \
        -t ${task.cpus} \
        -a ${task.cpus}

    mv proteomes/OrthoFinder/Results_Jun23 .
    """
}