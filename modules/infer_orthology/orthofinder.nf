process ORTHOFINDER {
    tag ''

    input:
    path proteomes

    output:
    tuple path("Orthogroups.tsv"), path("Orthogroups.GeneCount.tsv"), path("Orthogroups_UnassignedGenes.tsv"), emit: results
    path "SpeciesTree_rooted.txt", emit: species_tree
    path "orthofinder_results", emit: full_results

    script:
    """
    mkdir -p proteomes

    for proteome in ${proteomes}; do
        cp "\${proteome}" proteomes/
    done

    orthofinder \
        -f proteomes \
        -t ${task.cpus} \
        -a ${task.cpus}

    results_dir=\$(find proteomes/OrthoFinder -maxdepth 1 -type d -name 'Results_*' | sort | tail -n 1)

    if [[ -z "\${results_dir}" ]]; then
        echo "ERROR: OrthoFinder results directory was not found" >&2
        exit 1
    fi

    cp "\${results_dir}/Orthogroups/Orthogroups.tsv" Orthogroups.tsv
    cp "\${results_dir}/Orthogroups/Orthogroups.GeneCount.tsv" Orthogroups.GeneCount.tsv
    cp "\${results_dir}/Orthogroups/Orthogroups_UnassignedGenes.tsv" Orthogroups_UnassignedGenes.tsv
    cp "\${results_dir}/Species_Tree/SpeciesTree_rooted.txt" SpeciesTree_rooted.txt
    mv "\${results_dir}" orthofinder_results
    """

    stub:
    """
    orthofinder -h >/dev/null

    mkdir -p orthofinder_results/Orthogroups

    touch Orthogroups.tsv
    touch Orthogroups.GeneCount.tsv
    touch Orthogroups_UnassignedGenes.tsv
    touch SpeciesTree_rooted.txt

    cp Orthogroups.tsv orthofinder_results/Orthogroups/
    cp Orthogroups.GeneCount.tsv orthofinder_results/Orthogroups/
    cp Orthogroups_UnassignedGenes.tsv orthofinder_results/Orthogroups/
    cp SpeciesTree_rooted.txt orthofinder_results/Orthogroups/
    """
}
