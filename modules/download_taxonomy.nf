process DOWNLOAD_TAXONOMY {
    tag "${id}"

    input:
    tuple val(id), val(species), val(outgroups)

    output:
    tuple val(id), path("taxonomy_${id}.tsv")

    script:
    """
    cat ${species} | cut -f1 -d, > infos.txt
    tail -n +2 ${outgroups} | cut -f1 -d, >> infos.txt
    
    download_taxonomy.py -i infos.txt -o taxonomy_${id}.tsv
    """
}