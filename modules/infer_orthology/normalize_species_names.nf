process NORMALIZE_SPECIES_NAMES {
    tag ""
    stageInMode 'copy'

    input:
    path(file)
    path(taxonomy)

    output:
    path "${file.name}"

    script:
    """
    tail -n +2 "${taxonomy}" |
    while IFS=\$'\\t' read -r id name rest
    do
        [ -z "\${id}" ] && continue
        [ -z "\${name}" ] && continue

        if [[ "${file.name}" == *.txt ]]; then
            replacement="\\\"\${name}\\\""
        else
            replacement="\${name}"
        fi

        sed -i "1s|\${id}|\${replacement}|g" "${file.name}"
    done
    """

    stub:
    """
    touch "${file.name}"
    """
}