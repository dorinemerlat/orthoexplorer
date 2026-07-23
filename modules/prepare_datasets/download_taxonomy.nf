process DOWNLOAD_TAXONOMY {
    tag ""

    input:
    path ingroups_csv
    path outgroups_csv

    output:
    path "taxonomy.tsv"

    script:
    """
    # Combine both CSV files and generate a normalized dataset ID
    awk -F',' '
        BEGIN {
            OFS=","
            print "id", "name", "taxid"
        }

        FNR == 1 {
            next
        }

        {
            taxid = \$1
            name  = \$2

            id = tolower(name)
            gsub(/[^a-z0-9]+/, "-", id)
            gsub(/^-+|-+\$/, "", id)

            print id, name, taxid
        }
    ' ${ingroups_csv} ${outgroups_csv} > combined_species.csv

    download_taxonomy.py \
        --input combined_species.csv \
        --output taxonomy.tsv
    """

    stub:
    """
    download_taxonomy.py --help >/dev/null
    touch taxonomy.tsv
    """
}