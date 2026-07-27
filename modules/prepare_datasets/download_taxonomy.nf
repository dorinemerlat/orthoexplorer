process DOWNLOAD_TAXONOMY {
    tag ""

    input:
    path ingroups_csv
    path outgroups_csv

    output:
    path "taxonomy.tsv"

    script:
    """
    # Write the ingroup datasets.
    awk -F',' '
        BEGIN {
            OFS=","
            print "id", "name", "taxid", "group"
        }

        FNR == 1 {
            next
        }

        {
            taxid = \$1
            name  = \$2

            id = tolower(name)

            # Replace each non-alphanumeric character with one hyphen.
            gsub(/[^a-z0-9]/, "-", id)
            gsub(/^-+|-+\$/, "", id)

            print id, name, taxid, "ingroup"
        }
    ' "${ingroups_csv}" > combined_species.csv

    # Append the outgroup datasets.
    awk -F',' '
        BEGIN {
            OFS=","
        }

        FNR == 1 {
            next
        }

        {
            taxid = \$1
            name  = \$2

            id = tolower(name)

            # Replace each non-alphanumeric character with one hyphen.
            gsub(/[^a-z0-9]/, "-", id)
            gsub(/^-+|-+\$/, "", id)

            print id, name, taxid, "outgroup"
        }
    ' "${outgroups_csv}" >> combined_species.csv

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