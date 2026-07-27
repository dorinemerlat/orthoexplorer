process FILTER_LARGE_GENE_FAMILIES {
    tag "${max_copy_number}"

    input:
    path table
    path gene_count
    val max_copy_number

    output:
    path "orthogroups_annotation_cafe_filtered.tsv", emit: filtered
    path "large_gene_families.tsv", emit: removed
    path "Orthogroups.GeneCount.family_size_filtered.tsv", emit: small_family_gene_count
    path "Orthogroups.GeneCount.large_gene_families.tsv", emit: large_family_gene_count

    script:
    """
    awk -F'\t' '
    NR==1 {
        print > "orthogroups_annotation_cafe_filtered.tsv"
        print > "large_gene_families.tsv"
        next
    }

    \$11 <= ${max_copy_number} {
        print > "orthogroups_annotation_cafe_filtered.tsv"
        next
    }

    {
        print > "large_gene_families.tsv"
    }
    ' $table

    cut -f 1 orthogroups_annotation_cafe_filtered.tsv | tail -n +2 | sort > filtered_orthogroups.txt

    head -n 1 $gene_count > Orthogroups.GeneCount.family_size_filtered.tsv
    grep -f filtered_orthogroups.txt $gene_count >> Orthogroups.GeneCount.family_size_filtered.tsv

    head -n 1 $gene_count > Orthogroups.GeneCount.large_gene_families.tsv
    grep -vf filtered_orthogroups.txt $gene_count >> Orthogroups.GeneCount.large_gene_families.tsv
    """
}