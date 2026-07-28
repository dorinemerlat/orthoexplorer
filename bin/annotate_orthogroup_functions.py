#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


UNINFORMATIVE_PRODUCTS = {
    "",
    "hypothetical protein",
    "uncharacterized protein",
    "unknown protein",
}

OUTPUT_COLUMNS = [
    "orthogroup",
    "representative_product",
    "representative_product_count",
    "representative_product_fraction",
    "annotated_gene_count",
    "total_gene_count",
    "product_diversity",
    "species_count",
    "species_fraction",
    "min_copy_number",
    "max_copy_number",
    "max_copy_number_species",
    "copy_number_ratio",
    "all_products",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Generate representative functional annotations for OrthoFinder orthogroups.")
    parser.add_argument("--annotations", required=True, nargs="+", type=Path, help="GO annotation TSV files containing orthogroup assignments.")
    parser.add_argument("--output", required=True, type=Path, help="Output orthogroup annotation TSV file.")
    return parser.parse_args()


def clean_product(product):
    """Apply light normalization while preserving readable product names."""
    if pd.isna(product):
        return ""

    product = str(product).strip()
    product = re.sub(r"\s+", " ", product)
    product = product.rstrip(".;")

    return product


def product_key(product):
    """Return a normalized key used to compare product names."""
    return clean_product(product).casefold()


def read_annotations(annotation_paths):
    """Read and combine all per-species annotation tables."""
    tables = []

    required_columns = {
        "specie_name",
        "orthogroup",
        "id",
        "product",
    }

    for annotation_path in annotation_paths:
        table = pd.read_csv(
            annotation_path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )

        missing_columns = required_columns - set(table.columns)

        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"Missing required columns in {annotation_path}: {missing}")

        tables.append(table[list(required_columns)])

    if not tables:
        return pd.DataFrame(columns=sorted(required_columns))

    annotations = pd.concat(tables, ignore_index=True)

    # Genes without an OrthoFinder orthogroup cannot contribute to this table.
    annotations = annotations[
        annotations["orthogroup"].str.strip() != ""
    ].copy()

    annotations["product"] = annotations["product"].map(clean_product)
    annotations["product_key"] = annotations["product"].map(product_key)

    return annotations


def select_gene_products(annotations):
    """Keep one product annotation per gene despite repeated GO rows."""
    gene_columns = [
        "orthogroup",
        "specie_name",
        "id",
    ]

    # Check whether a gene has conflicting non-empty product annotations.
    product_counts = (
        annotations[annotations["product_key"] != ""]
        .groupby(gene_columns, sort=False)["product_key"]
        .nunique()
    )

    conflicting_genes = product_counts[product_counts > 1]

    if not conflicting_genes.empty:
        examples = [
            f"{orthogroup}:{species}:{gene_id}"
            for orthogroup, species, gene_id in conflicting_genes.index[:10]
        ]

        raise ValueError(
            f"{len(conflicting_genes)} genes have several different product annotations. "
            f"Examples: {', '.join(examples)}"
        )

    # Prefer a non-empty product if duplicate GO rows contain empty values.
    annotations["has_product"] = annotations["product_key"] != ""

    genes = (
        annotations
        .sort_values("has_product", ascending=False)
        .drop_duplicates(subset=gene_columns, keep="first")
        .drop(columns="has_product")
        .reset_index(drop=True)
    )

    return genes


def choose_display_product(product_group):
    """Choose a readable label for a normalized product group."""
    counts = product_group["product"].value_counts()

    return sorted(
        counts[counts == counts.max()].index,
        key=str.casefold,
    )[0]


def calculate_copy_number_statistics(orthogroup_genes):
    """Calculate copy-number statistics across species containing the family."""
    copy_numbers = orthogroup_genes.groupby("specie_name", sort=False).size()

    if copy_numbers.empty:
        return {
            "min_copy_number": 0,
            "max_copy_number": 0,
            "max_copy_number_species": "",
            "copy_number_ratio": 0.0,
        }

    min_copy_number = int(copy_numbers.min())
    max_copy_number = int(copy_numbers.max())

    max_copy_number_species = "; ".join(
        sorted(
            copy_numbers[copy_numbers == max_copy_number].index,
            key=str.casefold,
        )
    )

    return {
        "min_copy_number": min_copy_number,
        "max_copy_number": max_copy_number,
        "max_copy_number_species": max_copy_number_species,
        "copy_number_ratio": max_copy_number / min_copy_number,
    }


def summarize_products(orthogroup_genes, total_species_count):
    """Calculate functional annotation, presence and copy-number statistics."""
    total_gene_count = len(orthogroup_genes)
    species_count = orthogroup_genes["specie_name"].nunique()
    species_fraction = (
        species_count / total_species_count
        if total_species_count
        else 0.0
    )
    copy_number_statistics = calculate_copy_number_statistics(orthogroup_genes)

    annotated_genes = orthogroup_genes[
        orthogroup_genes["product_key"] != ""
    ].copy()

    annotated_gene_count = len(annotated_genes)
    product_diversity = annotated_genes["product_key"].nunique()

    product_rows = []

    for key, group in annotated_genes.groupby("product_key", sort=False):
        product_rows.append(
            {
                "product_key": key,
                "product": choose_display_product(group),
                "count": len(group),
                "species_count": group["specie_name"].nunique(),
            }
        )

    product_summary = pd.DataFrame(
        product_rows,
        columns=["product_key", "product", "count", "species_count"],
    )

    if product_summary.empty:
        return {
            "representative_product": "",
            "representative_product_count": 0,
            "representative_product_fraction": 0.0,
            "annotated_gene_count": annotated_gene_count,
            "total_gene_count": total_gene_count,
            "product_diversity": product_diversity,
            "species_count": species_count,
            "species_fraction": round(species_fraction, 4),
            "min_copy_number": copy_number_statistics["min_copy_number"],
            "max_copy_number": copy_number_statistics["max_copy_number"],
            "max_copy_number_species": copy_number_statistics["max_copy_number_species"],
            "copy_number_ratio": round(copy_number_statistics["copy_number_ratio"], 4),
            "all_products": "",
        }

    product_summary = product_summary.sort_values(
        by=["count", "species_count", "product_key"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    informative_products = product_summary[
        ~product_summary["product_key"].isin(UNINFORMATIVE_PRODUCTS)
    ]

    if informative_products.empty:
        representative_product = ""
        representative_product_count = 0
        representative_product_fraction = 0.0
    else:
        representative = informative_products.iloc[0]
        representative_product = representative["product"]
        representative_product_count = int(representative["count"])
        representative_product_fraction = (
            representative_product_count / annotated_gene_count
            if annotated_gene_count
            else 0.0
        )

    all_products = "; ".join(
        f"{row.product} [{row.count}]"
        for row in product_summary.itertuples()
    )

    return {
        "representative_product": representative_product,
        "representative_product_count": representative_product_count,
        "representative_product_fraction": round(representative_product_fraction, 4),
        "annotated_gene_count": annotated_gene_count,
        "total_gene_count": total_gene_count,
        "product_diversity": product_diversity,
        "species_count": species_count,
        "species_fraction": round(species_fraction, 4),
        "min_copy_number": copy_number_statistics["min_copy_number"],
        "max_copy_number": copy_number_statistics["max_copy_number"],
        "max_copy_number_species": copy_number_statistics["max_copy_number_species"],
        "copy_number_ratio": round(copy_number_statistics["copy_number_ratio"], 4),
        "all_products": all_products,
    }


def annotate_orthogroups(genes, total_species_count):
    """Generate one functional annotation row per orthogroup."""
    rows = []

    for orthogroup, orthogroup_genes in genes.groupby("orthogroup", sort=True):
        row = {
            "orthogroup": orthogroup,
            **summarize_products(
                orthogroup_genes,
                total_species_count,
            ),
        }
        rows.append(row)

    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def main():
    args = parse_args()

    annotations = read_annotations(args.annotations)
    total_species_count = annotations["specie_name"].nunique()

    genes = select_gene_products(annotations)
    results = annotate_orthogroups(
        genes,
        total_species_count,
    )

    results.to_csv(
        args.output,
        sep="\t",
        index=False,
        float_format="%.4f",
    )

    print(f"Read {len(args.annotations)} annotation tables")
    print(f"Found {total_species_count} species")
    print(f"Retained {len(genes)} genes assigned to orthogroups")
    print(f"Annotated {len(results)} orthogroups")
    print(f"Families with at least 100 copies in one species: {(results['max_copy_number'] >= 100).sum()}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()