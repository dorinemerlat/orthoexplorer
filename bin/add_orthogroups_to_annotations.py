#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Add OrthoFinder orthogroups to a GO annotation table.")
    parser.add_argument("--annotations", required=True, type=Path, help="Input GO annotation TSV file.")
    parser.add_argument("--orthogroups", required=True, type=Path, help="OrthoFinder Orthogroups.tsv file.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV file.")
    return parser.parse_args()


def normalize_species_name(name):
    """Convert a species name to the identifier format used by OrthoFinder."""
    name = str(name).strip().lower()
    name = re.sub(r"[^a-z0-9]+", "-", name)
    return name.strip("-")


def read_annotations(annotation_path):
    """Read the annotation table and validate required columns."""
    annotations = pd.read_csv(annotation_path, sep="\t", dtype=str, keep_default_na=False)

    required_columns = {"specie_name", "id"}
    missing_columns = required_columns - set(annotations.columns)

    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Missing required columns in {annotation_path}: {missing}")

    if "orthogroup" in annotations.columns:
        raise ValueError(f"Column 'orthogroup' already exists in annotation file: {annotation_path}")

    annotations["specie"] = annotations["specie_name"].map(normalize_species_name)

    return annotations


def reshape_orthogroups(orthogroup_path):
    """Convert Orthogroups.tsv to one row per species, protein and orthogroup."""
    orthogroups = pd.read_csv(orthogroup_path, sep="\t", dtype=str)

    if "Orthogroup" not in orthogroups.columns:
        raise ValueError(f"Missing required column 'Orthogroup' in file: {orthogroup_path}")

    # Convert species columns into a long-format table.
    orthogroups = orthogroups.melt(
        id_vars="Orthogroup",
        var_name="specie",
        value_name="id",
    )

    # Create one row per comma-separated protein identifier.
    orthogroups["id"] = orthogroups["id"].fillna("").str.split(",")
    orthogroups = orthogroups.explode("id")
    orthogroups["id"] = orthogroups["id"].str.strip()
    orthogroups["specie"] = orthogroups["specie"].map(normalize_species_name)

    orthogroups = (
        orthogroups[orthogroups["id"] != ""]
        .rename(columns={"Orthogroup": "orthogroup"})
        [["orthogroup", "specie", "id"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # A protein identifier must belong to only one orthogroup within one species.
    duplicated = orthogroups.duplicated(
        subset=["specie", "id"],
        keep=False,
    )

    if duplicated.any():
        examples = (
            orthogroups.loc[duplicated, ["specie", "id"]]
            .drop_duplicates()
            .head(10)
            .apply(lambda row: f"{row['specie']}:{row['id']}", axis=1)
            .tolist()
        )

        raise ValueError(
            f"{duplicated.sum()} species/protein combinations belong to several orthogroups. "
            f"Examples: {', '.join(examples)}"
        )

    return orthogroups


def add_orthogroups(annotations, orthogroups):
    """Add orthogroups using the species name and protein identifier."""
    annotations = annotations.merge(
        orthogroups,
        on=["specie", "id"],
        how="left",
        validate="many_to_one",
    )

    annotations["orthogroup"] = annotations["orthogroup"].fillna("")
    annotations = annotations.drop(columns="specie")

    # Place orthogroup directly after specie_name.
    columns = annotations.columns.tolist()
    columns.remove("orthogroup")
    columns.insert(columns.index("specie_name") + 1, "orthogroup")

    return annotations[columns]


def main():
    args = parse_args()

    annotations = read_annotations(args.annotations)
    orthogroups = reshape_orthogroups(args.orthogroups)
    annotations = add_orthogroups(annotations, orthogroups)

    annotated_count = annotations["orthogroup"].ne("").sum()
    missing_count = annotations["orthogroup"].eq("").sum()

    annotations.to_csv(args.output, sep="\t", index=False)

    print(f"Read {len(orthogroups)} protein-to-orthogroup associations")
    print(f"Assigned an orthogroup to {annotated_count} annotation rows")
    print(f"No orthogroup found for {missing_count} annotation rows")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()