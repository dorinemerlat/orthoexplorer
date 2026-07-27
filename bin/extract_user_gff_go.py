#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import gffutils
import pandas as pd


TRANSCRIPT_FEATURE_TYPES = ("mRNA", "transcript")
GO_PATTERN = re.compile(r"GO:\d{7}")

OUTPUT_COLUMNS = [
    "specie_name",
    "id",
    "gene_id",
    "gene_name",
    "product",
    "protein_id",
    "go_id",
    "go_term",
    "category",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Extract GO annotations from user-provided GFF3 transcripts.")
    parser.add_argument("--gff", required=True, type=Path, help="Input GFF3 file.")
    parser.add_argument("--species-name", required=True, help="Scientific species name.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV file.")
    return parser.parse_args()


def get_attribute(feature, name, default=""):
    """Return the first value of a GFF3 attribute."""
    values = feature.attributes.get(name, [])
    return values[0] if values else default


def extract_go_terms(feature):
    """Return GO identifiers associated with a transcript."""
    go_terms = set()

    for value in feature.attributes.get("Ontology_term", []):
        go_terms.update(GO_PATTERN.findall(value))

    return sorted(go_terms)


def extract_annotations(gff_path, species_name):
    """Extract one row per transcript and GO term."""

    # Build an in-memory GFF database for easy feature traversal.
    db = gffutils.create_db(
        gff_path,
        ":memory:",
        merge_strategy="create_unique",
        keep_order=True,
    )

    rows = []

    # Only transcript features are expected to carry functional annotations.
    for feature_type in TRANSCRIPT_FEATURE_TYPES:
        for transcript in db.features_of_type(feature_type, order_by=("seqid", "start")):
            protein_id = get_attribute(transcript, "ID", transcript.id)
            gene_id = get_attribute(transcript, "Parent")
            gene_name = get_attribute(transcript, "Name")
            product = get_attribute(transcript, "product")
            go_ids = extract_go_terms(transcript)

            # Write one output row per GO term.
            for go_id in go_ids:
                rows.append(
                    {
                        "specie_name": species_name,
                        "id": protein_id,
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "product": product,
                        "protein_id": protein_id,
                        "go_id": go_id,
                        "go_term": "",
                        "category": "",
                    }
                )

    return rows


def main():
    args = parse_args()

    rows = extract_annotations(args.gff, args.species_name)

    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(
        args.output,
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()