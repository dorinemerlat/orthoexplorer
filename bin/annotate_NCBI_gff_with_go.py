#!/usr/bin/env python3

import argparse
from pathlib import Path

import gffutils
import pandas as pd


TRANSCRIPT_FEATURE_TYPES = ("mRNA", "transcript")

OUTPUT_COLUMNS = [
    "specie_name",
    "id",
    "protein_id",
    "gene_id",
    "gene_name",
    "product",
    "go_id",
    "go_term",
    "category",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Extract NCBI GO annotations for transcripts from a GFF3 file.")
    parser.add_argument("--gff", required=True, type=Path, help="Input NCBI GFF3 file.")
    parser.add_argument("--gene2go", required=True, type=Path, help="NCBI gene2go file, optionally gzipped.")
    parser.add_argument("--taxid", required=True, help="NCBI taxonomy identifier.")
    parser.add_argument("--species-name", required=True, help="Scientific species name.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV file.")
    return parser.parse_args()


def get_attribute(feature, name, default=""):
    """Return the first value of a GFF3 attribute."""
    values = feature.attributes.get(name, [])
    return values[0] if values else default


def extract_gene_id(feature):
    """Return the NCBI GeneID from the Dbxref attribute."""
    for dbxref in feature.attributes.get("Dbxref", []):
        if dbxref.startswith("GeneID:"):
            return dbxref.removeprefix("GeneID:")

    return ""


def keep_most_complete_rows(df, subset, columns):
    """Keep the most complete row for each identifier."""
    df = df.copy()
    df["_score"] = df[columns].ne("").sum(axis=1)

    return (
        df.sort_values("_score", ascending=False)
        .drop_duplicates(subset=subset, keep="first")
        .drop(columns="_score")
        .reset_index(drop=True)
    )


def read_gff_transcripts(gff_path):
    """Read transcript annotations from an NCBI GFF3 file."""

    # Build an in-memory database to query transcript features.
    db = gffutils.create_db(
        str(gff_path),
        ":memory:",
        merge_strategy="create_unique",
        keep_order=True,
    )

    rows = []

    for feature_type in TRANSCRIPT_FEATURE_TYPES:
        for transcript in db.features_of_type(feature_type, order_by=("seqid", "start")):
            rows.append(
                {
                    "id": get_attribute(transcript, "ID", transcript.id),
                    "protein_id": get_attribute(transcript, "protein_id"),
                    "gene_id": extract_gene_id(transcript),
                    "gene_name": get_attribute(transcript, "gene"),
                    "product": get_attribute(transcript, "product"),
                }
            )

    transcripts = pd.DataFrame(
        rows,
        columns=["id", "protein_id", "gene_id", "gene_name", "product"],
    )

    # NCBI GFF files may contain several transcripts for the same GeneID.
    return keep_most_complete_rows(
        transcripts,
        subset="gene_id",
        columns=["id", "protein_id", "gene_name", "product"],
    )


def load_taxid_annotations(gene2go_path, taxid, requested_gene_ids):
    """Load GO annotations for the requested taxon and GeneIDs."""
    columns = [
        "taxid",
        "gene_id",
        "go_id",
        "evidence",
        "qualifier",
        "go_term",
        "pubmed",
        "category",
    ]

    chunks = []
    requested_gene_ids = set(requested_gene_ids)

    # Read gene2go in chunks because the complete file can be large.
    for chunk in pd.read_csv(
        gene2go_path,
        sep="\t",
        names=columns,
        header=0,
        dtype=str,
        compression="infer",
        chunksize=500_000,
    ):
        chunk = chunk[
            (chunk["taxid"] == str(taxid))
            & chunk["gene_id"].isin(requested_gene_ids)
            & chunk["go_id"].str.fullmatch(r"GO:\d{7}", na=False)
        ].copy()

        if chunk.empty:
            continue

        # Exclude annotations explicitly negated with the NOT qualifier.
        qualifiers = chunk["qualifier"].fillna("").str.upper()
        chunk = chunk[
            ~qualifiers.str.split("|", regex=False).apply(lambda values: "NOT" in values)
        ]

        if not chunk.empty:
            chunks.append(chunk[["gene_id", "go_id", "go_term", "category"]])

    if not chunks:
        return pd.DataFrame(columns=["gene_id", "go_id", "go_term", "category"])

    return (
        pd.concat(chunks, ignore_index=True)
        .drop_duplicates(subset=["gene_id", "go_id", "go_term", "category"])
        .reset_index(drop=True)
    )


def main():
    args = parse_args()

    transcripts = read_gff_transcripts(args.gff)
    print(f"Read {len(transcripts)} genes from GFF: {args.gff}")

    gene2go_annotations = load_taxid_annotations(
        args.gene2go,
        args.taxid,
        transcripts["gene_id"].dropna(),
    )
    print(f"Loaded {len(gene2go_annotations)} GO annotations for taxid {args.taxid}")

    # Keep every GFF gene, including genes without a gene2go annotation.
    annotations = transcripts.merge(
        gene2go_annotations,
        on="gene_id",
        how="left",
    )

    annotations.insert(0, "specie_name", args.species_name)
    annotations = annotations[OUTPUT_COLUMNS].fillna("")

    annotations.to_csv(args.output, sep="\t", index=False)
    print(f"GO annotations written to: {args.output}")


if __name__ == "__main__":
    main()
