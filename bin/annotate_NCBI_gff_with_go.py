#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import gffutils
import pandas as pd


TRANSCRIPT_FEATURE_TYPES = ("mRNA", "transcript")


def parse_args():
    parser = argparse.ArgumentParser(description="Add NCBI Gene Ontology annotations to GFF3 transcripts.")
    parser.add_argument("--gff", required=True, type=Path, help="Input NCBI GFF3 file.")
    parser.add_argument("--gene2go", required=True, type=Path, help="NCBI gene2go file, optionally gzipped.")
    parser.add_argument("--taxid", required=True, help="NCBI taxonomy identifier.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV file.")
    return parser.parse_args()


def get_attribute(feature, name, default=""):
    """Return a GFF3 attribute value."""
    return feature.attributes.get(name, [default])[0]


def extract_gene_id(feature):
    """Return the NCBI GeneID from the Dbxref attribute."""
    for dbxref in feature.attributes.get("Dbxref", []):
        if dbxref.startswith("GeneID:"):
            return dbxref.removeprefix("GeneID:")
    return ""


def keep_most_complete_rows(df, subset, columns):
    """Keep the most complete row for each group."""
    df = df.copy()

    df["_score"] = df[columns].ne("").sum(axis=1)

    df = (
        df.sort_values("_score", ascending=False)
        .drop_duplicates(subset=subset, keep="first")
        .drop(columns="_score")
        .reset_index(drop=True)
    )

    return df


def read_gff_transcripts(gff_path):
    """Read transcript annotations from an NCBI GFF3 file."""
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
                    "gene_id": extract_gene_id(transcript),
                    "gene": get_attribute(transcript, "gene"),
                    "product": get_attribute(transcript, "product"),
                    "protein_id": get_attribute(transcript, "protein_id"),
                }
            )

    df = pd.DataFrame(
        rows,
        columns=["gene_id", "gene", "product", "protein_id"],
    )

    df = keep_most_complete_rows(
        df,
        subset="gene_id",
        columns=["gene", "product", "protein_id"],
    )

    return df


def load_taxid_annotations(gene2go_path, taxid):
    """Load all gene2go annotations for the requested taxon."""
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

    for chunk in pd.read_csv(
        gene2go_path,
        sep="\t",
        names=columns,
        header=0,
        dtype=str,
        compression="infer",
        chunksize=500_000,
    ):
        chunk = chunk[chunk["taxid"] == str(taxid)]

        if not chunk.empty:
            chunks.append(
                chunk[["taxid", "gene_id", "go_id", "go_term", "category"]]
            )

    if not chunks:
        return pd.DataFrame(
            columns=["taxid", "gene_id", "go_id", "go_term", "category"]
        )

    return pd.concat(chunks, ignore_index=True)


def load_go_terms(gene2go_taxid, requested_gene_ids):
    """Extract GO terms for the requested GeneIDs."""

    chunk = gene2go_taxid[
        gene2go_taxid["gene_id"].isin(requested_gene_ids)
        & gene2go_taxid["go_id"].str.fullmatch(r"GO:\d{7}", na=False)
    ].copy()

    qualifiers = chunk["qualifier"].fillna("").str.upper()
    chunk = chunk[
        ~qualifiers.str.split("|", regex=False).apply(lambda values: "NOT" in values)
    ]

    return (
        chunk[["gene_id", "go_id"]]
        .drop_duplicates()
        .groupby("gene_id", as_index=False)["go_id"]
        .agg(lambda values: ",".join(sorted(values)))
        .rename(columns={"go_id": "go_terms"})
    )


def main():
    args = parse_args()

    transcripts = read_gff_transcripts(args.gff)
    print(f"Read {len(transcripts)} transcripts from GFF: {args.gff}")

    requested_gene_ids = transcripts["gene_id"].tolist()

    gene2go_taxid = load_taxid_annotations(args.gene2go, args.taxid)
    print(f"Loaded {len(gene2go_taxid)} gene2go annotations for taxid {args.taxid}")

    go_annotations = transcripts.merge(
        gene2go_taxid,
        on="gene_id",
        how="left",
    )
    
    go_annotations.to_csv(args.output, sep="\t", index=False)
    print(f"GO annotations written to: {args.output}")

if __name__ == "__main__":
    main()