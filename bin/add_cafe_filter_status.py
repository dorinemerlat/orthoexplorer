#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "representative_product",
    "copy_number_ratio",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Add CAFE filtering recommendations to an orthogroup annotation table.")
    parser.add_argument("--annotations", required=True, type=Path, help="Input orthogroup annotation TSV file.")
    parser.add_argument("--blacklist", required=True, type=Path, help="Text file containing one blacklist term per line.")
    parser.add_argument("--ratio-threshold", required=True, type=float, help="Copy-number ratio above which families are evaluated.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV file.")
    return parser.parse_args()


def read_blacklist(path):
    """Read non-empty, uncommented blacklist terms."""
    with path.open() as handle:
        terms = [
            line.strip().lower()
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]

    if not terms:
        raise ValueError(f"No blacklist terms were found in '{path}'.")

    return terms


def validate_columns(table):
    """Check that the required columns are present."""
    missing_columns = REQUIRED_COLUMNS - set(table.columns)

    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Missing required column(s): {missing}")


def find_matching_keyword(product, blacklist):
    """Return the first blacklist term found in the product name."""
    if pd.isna(product):
        return ""

    normalized_product = str(product).strip().lower()

    for keyword in blacklist:
        if keyword in normalized_product:
            return keyword

    return ""


def assign_cafe_status(row, ratio_threshold):
    """Assign the recommended CAFE filtering status."""
    ratio = row["copy_number_ratio"]
    matched_keyword = row["cafe_matched_keyword"]

    if pd.isna(ratio):
        return pd.Series(
            {
                "cafe_status": "suspect",
                "cafe_filter_reason": "missing_copy_number_ratio",
            }
        )

    if ratio <= ratio_threshold:
        return pd.Series(
            {
                "cafe_status": "keep",
                "cafe_filter_reason": "",
            }
        )

    if matched_keyword:
        return pd.Series(
            {
                "cafe_status": "remove",
                "cafe_filter_reason": "high_copy_number_ratio_and_blacklisted_product",
            }
        )

    return pd.Series(
        {
            "cafe_status": "suspect",
            "cafe_filter_reason": "high_copy_number_ratio",
        }
    )


def main():
    args = parse_args()

    if args.ratio_threshold <= 0:
        raise ValueError("--ratio-threshold must be greater than 0.")

    blacklist = read_blacklist(args.blacklist)
    table = pd.read_csv(args.annotations, sep="\t")
    validate_columns(table)

    table["copy_number_ratio"] = pd.to_numeric(table["copy_number_ratio"], errors="coerce")
    table["cafe_matched_keyword"] = table["representative_product"].apply(
        find_matching_keyword,
        blacklist=blacklist,
    )

    statuses = table.apply(assign_cafe_status, axis=1, ratio_threshold=args.ratio_threshold)
    table["cafe_status"] = statuses["cafe_status"]
    table["cafe_filter_reason"] = statuses["cafe_filter_reason"]

    table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()