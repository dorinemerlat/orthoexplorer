#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Build orthogroup-to-GO associations from annotated orthogroup tables.")
    parser.add_argument("--annotations", required=True, nargs="+", type=Path, help="Input GO annotation TSV files.")
    parser.add_argument("--output", required=True, type=Path, help="Output orthogroup-to-GO TSV file.")
    return parser.parse_args()


def main():
    args = parse_args()

    annotation_tables = []

    for annotation_file in args.annotations:
        table = pd.read_csv(annotation_file, sep="\t", dtype=str, keep_default_na=False)

        required_columns = {"orthogroup", "go_id"}
        missing_columns = required_columns - set(table.columns)

        if missing_columns:
            raise ValueError(f"{annotation_file}: missing columns: {', '.join(sorted(missing_columns))}")

        annotation_tables.append(table[["orthogroup", "go_id"]])

    annotations = pd.concat(annotation_tables, ignore_index=True)

    annotations = annotations[
        annotations["orthogroup"].str.startswith("OG")
        & annotations["go_id"].str.match(r"^GO:\d{7}$")
    ].drop_duplicates()

    orthogroup2go = (
        annotations.groupby("orthogroup")["go_id"]
        .apply(lambda go_ids: ";".join(sorted(set(go_ids))))
        .reset_index()
    )

    orthogroup2go.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()