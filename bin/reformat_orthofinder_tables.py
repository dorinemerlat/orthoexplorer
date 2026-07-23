#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
import re
import shutil
from typing import Dict, Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rename OrthoFinder species columns using the species and "
            "outgroup CSV tables. Only table headers are modified."
        )
    )
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--gene-counts", required=True, type=Path)
    parser.add_argument("--unassigned", type=Path)
    parser.add_argument("--species", required=True, type=Path)
    parser.add_argument("--outgroups", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser.parse_args()


def slugify(name: str) -> str:
    slug = re.sub(r"\s+", "-", name.strip())
    slug = slug.replace(".", "-")
    slug = re.sub(r"[^A-Za-z0-9_-]", "-", slug)
    slug = re.sub(r"-+", "-", slug)
    return slug.strip("-").lower()


def fasta_stem(path_value: str) -> str:
    name = Path(path_value).name
    lower = name.lower()

    for suffix in (
        ".fasta.gz", ".faa.gz", ".fa.gz", ".fas.gz",
        ".fasta", ".faa", ".fa", ".fas",
    ):
        if lower.endswith(suffix):
            return name[: -len(suffix)]

    return Path(name).stem


def read_mapping(csv_paths: Iterable[Path]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}

    for csv_path in csv_paths:
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)

            required = {"taxid", "name"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"{csv_path} is missing required columns: "
                    f"{', '.join(sorted(missing))}"
                )

            for line_number, row in enumerate(reader, start=2):
                taxid = (row.get("taxid") or "").strip()
                name = (row.get("name") or "").strip()
                supplied_file = (row.get("file") or "").strip()

                if not taxid or not name:
                    raise ValueError(
                        f"{csv_path}:{line_number}: taxid and name are required"
                    )

                source_name = (
                    fasta_stem(supplied_file)
                    if supplied_file
                    else slugify(name)
                )

                previous = mapping.get(source_name)
                if previous and previous != name:
                    raise ValueError(
                        f"OrthoFinder name '{source_name}' maps to both "
                        f"'{previous}' and '{name}'"
                    )

                mapping[source_name] = name

    return mapping


def rename_header(
    input_path: Path,
    output_path: Path,
    mapping: Dict[str, str],
) -> None:
    with input_path.open(newline="", encoding="utf-8") as source:
        reader = csv.reader(source, delimiter="\t")
        rows = iter(reader)

        try:
            header = next(rows)
        except StopIteration:
            output_path.touch()
            return

        renamed_header = [mapping.get(column, column) for column in header]

        with output_path.open("w", newline="", encoding="utf-8") as target:
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(renamed_header)
            writer.writerows(rows)


def output_path(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}.{suffix}")


def main() -> None:
    args = parse_args()
    mapping = read_mapping((args.species, args.outgroups))

    rename_header(
        args.orthogroups,
        output_path(args.output_prefix, "Orthogroups.tsv"),
        mapping,
    )
    rename_header(
        args.gene_counts,
        output_path(args.output_prefix, "Orthogroups.GeneCount.tsv"),
        mapping,
    )

    if args.unassigned:
        rename_header(
            args.unassigned,
            output_path(
                args.output_prefix,
                "Orthogroups_UnassignedGenes.tsv",
            ),
            mapping,
        )


if __name__ == "__main__":
    main()
