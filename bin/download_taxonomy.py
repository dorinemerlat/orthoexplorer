#!/usr/bin/env python3

import argparse
import csv
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


RANK_ORDER = [
    "superkingdom",
    "kingdom",
    "phylum",
    "subphylum",
    "class",
    "subclass",
    "infraclass",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "genus",
    "species",
    "subspecies",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Fetch NCBI taxonomy information for datasets listed in a CSV file."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=(
            "Input CSV file containing the columns: "
            "id, name and taxid."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        default="taxonomy.tsv",
        help="Output taxonomy TSV file.",
    )

    parser.add_argument(
        "--sleep",
        type=float,
        default=1.0,
        help="Delay in seconds between NCBI API requests.",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="HTTP request timeout in seconds.",
    )

    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Maximum number of attempts for each HTTP request.",
    )

    return parser.parse_args()


def fetch_url(url, timeout=30, retries=3):
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "OrthoExplorer-taxonomy-fetcher/1.0",
        },
    )

    last_error = None

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(
                request,
                timeout=timeout,
            ) as response:
                return response.read()

        except (
            urllib.error.HTTPError,
            urllib.error.URLError,
            TimeoutError,
        ) as error:
            last_error = error

            if attempt < retries:
                time.sleep(attempt * 2)

    raise RuntimeError(
        f"NCBI request failed after {retries} attempts: {last_error}"
    )


def fetch_taxonomy(taxid, timeout=30, retries=3):
    encoded_taxid = urllib.parse.quote(str(taxid))

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        f"efetch.fcgi?db=taxonomy&id={encoded_taxid}&retmode=xml"
    )

    data = fetch_url(
        url,
        timeout=timeout,
        retries=retries,
    )

    root = ET.fromstring(data)
    taxon = root.find(".//Taxon")

    if taxon is None:
        raise ValueError(
            f"No taxonomy found for taxid: {taxid}"
        )

    scientific_name = taxon.findtext(
        "ScientificName",
        default="",
    ).strip()

    resolved_taxid = taxon.findtext(
        "TaxId",
        default=str(taxid),
    ).strip()

    lineage = {}

    for lineage_taxon in taxon.findall(".//LineageEx/Taxon"):
        rank = lineage_taxon.findtext(
            "Rank",
            default="",
        ).strip().lower()

        scientific_lineage_name = lineage_taxon.findtext(
            "ScientificName",
            default="",
        ).strip()

        if not rank or rank in {"no rank", "clade"}:
            continue

        lineage[rank] = scientific_lineage_name

    current_rank = taxon.findtext(
        "Rank",
        default="",
    ).strip().lower()

    if current_rank and current_rank not in {"no rank", "clade"}:
        lineage[current_rank] = scientific_name

    return {
        "specie_name": scientific_name,
        "taxid": resolved_taxid,
        "lineage": lineage,
    }


def read_input_csv(input_file):
    datasets = []
    seen_ids = set()

    with open(
        input_file,
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)

        required_columns = {
            "id",
            "name",
            "taxid",
            "group",
        }

        missing_columns = required_columns - set(
            reader.fieldnames or []
        )

        if missing_columns:
            raise ValueError(
                "Input CSV is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for line_number, row in enumerate(reader, start=2):
            dataset_id = (row.get("id") or "").strip()
            name = (row.get("name") or "").strip()
            taxid = (row.get("taxid") or "").strip()
            group = (row.get("group") or "").strip().lower()

            if not dataset_id:
                raise ValueError(
                    f"Line {line_number}: missing dataset ID"
                )

            if dataset_id in seen_ids:
                raise ValueError(
                    f"Line {line_number}: duplicate dataset ID "
                    f"'{dataset_id}'"
                )

            if not name:
                raise ValueError(
                    f"Line {line_number}: missing dataset name"
                )

            if not taxid:
                raise ValueError(
                    f"Line {line_number}: missing taxid for "
                    f"'{dataset_id}'"
                )

            if not taxid.isdigit():
                raise ValueError(
                    f"Line {line_number}: invalid taxid "
                    f"'{taxid}' for '{dataset_id}'"
                )

            if group not in {"ingroup", "outgroup"}:
                raise ValueError(
                    f"Line {line_number}: invalid group "
                    f"'{group}' for '{dataset_id}'. "
                    "Expected 'ingroup' or 'outgroup'."
                )

            seen_ids.add(dataset_id)

            datasets.append(
                {
                    "id": dataset_id,
                    "name": name,
                    "taxid": taxid,
                    "group": group,
                }
            )

    return datasets


def determine_used_ranks(results):
    return [
        rank
        for rank in RANK_ORDER
        if any(
            rank in result["lineage"]
            for result in results
        )
    ]


def write_output(output_file, results):
    used_ranks = determine_used_ranks(results)

    fieldnames = [
        "id",
        "name",
        "specie_name",
        "taxid",
        "group",
        *used_ranks,
    ]

    with open(
        output_file,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for result in results:
            row = {
                "id": result["id"],
                "name": result["name"],
                "specie_name": result["specie_name"],
                "taxid": result["taxid"],
                "group": result["group"],
            }

            for rank in used_ranks:
                row[rank] = result["lineage"].get(
                    rank,
                    "",
                )

            writer.writerow(row)


def main():
    args = parse_args()

    try:
        datasets = read_input_csv(args.input)
    except (OSError, ValueError) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1

    results = []
    failed_queries = 0

    for index, dataset in enumerate(datasets):
        dataset_id = dataset["id"]
        input_name = dataset["name"]
        input_taxid = dataset["taxid"]
        input_group = dataset["group"]

        try:
            taxonomy = fetch_taxonomy(
                input_taxid,
                timeout=args.timeout,
                retries=args.retries,
            )

            results.append(
                {
                    "id": dataset_id,
                    "name": input_name,
                    "specie_name": taxonomy["specie_name"],
                    "taxid": taxonomy["taxid"],
                    "group": input_group,
                    "lineage": taxonomy["lineage"],
                }
            )

            print(
                f"OK: {dataset_id} -> "
                f"{taxonomy['specie_name']} "
                f"({taxonomy['taxid']})"
            )

        except Exception as error:
            failed_queries += 1

            print(
                f"WARNING: failed for {dataset_id} "
                f"(taxid {input_taxid}): {error}",
                file=sys.stderr,
            )

            results.append(
                {
                    "id": dataset_id,
                    "name": input_name,
                    "specie_name": "",
                    "taxid": input_taxid,
                    "group": input_group,
                    "lineage": {},
                }
            )

        if index < len(datasets) - 1:
            time.sleep(args.sleep)

    write_output(
        args.output,
        results,
    )

    print(f"Written: {args.output}")

    if failed_queries:
        print(
            f"WARNING: taxonomy retrieval failed for "
            f"{failed_queries} dataset(s).",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())