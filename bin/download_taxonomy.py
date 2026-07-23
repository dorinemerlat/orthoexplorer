#!/usr/bin/env python3

import argparse
import csv
import time
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


def fetch_url(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "taxonomy-fetcher/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read()


def is_taxid(value):
    return value.strip().isdigit()


def resolve_name_to_taxid(name):
    query = urllib.parse.quote(name)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=taxonomy&term={query}&retmode=xml"
    )

    data = fetch_url(url)
    root = ET.fromstring(data)

    ids = [x.text for x in root.findall(".//IdList/Id")]

    if not ids:
        raise ValueError(f"No taxid found for name: {name}")

    return ids[0]


def fetch_taxonomy(taxid):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=taxonomy&id={taxid}&retmode=xml"
    )

    data = fetch_url(url)
    root = ET.fromstring(data)

    taxon = root.find(".//Taxon")
    if taxon is None:
        raise ValueError(f"No taxonomy found for taxid: {taxid}")

    scientific_name = taxon.findtext("ScientificName", default="")
    taxid = taxon.findtext("TaxId", default=taxid)

    lineage = {}

    for lineage_taxon in taxon.findall(".//LineageEx/Taxon"):
        rank = lineage_taxon.findtext("Rank", default="").strip()
        name = lineage_taxon.findtext("ScientificName", default="").strip()

        if not rank:
            continue

        rank_lower = rank.lower()

        if rank_lower in {"no rank", "clade"}:
            continue

        lineage[rank_lower] = name

    current_rank = taxon.findtext("Rank", default="").strip().lower()
    if current_rank and current_rank not in {"no rank", "clade"}:
        lineage[current_rank] = scientific_name

    return {
        "specie": scientific_name,
        "taxid": taxid,
        "lineage": lineage,
    }


def read_input(input_file):
    values = []

    with open(input_file, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.lower() in {"species", "specie", "name", "taxid"}:
                continue

            values.append(line)

    return values


def main():
    parser = argparse.ArgumentParser(
        description="Fetch taxonomy table from taxids or species names using NCBI Taxonomy API."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file with one taxid or species name per line."
    )

    parser.add_argument(
        "-o", "--output",
        default="taxonomy_table.tsv",
        help="Output TSV file."
    )

    parser.add_argument(
        "--sleep",
        type=float,
        default=1,
        help="Delay between API calls, to avoid hammering NCBI."
    )

    args = parser.parse_args()

    queries = read_input(args.input)
    results = []

    for query in queries:
        query = query.strip()

        try:
            if is_taxid(query):
                taxid = query
            else:
                taxid = resolve_name_to_taxid(query)
                time.sleep(args.sleep)

            taxonomy = fetch_taxonomy(taxid)
            results.append(taxonomy)

            print(f"OK: {query} -> {taxonomy['specie']} ({taxonomy['taxid']})")

        except Exception as e:
            print(f"WARNING: failed for {query}: {e}")
            results.append({
                "specie": query,
                "taxid": "NA",
                "lineage": {},
            })

        time.sleep(args.sleep)

    used_ranks = []
    for rank in RANK_ORDER:
        if any(rank in result["lineage"] for result in results):
            used_ranks.append(rank)

    fieldnames = ["specie", "taxid"] + used_ranks

    with open(args.output, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for result in results:
            row = {
                "specie": result["specie"],
                "taxid": result["taxid"],
            }

            for rank in used_ranks:
                row[rank] = result["lineage"].get(rank, "")

            writer.writerow(row)

    print(f"\nWritten: {args.output}")


if __name__ == "__main__":
    main()