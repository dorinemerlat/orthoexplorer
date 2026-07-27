#!/usr/bin/env python3

import argparse
import csv
import re
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml
from Bio import Phylo
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter


DEFAULT_RANKS = ["phylum", "subphylum", "class", "subclass", "order"]


def parse_args():
    parser = argparse.ArgumentParser(description="Assign genes to conservation ranks from OrthoFinder results.")
    parser.add_argument("--gene-counts", required=True, type=Path, help="Orthogroups.GeneCount.tsv")
    parser.add_argument("--orthogroups", required=True, type=Path, help="Orthogroups.tsv")
    parser.add_argument("--unassigned", required=True, type=Path, help="Orthogroups_UnassignedGenes.tsv")
    parser.add_argument("--tree", required=True, type=Path, help="Newick species tree")
    parser.add_argument("--taxonomy", required=True, type=Path, help="Taxonomy TSV containing name, group and taxonomy ranks")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--ranks", default=",".join(DEFAULT_RANKS), help="Comma-separated taxonomy ranks to use")
    parser.add_argument("--threshold", type=float, default=75.0, help="Minimum conservation percentage used for rank assignment")
    parser.add_argument("--max-outgroup-rank", default="subphylum", help="Most specific rank tested for outgroup species")
    parser.add_argument("--colors", type=Path, help="Optional flat YAML file mapping clade names to colors")
    parser.add_argument("--clades", default="", help="Optional comma-separated clade names found among taxonomy values, for example Dignatha,Progoneata")
    return parser.parse_args()


def safe_name(value):
    """Return a filesystem-safe file name."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())


def split_genes(cell):
    """Split an OrthoFinder gene-list cell into gene identifiers."""
    if not isinstance(cell, str) or not cell.strip():
        return []
    return [gene.strip() for gene in cell.split(",") if gene.strip()]


def thousands(value, _position):
    """Format axis labels with spaces as thousands separators."""
    return f"{int(value):,}".replace(",", " ")


def read_table(path):
    """Read a tab-separated table as strings and replace missing values with empty strings."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def read_gene_counts(path):
    """Read Orthogroups.GeneCount.tsv and convert species columns to integers."""
    table = read_table(path)
    table = table.rename(columns={table.columns[0]: "Orthogroup"})
    table = table.drop(columns=["Total"], errors="ignore")
    species_columns = [column for column in table.columns if column != "Orthogroup"]

    for species in species_columns:
        table[species] = pd.to_numeric(table[species], errors="coerce").fillna(0).astype(int)

    return table, species_columns


def read_orthogroups(path):
    """Read Orthogroups.tsv or Orthogroups_UnassignedGenes.tsv."""
    table = read_table(path)
    table = table.rename(columns={table.columns[0]: "Orthogroup"})
    return table.drop(columns=["Total"], errors="ignore")


def read_taxonomy(path, ranks):
    """Read and validate the taxonomy table."""
    taxonomy = read_table(path)
    taxonomy.columns = [column.strip() for column in taxonomy.columns]

    required_columns = {"name", "group"}
    missing_columns = required_columns - set(taxonomy.columns)
    if missing_columns:
        raise ValueError(f"Taxonomy table is missing required columns: {', '.join(sorted(missing_columns))}")

    missing_ranks = [rank for rank in ranks if rank not in taxonomy.columns]
    if missing_ranks:
        raise ValueError(f"Taxonomy table is missing requested rank columns: {', '.join(missing_ranks)}")

    taxonomy["name"] = taxonomy["name"].astype(str).str.strip()
    taxonomy["group"] = taxonomy["group"].astype(str).str.strip().str.lower()

    invalid_groups = sorted(set(taxonomy["group"]) - {"ingroup", "outgroup"})
    if invalid_groups:
        raise ValueError(f"Invalid values in taxonomy group column: {', '.join(invalid_groups)}")

    duplicated_names = taxonomy.loc[taxonomy["name"].duplicated(keep=False), "name"].unique()
    if len(duplicated_names):
        raise ValueError(f"Duplicate names in taxonomy table: {', '.join(sorted(duplicated_names))}")

    return taxonomy

def read_tree_order(path):
    """Return terminal names in left-to-right Newick order."""
    newick = path.read_text(encoding="utf-8").strip()
    newick = newick.replace('""', "")
    newick = re.sub(r'"([^"]+)"', r"'\1'", newick)
    tree = Phylo.read(StringIO(newick), "newick")
    names = []

    for terminal in tree.get_terminals():
        if terminal.name is None:
            raise ValueError("The tree contains an unnamed terminal.")

        name = terminal.name.strip().strip("'\"")
        if name in names:
            raise ValueError(f"Duplicate species in tree: {name}")

        names.append(name)

    return names

def validate_species(tree_order, taxonomy, species_columns):
    """Ensure OrthoFinder species are represented in both the tree and taxonomy table."""
    taxonomy_names = set(taxonomy["name"])
    tree_names = set(tree_order)
    species_set = set(species_columns)

    missing_from_taxonomy = sorted(species_set - taxonomy_names)
    if missing_from_taxonomy:
        raise ValueError(f"Species missing from taxonomy table: {', '.join(missing_from_taxonomy)}")

    missing_from_tree = sorted(species_set - tree_names)
    if missing_from_tree:
        raise ValueError(f"Species missing from tree: {', '.join(missing_from_tree)}")

    return [species for species in tree_order if species in species_set]


def load_colors(path):
    """Read a flat YAML mapping clade names to colors."""
    if not path:
        return {}

    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    if not isinstance(data, dict):
        raise ValueError("The colors YAML file must contain a flat key-value mapping.")

    return {str(key): str(value) for key, value in data.items()}


def build_clade_species(taxonomy, clade_columns, species_order):
    """Map each clade name to its member species in tree order."""
    species_set = set(species_order)
    clade_members = defaultdict(set)

    for _, row in taxonomy.iterrows():
        species = row["name"]
        if species not in species_set:
            continue

        for column in clade_columns:
            clade = str(row.get(column, "")).strip()
            if clade and clade.lower() not in {"no rank", "clade", "nan"}:
                clade_members[clade].add(species)

    return {
        clade: [species for species in species_order if species in members]
        for clade, members in clade_members.items()
    }


def build_species_lineage(taxonomy, ranks, extra_clades, species_order, max_outgroup_rank, clade_species):
    """Build the ordered clades tested for each species."""
    if max_outgroup_rank not in ranks:
        raise ValueError(f"--max-outgroup-rank '{max_outgroup_rank}' is not present in --ranks")

    missing_clades = [clade for clade in extra_clades if clade not in clade_species]
    if missing_clades:
        raise ValueError(f"Requested clades were not found among taxonomy values: {', '.join(missing_clades)}")

    taxonomy_by_species = taxonomy.set_index("name", drop=False)
    max_outgroup_index = ranks.index(max_outgroup_rank)
    lineage = {}

    for species in species_order:
        row = taxonomy_by_species.loc[species]
        ranks_to_use = ranks[:max_outgroup_index + 1] if row["group"] == "outgroup" else ranks
        clades = []

        for rank in ranks_to_use:
            clade = str(row.get(rank, "")).strip()
            if not clade or clade.lower() in {"no rank", "clade", "nan"}:
                continue
            if len(clade_species.get(clade, [])) <= 1:
                continue

            current_members = set(clade_species[clade])
            if any(current_members == set(clade_species[previous]) for previous in clades):
                continue

            clades.append(clade)

        for clade in extra_clades:
            if species not in clade_species[clade]:
                continue
            if len(clade_species[clade]) <= 1:
                continue

            current_members = set(clade_species[clade])
            if any(current_members == set(clade_species[previous]) for previous in clades):
                continue

            clades.append(clade)

        lineage[species] = clades

    return lineage

def precompute_orthogroup_percentages(gene_counts, species_order, clade_species):
    """Precompute presence percentages for every orthogroup and clade."""
    percentages = {}
    present_species_by_orthogroup = {}

    for _, row in gene_counts.iterrows():
        orthogroup = row["Orthogroup"]
        present_species = [species for species in species_order if row[species] > 0]
        present_set = set(present_species)
        present_species_by_orthogroup[orthogroup] = present_species
        percentages[orthogroup] = {}

        for clade, members in clade_species.items():
            percentages[orthogroup][clade] = 100.0 * sum(species in present_set for species in members) / len(members)

    return percentages, present_species_by_orthogroup


def classify_gene(species, orthogroup, species_lineage, percentages, present_species_by_orthogroup, threshold, unassigned=False):
    """Assign one gene to a conservation category."""
    if unassigned:
        return {"species_specific": False, "unassigned": True, "final_rank": "Species-specific", "core_shared": ""}

    present_species = present_species_by_orthogroup.get(orthogroup, [])
    if present_species == [species]:
        return {"species_specific": True, "unassigned": False, "final_rank": "Species-specific", "core_shared": ""}

    for clade in species_lineage.get(species, []):
        percentage = percentages.get(orthogroup, {}).get(clade, 0.0)
        if percentage >= threshold:
            return {
                "species_specific": False,
                "unassigned": False,
                "final_rank": clade,
                "core_shared": "core" if percentage == 100.0 else "shared",
            }

    return {"species_specific": False, "unassigned": False, "final_rank": "Others", "core_shared": ""}


def init_species_writers(outdir, species_lineage, species_order):
    """Create one gene-level output table per species."""
    species_dir = outdir / "per_species_gene_tables"
    species_dir.mkdir(parents=True, exist_ok=True)
    handles = {}
    writers = {}

    for species in species_order:
        clades = species_lineage[species]
        handle = (species_dir / f"{safe_name(species)}.tsv").open("w", newline="", encoding="utf-8")
        fields = ["species", "gene", "orthogroup", *clades, "species_specific", "unassigned", "final_rank", "core_shared"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        handles[species] = handle
        writers[species] = writer

    return handles, writers


def close_handles(handles):
    """Close all per-species output files."""
    for handle in handles.values():
        handle.close()


def make_gene_row(species, gene, orthogroup, species_lineage, percentages, classification):
    """Create one gene-level output row."""
    row = {
        "species": species,
        "gene": gene,
        "orthogroup": orthogroup,
        "species_specific": str(classification["species_specific"]),
        "unassigned": str(classification["unassigned"]),
        "final_rank": classification["final_rank"],
        "core_shared": classification["core_shared"],
    }

    for clade in species_lineage.get(species, []):
        percentage = percentages.get(orthogroup, {}).get(clade, "")
        row[clade] = "" if percentage == "" else round(percentage, 2)

    return row


def category_from_classification(classification):
    """Convert a gene classification into a summary category."""
    if classification["final_rank"] in {"Species-specific", "Others"}:
        return classification["final_rank"]
    return f"{classification['final_rank']}_{classification['core_shared']}"


def process_orthogroups(orthogroups, species_order, species_lineage, percentages, present_species_by_orthogroup, threshold, writers):
    """Classify genes from assigned orthogroups."""
    summary = defaultdict(Counter)

    for _, row in orthogroups.iterrows():
        orthogroup = row["Orthogroup"]

        for species in species_order:
            genes = split_genes(row.get(species, ""))
            if not genes:
                continue

            classification = classify_gene(
                species,
                orthogroup,
                species_lineage,
                percentages,
                present_species_by_orthogroup,
                threshold,
            )
            category = category_from_classification(classification)

            for gene in genes:
                writers[species].writerow(
                    make_gene_row(species, gene, orthogroup, species_lineage, percentages, classification)
                )
                summary[species][category] += 1

    return summary


def process_unassigned(unassigned, species_order, species_lineage, writers):
    """Classify unassigned genes as species-specific."""
    summary = defaultdict(Counter)

    for _, row in unassigned.iterrows():
        orthogroup = row["Orthogroup"]

        for species in species_order:
            genes = split_genes(row.get(species, ""))
            if not genes:
                continue

            classification = classify_gene(species, orthogroup, species_lineage, {}, {}, 0, unassigned=True)

            for gene in genes:
                writers[species].writerow(
                    make_gene_row(species, gene, orthogroup, species_lineage, {}, classification)
                )
                summary[species]["Species-specific"] += 1

    return summary


def merge_summaries(first, second):
    """Merge two species/category count dictionaries."""
    for species, counts in second.items():
        first[species].update(counts)
    return first


def build_summary_table(summary, species_order, species_lineage):
    """Build the final species-by-category summary table in tree order."""
    categories = []

    for species in species_order:
        for clade in species_lineage[species]:
            for status in ("core", "shared"):
                category = f"{clade}_{status}"
                if category not in categories:
                    categories.append(category)

    categories.extend(["Others", "Species-specific"])
    rows = []

    for species in species_order:
        row = {"species": species}
        row.update({category: summary[species].get(category, 0) for category in categories})
        row["total"] = sum(row[category] for category in categories)
        rows.append(row)

    return pd.DataFrame(rows)


def parse_category(category):
    """Split a summary category into clade and core/shared status."""
    if category.endswith("_core"):
        return category[:-5], "core"
    if category.endswith("_shared"):
        return category[:-7], "shared"
    return category, ""


def get_color_for_clade(clade, colors, auto_colors, colormap):
    """Return the configured or automatically assigned color for a clade."""
    if clade in colors:
        return colors[clade]
    if clade not in auto_colors:
        auto_colors[clade] = colormap(len(auto_colors) % 20)
    return auto_colors[clade]


def plot_stacked(summary, colors, output_prefix, threshold):
    """Create horizontal stacked barplots in tree order."""
    categories = [column for column in summary.columns if column not in {"species", "total"}]
    figure_height = max(7, 0.35 * len(summary))
    figure, axis = plt.subplots(figsize=(16, figure_height))
    positions = range(len(summary))
    left = [0] * len(summary)
    colormap = plt.get_cmap("tab20")
    auto_colors = {}
    clades_in_plot = []

    for category in categories:
        values = summary[category].values
        clade, status = parse_category(category)

        if values.sum() == 0:
            continue

        if clade not in clades_in_plot:
            clades_in_plot.append(clade)

        color = get_color_for_clade(clade, colors, auto_colors, colormap)
        hatch = "///" if status == "core" else None
        edgecolor = "white" if status == "core" else "none"

        axis.barh(positions, values, left=left, color=color, edgecolor=edgecolor, linewidth=0.0, hatch=hatch)
        left = [current + value for current, value in zip(left, values)]

    axis.set_yticks(list(positions))
    axis.set_yticklabels(summary["species"], fontsize=9, fontstyle="italic")
    axis.invert_yaxis()
    axis.set_xlabel("Number of genes/proteins", fontsize=12, fontweight="bold")
    axis.set_title("Proteome composition by conservation level", fontsize=15, fontweight="bold")
    axis.xaxis.set_major_formatter(FuncFormatter(thousands))
    axis.grid(axis="x", linestyle="--", alpha=0.3)

    for spine in axis.spines.values():
        spine.set_visible(False)

    legend_handles = []

    for clade in colors:
        if clade in clades_in_plot:
            legend_handles.append(Patch(facecolor=colors[clade], edgecolor="none", label=clade))

    for clade in clades_in_plot:
        if clade not in colors:
            legend_handles.append(
                Patch(facecolor=get_color_for_clade(clade, colors, auto_colors, colormap), edgecolor="none", label=clade)
            )

    legend_handles.extend(
        [
            Patch(facecolor="black", edgecolor="black", label=f"Shared (≥{threshold:g}%)"),
            Patch(facecolor="black", edgecolor="white", hatch="///", linewidth=0.8, label="Core (100%)"),
        ]
    )

    axis.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)
    figure.tight_layout()
    figure.savefig(f"{output_prefix}.pdf")
    figure.savefig(f"{output_prefix}.svg")
    figure.savefig(f"{output_prefix}.png", dpi=300)
    plt.close(figure)


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    ranks = [rank.strip() for rank in args.ranks.split(",") if rank.strip()]
    extra_clades = [clade.strip() for clade in args.clades.split(",") if clade.strip()]

    gene_counts, species_columns = read_gene_counts(args.gene_counts)
    orthogroups = read_orthogroups(args.orthogroups)
    unassigned = read_orthogroups(args.unassigned)
    taxonomy = read_taxonomy(args.taxonomy, ranks)
    tree_order = read_tree_order(args.tree)
    species_order = validate_species(tree_order, taxonomy, species_columns)

    taxonomy = taxonomy.set_index("name").loc[species_order].reset_index()
    taxonomy.to_csv(args.outdir / "taxonomy.tsv", sep="\t", index=False)

    clade_columns = ranks
    clade_species = build_clade_species(taxonomy, clade_columns, species_order)
    species_lineage = build_species_lineage(
        taxonomy,
        ranks,
        extra_clades,
        species_order,
        args.max_outgroup_rank,
        clade_species,
    )

    percentages, present_species_by_orthogroup = precompute_orthogroup_percentages(
        gene_counts,
        species_order,
        clade_species,
    )

    handles, writers = init_species_writers(args.outdir, species_lineage, species_order)

    try:
        summary = process_orthogroups(
            orthogroups,
            species_order,
            species_lineage,
            percentages,
            present_species_by_orthogroup,
            args.threshold,
            writers,
        )
        summary = merge_summaries(
            summary,
            process_unassigned(unassigned, species_order, species_lineage, writers),
        )
    finally:
        close_handles(handles)

    summary_table = build_summary_table(summary, species_order, species_lineage)
    summary_path = args.outdir / "proteome_composition_summary.tsv"
    summary_table.to_csv(summary_path, sep="\t", index=False)

    colors = load_colors(args.colors)
    plot_stacked(
        summary_table,
        colors,
        args.outdir / "proteome_composition_stacked_barplot",
        args.threshold,
    )

    print(f"Done. Results written in: {args.outdir}")
    print(f"Summary table: {summary_path}")


if __name__ == "__main__":
    main()