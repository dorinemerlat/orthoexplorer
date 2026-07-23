#!/usr/bin/env python3

import argparse
import csv
import os
import re
from collections import Counter, defaultdict

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

try:
    import yaml
except ImportError:
    yaml = None


DEFAULT_RANKS = ["phylum", "subphylum", "class", "subclass", "order"]


def safe_name(value):
    """Return a filesystem-safe name."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())


def split_genes(cell):
    """Split an OrthoFinder gene list cell into individual gene IDs."""
    if not isinstance(cell, str) or not cell.strip():
        return []
    return [gene.strip() for gene in cell.split(",") if gene.strip()]


def thousands(x, pos):
    """Format axis labels with spaces as thousands separators."""
    return f"{int(x):,}".replace(",", " ")


def read_table(path):
    """Read a tab-separated table as strings and replace missing values by empty strings."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def read_gene_counts(path):
    """Read Orthogroups.GeneCount.tsv and convert species columns to integers."""
    df = read_table(path)
    df = df.rename(columns={df.columns[0]: "Orthogroup"})
    df = df.drop(columns=["Total"], errors="ignore")

    species_cols = [col for col in df.columns if col != "Orthogroup"]

    for species in species_cols:
        df[species] = pd.to_numeric(df[species], errors="coerce").fillna(0).astype(int)

    return df, species_cols


def read_orthogroups(path):
    """Read Orthogroups.tsv or Orthogroups_UnassignedGenes.tsv."""
    df = read_table(path)
    df = df.rename(columns={df.columns[0]: "Orthogroup"})
    df = df.drop(columns=["Total"], errors="ignore")
    return df


def read_taxonomy(path, ranks):
    """Read taxonomy table and ensure requested rank columns exist."""
    taxonomy = read_table(path)
    taxonomy.columns = [col.strip() for col in taxonomy.columns]

    if "specie" not in taxonomy.columns:
        raise ValueError("The taxonomy table must contain a 'specie' column.")

    if "taxid" not in taxonomy.columns:
        taxonomy["taxid"] = ""

    for rank in ranks:
        if rank not in taxonomy.columns:
            taxonomy[rank] = ""

    return taxonomy


def read_outgroups(path):
    """Read one outgroup taxid per line."""
    if not path:
        return set()

    with open(path, encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def load_colors(path):
    """
    Read a flat colors.yaml file.

    Example:
    Arthropoda: "#a6cee3"
    Myriapoda: "#1f78b4"
    Others: "#808080"
    """
    if not path:
        return {}

    if yaml is not None:
        with open(path, encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
        if not data:
            return {}
        return {str(key): str(value) for key, value in data.items()}

    # Fallback if PyYAML is not installed.
    # This expects simple lines such as: Clade: "#a6cee3"
    colors = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue

            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()

            if " #" in value:
                value = value.split(" #", 1)[0].strip()

            value = value.strip('"').strip("'")

            if key and value:
                colors[key] = value

    return colors


def build_clade_species(taxonomy, ranks, species_cols):
    """
    Build a dictionary mapping each clade name to the list of species it contains.

    Example:
    {
        "Arthropoda": ["A", "B", "C"],
        "Myriapoda": ["A", "B"],
        "Diplopoda": ["A"]
    }
    """
    species_set = set(species_cols)
    clade_species = defaultdict(list)

    for _, row in taxonomy.iterrows():
        species = row["specie"]

        if species not in species_set:
            continue

        for rank in ranks:
            clade = str(row.get(rank, "")).strip()

            if not clade:
                continue

            if clade.lower() in {"no rank", "clade", "nan"}:
                continue

            if species not in clade_species[clade]:
                clade_species[clade].append(species)

    return clade_species


def build_species_lineage(
    taxonomy,
    ranks,
    species_cols,
    outgroup_taxids,
    max_outgroup_rank,
    clade_species,
):
    """
    Build the ordered list of clades to test for each species.

    Rules:
    - ranks are tested from broad to specific;
    - outgroup species can be stopped at --max-outgroup-rank;
    - clades represented by only one species are removed;
    - if two ranks contain exactly the same species, keep only the broadest one.
    """
    lineage = {}
    tax_by_species = taxonomy.set_index("specie", drop=False)

    for species in species_cols:
        if species not in tax_by_species.index:
            lineage[species] = []
            continue

        row = tax_by_species.loc[species]
        taxid = str(row.get("taxid", "")).strip()

        if taxid in outgroup_taxids and max_outgroup_rank in ranks:
            max_index = ranks.index(max_outgroup_rank)
            ranks_to_use = ranks[:max_index + 1]
        else:
            ranks_to_use = ranks

        clades = []

        for rank in ranks_to_use:
            clade = str(row.get(rank, "")).strip()

            if not clade:
                continue

            if clade.lower() in {"no rank", "clade", "nan"}:
                continue

            if len(clade_species.get(clade, [])) <= 1:
                continue

            current_species_set = set(clade_species.get(clade, []))

            duplicate_of_previous = False
            for previous_clade in clades:
                previous_species_set = set(clade_species.get(previous_clade, []))

                if current_species_set == previous_species_set:
                    duplicate_of_previous = True
                    break

            if duplicate_of_previous:
                continue

            clades.append(clade)

        lineage[species] = clades

    return lineage


def precompute_orthogroup_percentages(gene_counts, species_cols, clade_species):
    """
    Precompute, for each orthogroup, the percentage of species in each clade
    where the orthogroup is present.

    This avoids recalculating the same value for each paralogous gene.
    """
    percentages = {}
    present_species_by_og = {}

    for _, row in gene_counts.iterrows():
        orthogroup = row["Orthogroup"]

        present_species = [species for species in species_cols if row[species] > 0]
        present_set = set(present_species)

        present_species_by_og[orthogroup] = present_species
        percentages[orthogroup] = {}

        for clade, species_list in clade_species.items():
            if not species_list:
                percentages[orthogroup][clade] = 0.0
                continue

            n_present = sum(species in present_set for species in species_list)
            percentages[orthogroup][clade] = 100.0 * n_present / len(species_list)

    return percentages, present_species_by_og


def classify_gene(
    species,
    orthogroup,
    species_lineage,
    og_percentages,
    present_species_by_og,
    threshold,
    unassigned=False,
):
    """
    Assign one gene to a conservation category.

    Priority:
    1. unassigned genes -> Species-specific
    2. orthogroups present only in this species -> Species-specific
    3. first clade in the species lineage with conservation >= threshold
    4. otherwise -> Others
    """
    if unassigned:
        return {
            "species_specific": False,
            "unassigned": True,
            "final_rank": "Species-specific",
            "core_shared": "",
        }

    present_species = present_species_by_og.get(orthogroup, [])
    species_specific = len(present_species) == 1 and present_species[0] == species

    if species_specific:
        return {
            "species_specific": True,
            "unassigned": False,
            "final_rank": "Species-specific",
            "core_shared": "",
        }

    for clade in species_lineage.get(species, []):
        percentage = og_percentages.get(orthogroup, {}).get(clade, 0.0)

        if percentage >= threshold:
            return {
                "species_specific": False,
                "unassigned": False,
                "final_rank": clade,
                "core_shared": "core" if percentage == 100.0 else "shared",
            }

    return {
        "species_specific": False,
        "unassigned": False,
        "final_rank": "Others",
        "core_shared": "",
    }


def init_species_writers(outdir, species_lineage):
    """Create one gene-level output table per species."""
    species_dir = os.path.join(outdir, "per_species_gene_tables")
    os.makedirs(species_dir, exist_ok=True)

    handles = {}
    writers = {}

    for species, clades in species_lineage.items():
        path = os.path.join(species_dir, f"{safe_name(species)}.tsv")
        handle = open(path, "w", newline="", encoding="utf-8")

        header = (
            ["specie", "gene", "orthogroup"]
            + clades
            + ["species_specific", "unassigned", "final_rank", "core_shared"]
        )

        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()

        handles[species] = handle
        writers[species] = writer

    return handles, writers


def close_handles(handles):
    """Close all open per-species output files."""
    for handle in handles.values():
        handle.close()


def make_gene_row(species, gene, orthogroup, species_lineage, og_percentages, classification):
    """Create one row for the per-species gene-level table."""
    row = {
        "specie": species,
        "gene": gene,
        "orthogroup": orthogroup,
        "species_specific": str(classification["species_specific"]),
        "unassigned": str(classification["unassigned"]),
        "final_rank": classification["final_rank"],
        "core_shared": classification["core_shared"],
    }

    for clade in species_lineage.get(species, []):
        percentage = og_percentages.get(orthogroup, {}).get(clade, "")
        row[clade] = "" if percentage == "" else round(percentage, 2)

    return row


def category_from_classification(classification):
    """Convert a gene-level classification into a stacked-bar category."""
    final_rank = classification["final_rank"]
    core_shared = classification["core_shared"]

    if final_rank == "Species-specific":
        return "Species-specific"

    if final_rank == "Others":
        return "Others"

    return f"{final_rank}_{core_shared}"


def process_orthogroups(
    orthogroups,
    species_cols,
    species_lineage,
    og_percentages,
    present_species_by_og,
    threshold,
    writers,
):
    """Process assigned orthogroups and write gene-level rows."""
    summary = defaultdict(Counter)

    for _, row in orthogroups.iterrows():
        orthogroup = row["Orthogroup"]

        for species in species_cols:
            if species not in row.index:
                continue

            genes = split_genes(row[species])

            if not genes:
                continue

            classification = classify_gene(
                species=species,
                orthogroup=orthogroup,
                species_lineage=species_lineage,
                og_percentages=og_percentages,
                present_species_by_og=present_species_by_og,
                threshold=threshold,
                unassigned=False,
            )

            category = category_from_classification(classification)

            for gene in genes:
                out_row = make_gene_row(
                    species=species,
                    gene=gene,
                    orthogroup=orthogroup,
                    species_lineage=species_lineage,
                    og_percentages=og_percentages,
                    classification=classification,
                )

                writers[species].writerow(out_row)
                summary[species][category] += 1

    return summary


def process_unassigned(unassigned, species_cols, species_lineage, writers):
    """Process unassigned genes and count them as Species-specific."""
    summary = defaultdict(Counter)

    if unassigned is None:
        return summary

    for _, row in unassigned.iterrows():
        orthogroup = row["Orthogroup"]

        for species in species_cols:
            if species not in row.index:
                continue

            genes = split_genes(row[species])

            if not genes:
                continue

            classification = classify_gene(
                species=species,
                orthogroup=orthogroup,
                species_lineage=species_lineage,
                og_percentages={},
                present_species_by_og={},
                threshold=0,
                unassigned=True,
            )

            for gene in genes:
                out_row = make_gene_row(
                    species=species,
                    gene=gene,
                    orthogroup=orthogroup,
                    species_lineage=species_lineage,
                    og_percentages={},
                    classification=classification,
                )

                writers[species].writerow(out_row)
                summary[species]["Species-specific"] += 1

    return summary


def merge_summaries(first, second):
    """Merge two species -> category -> count dictionaries."""
    for species, counts in second.items():
        for category, n in counts.items():
            first[species][category] += n
    return first


def build_summary_table(summary, taxonomy, species_cols, species_lineage):
    """
    Build the final species x category table.

    Each row total should correspond to the number of proteins for the species.
    """
    categories = []

    for species in species_cols:
        for clade in species_lineage.get(species, []):
            for status in ["core", "shared"]:
                category = f"{clade}_{status}"
                if category not in categories:
                    categories.append(category)

    categories += ["Others", "Species-specific"]

    rows = []
    taxonomy_order = [species for species in taxonomy["specie"].tolist() if species in species_cols]

    for species in taxonomy_order:
        row = {"specie": species}
        total = 0

        for category in categories:
            n = summary[species].get(category, 0)
            row[category] = n
            total += n

        row["total"] = total
        rows.append(row)

    return pd.DataFrame(rows)


def parse_category(category):
    """Return the clade and core/shared status from a category name."""
    if category.endswith("_core"):
        return category.replace("_core", ""), "core"

    if category.endswith("_shared"):
        return category.replace("_shared", ""), "shared"

    return category, ""


def get_color_for_clade(clade, colors, auto_colors, cmap):
    """Return the color of a clade, prioritizing the user-provided colors.yaml."""
    if clade in colors:
        return colors[clade]

    if clade not in auto_colors:
        auto_colors[clade] = cmap(len(auto_colors) % 20)

    return auto_colors[clade]


def plot_stacked(summary_df, colors, output_prefix, threshold):
    """
    Make a horizontal stacked barplot.

    Plot style:
    - same color for core/shared categories of the same clade;
    - shared categories are shown with white hatches;
    - legend has one color entry per clade and a separate core/shared style legend;
    - no enclosing frame around the plot.
    """
    categories = [col for col in summary_df.columns if col not in {"specie", "total"}]

    fig_height = max(7, 0.35 * len(summary_df))
    fig, ax = plt.subplots(figsize=(16, fig_height))

    y = range(len(summary_df))
    left = [0] * len(summary_df)

    cmap = plt.get_cmap("tab20")
    auto_colors = {}
    clades_in_plot = []

    for category in categories:
        values = summary_df[category].values
        clade, status = parse_category(category)

        if values.sum() == 0:
            continue

        if clade not in clades_in_plot:
            clades_in_plot.append(clade)

        color = get_color_for_clade(clade, colors, auto_colors, cmap)

        # Shared = same clade color + white hatches.
        # Core and non-core special categories = solid fill.
        if category in {"Others", "Species-specific"}:
            hatch = None
            edgecolor = "none"

        elif status == "core":
            hatch = "///"
            edgecolor = "white"

        else:  # shared
            hatch = None
            edgecolor = "none"

        ax.barh(
            y,
            values,
            left=left,
            color=color,
            edgecolor=edgecolor,
            linewidth=0.0,
            hatch=hatch,
        )

        left = [old + value for old, value in zip(left, values)]

    ax.set_yticks(list(y))
    ax.set_yticklabels(summary_df["specie"], fontsize=9, fontstyle="italic")
    ax.invert_yaxis()

    ax.set_xlabel("Number of genes/proteins", fontsize=12, fontweight="bold")
    ax.set_title(
        "Proteome composition by conservation level",
        fontsize=15,
        fontweight="bold",
    )

    ax.xaxis.set_major_formatter(FuncFormatter(thousands))
    ax.grid(axis="x", linestyle="--", alpha=0.3)

    # Remove the enclosing frame around the plotting area.
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Build one combined legend:
    # - one colored patch per clade
    # - then one solid white patch for core
    # - and one hatched white patch for shared
    legend_handles = []

    # Clades dans l'ordre du YAML
    for clade in colors:
        if clade in clades_in_plot:
            legend_handles.append(
                Patch(
                    facecolor=colors[clade],
                    edgecolor="none",
                    label=clade,
                )
            )

    # Clades non définis dans le YAML (sécurité)
    for clade in clades_in_plot:
        if clade not in colors:
            legend_handles.append(
                Patch(
                    facecolor=get_color_for_clade(clade, colors, auto_colors, cmap),
                    edgecolor="none",
                    label=clade,
                )
            )
        
    legend_handles.extend([
        Patch(
            facecolor="black",
            edgecolor="black",
            label="Shared (≥75%)"
        ),
        Patch(
            facecolor="black",
            edgecolor="white",
            hatch="///",
            linewidth=0.8,
            label="Core (100%)"
        ),
    ])

    ax.legend(
        handles=legend_handles,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
        frameon=False,
    )

    plt.tight_layout()

    fig.savefig(output_prefix + ".pdf")
    fig.savefig(output_prefix + ".svg")
    fig.savefig(output_prefix + ".png", dpi=300)

    plt.close(fig)


def reorder_taxonomy_by_phylogeny(taxonomy):
    """
    Reorder taxonomy rows alphabetically from kingdom onward,
    using columns from left to right to reflect taxonomy hierarchy.
    """
    columns = list(taxonomy.columns)

    if "kingdom" not in columns:
        raise ValueError("Cannot reorder taxonomy: no 'kingdom' column found.")

    start = columns.index("kingdom")
    sort_columns = columns[start:]

    # Keep only informative taxonomy columns
    sort_columns = [
        col for col in sort_columns
        if col not in {"taxid", "specie"}
    ]

    if "specie" in taxonomy.columns:
        sort_columns.append("specie")

    return taxonomy.sort_values(
        by=sort_columns,
        key=lambda col: col.astype(str).str.lower(),
        na_position="last"
    ).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(
        description="Assign each gene to a conservation rank from OrthoFinder results."
    )

    parser.add_argument("--gene-counts", required=True, help="Orthogroups.GeneCount.tsv")
    parser.add_argument("--orthogroups", required=True, help="Orthogroups.tsv")
    parser.add_argument("--unassigned", required=True, help="Orthogroups_UnassignedGenes.tsv")
    parser.add_argument("--taxonomy", required=True, help="taxonomy_table.tsv")
    parser.add_argument("--outdir", required=True, help="Output directory")

    parser.add_argument(
        "--ranks",
        default=",".join(DEFAULT_RANKS),
        help="Comma-separated taxonomy ranks to use"
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=75.0,
        help="Minimum conservation percentage for final_rank assignment"
    )

    parser.add_argument(
        "--outgroups",
        default=None,
        help="File with one outgroup taxid per line"
    )

    parser.add_argument(
        "--max-outgroup-rank",
        default="subphylum",
        help="Maximum taxonomy rank used for outgroup species"
    )

    parser.add_argument(
        "--colors-yaml",
        default=None,
        help="Optional flat colors.yaml file"
    )
    
    parser.add_argument(
        "--reorder-table",
        action="store_true",
        help="Reorder taxonomy/species order alphabetically from kingdom onward"
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    ranks = [rank.strip() for rank in args.ranks.split(",") if rank.strip()]

    taxonomy = read_taxonomy(args.taxonomy, ranks)
    if args.reorder_table:
        taxonomy = reorder_taxonomy_by_phylogeny(taxonomy)

    taxonomy_path = os.path.join(args.outdir, "taxonomy.tsv")
    taxonomy.to_csv(taxonomy_path, sep="\t", index=False)

    gene_counts, species_cols = read_gene_counts(args.gene_counts)
    orthogroups = read_orthogroups(args.orthogroups)
    unassigned = read_orthogroups(args.unassigned)

    colors = load_colors(args.colors_yaml)
    outgroup_taxids = read_outgroups(args.outgroups)

    clade_species = build_clade_species(
        taxonomy=taxonomy,
        ranks=ranks,
        species_cols=species_cols,
    )

    species_lineage = build_species_lineage(
        taxonomy=taxonomy,
        ranks=ranks,
        species_cols=species_cols,
        outgroup_taxids=outgroup_taxids,
        max_outgroup_rank=args.max_outgroup_rank,
        clade_species=clade_species,
    )

    og_percentages, present_species_by_og = precompute_orthogroup_percentages(
        gene_counts=gene_counts,
        species_cols=species_cols,
        clade_species=clade_species,
    )

    handles, writers = init_species_writers(
        outdir=args.outdir,
        species_lineage=species_lineage,
    )

    try:
        assigned_summary = process_orthogroups(
            orthogroups=orthogroups,
            species_cols=species_cols,
            species_lineage=species_lineage,
            og_percentages=og_percentages,
            present_species_by_og=present_species_by_og,
            threshold=args.threshold,
            writers=writers,
        )

        unassigned_summary = process_unassigned(
            unassigned=unassigned,
            species_cols=species_cols,
            species_lineage=species_lineage,
            writers=writers,
        )

        summary = merge_summaries(assigned_summary, unassigned_summary)

    finally:
        close_handles(handles)

    summary_df = build_summary_table(
        summary=summary,
        taxonomy=taxonomy,
        species_cols=species_cols,
        species_lineage=species_lineage,
    )

    summary_path = os.path.join(args.outdir, "proteome_composition_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)

    plot_stacked(
        summary_df=summary_df,
        colors=colors,
        output_prefix=os.path.join(args.outdir, "proteome_composition_stacked_barplot"),
        threshold=args.threshold,
    )

    print(f"Done. Results written in: {args.outdir}")
    print(f"Summary table: {summary_path}")


if __name__ == "__main__":
    main()