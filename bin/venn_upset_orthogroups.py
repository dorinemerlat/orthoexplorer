#!/usr/bin/env python3

import argparse
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import yaml


def read_table(path):
    """Read a tab-separated table as strings and replace missing values with empty strings."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def safe_name(value):
    """Return a filesystem-safe name."""
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = value.lower()
    return value.strip("_")


def split_genes(cell):
    """Split an OrthoFinder cell containing comma-separated gene IDs."""
    if not isinstance(cell, str) or not cell.strip():
        return []
    return [gene.strip() for gene in cell.split(",") if gene.strip()]


def read_taxonomy(path, clades):
    """
    Read the taxonomy table and return clade -> species mappings.

    Each requested clade is searched across all taxonomy columns except the
    species-name column. A species cannot belong to more than one requested clade.
    """
    taxonomy = read_table(path)
    taxonomy.columns = [column.strip() for column in taxonomy.columns]

    if "name" not in taxonomy.columns:
        raise ValueError("The taxonomy table must contain a 'name' column.")

    taxonomy["name"] = taxonomy["name"].astype(str).str.strip()
    taxonomy_columns = [column for column in taxonomy.columns if column != "name"]
    clade_to_species = {}

    for clade in clades:
        mask = taxonomy[taxonomy_columns].astype(str).eq(clade).any(axis=1)
        species = taxonomy.loc[mask, "name"].tolist()

        if not species:
            raise ValueError(f"No species found for clade: {clade}")

        clade_to_species[clade] = species

    species_to_clades = {}
    for clade, species_list in clade_to_species.items():
        for species in species_list:
            species_to_clades.setdefault(species, []).append(clade)

    duplicated = {
        species: assigned_clades
        for species, assigned_clades in species_to_clades.items()
        if len(assigned_clades) > 1
    }

    if duplicated:
        message = "\n".join(
            f"{species}: {', '.join(assigned_clades)}"
            for species, assigned_clades in duplicated.items()
        )
        raise ValueError(
            "Some species belong to more than one requested clade:\n" + message
        )

    return clade_to_species


def read_gene_counts(path):
    """Read Orthogroups.GeneCount.tsv and convert species columns to integers."""
    gene_counts = read_table(path)
    gene_counts = gene_counts.rename(columns={gene_counts.columns[0]: "Orthogroup"})
    gene_counts = gene_counts.drop(columns=["Total"], errors="ignore")

    for column in gene_counts.columns:
        if column != "Orthogroup":
            gene_counts[column] = pd.to_numeric(
                gene_counts[column], errors="coerce"
            ).fillna(0).astype(int)

    return gene_counts


def read_orthogroups(path):
    """Read Orthogroups.tsv as a gene table."""
    orthogroups = read_table(path)
    orthogroups = orthogroups.rename(columns={orthogroups.columns[0]: "Orthogroup"})
    return orthogroups.drop(columns=["Total"], errors="ignore")


def validate_species_columns(gene_counts, orthogroups, clade_to_species):
    """Check that all selected species are present in both OrthoFinder tables."""
    selected_species = {
        species
        for species_list in clade_to_species.values()
        for species in species_list
    }

    missing_gene_counts = sorted(selected_species - set(gene_counts.columns))
    missing_orthogroups = sorted(selected_species - set(orthogroups.columns))

    messages = []
    if missing_gene_counts:
        messages.append(
            "Species missing from the gene-count table:\n"
            + "\n".join(missing_gene_counts)
        )
    if missing_orthogroups:
        messages.append(
            "Species missing from the orthogroups table:\n"
            + "\n".join(missing_orthogroups)
        )

    if messages:
        raise ValueError("\n\n".join(messages))


def build_presence_absence_table(gene_counts, clade_to_species):
    """
    Build an orthogroup x clade presence/absence table.

    A clade is present when an orthogroup occurs in at least one species from
    that clade. Paralogs in one species count as one species-level presence.
    """
    rows = []

    for _, row in gene_counts.iterrows():
        output_row = {"Orthogroup": row["Orthogroup"]}
        present_clades = []
        total_present_species = 0

        for clade, species_list in clade_to_species.items():
            n_present = sum(row[species] > 0 for species in species_list)
            total_present_species += n_present
            is_present = int(n_present >= 1)

            output_row[clade] = is_present
            output_row[f"{clade}_n_species"] = n_present

            if is_present:
                present_clades.append(clade)

        output_row["n_present_species_in_selected_clades"] = total_present_species
        output_row["Intersection"] = "&".join(present_clades) if present_clades else "None"
        rows.append(output_row)

    return pd.DataFrame(rows)


def remove_single_species_orthogroups(presence_absence):
    """Remove orthogroups present in only one species among the selected clades."""
    return (
        presence_absence[
            presence_absence["n_present_species_in_selected_clades"] > 1
        ]
        .copy()
        .reset_index(drop=True)
    )


def count_intersections(presence_absence):
    """Count orthogroups per intersection."""
    return (
        presence_absence["Intersection"]
        .value_counts()
        .rename_axis("Intersection")
        .reset_index(name="Orthogroups")
    )


def remove_none_intersection(intersection_counts):
    """Remove orthogroups not assigned to any requested clade."""
    return (
        intersection_counts[intersection_counts["Intersection"] != "None"]
        .copy()
        .reset_index(drop=True)
    )


def load_colors(path):
    """Read a flat colors.yaml file."""
    if not path:
        return {}

    with open(path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if not data:
        return {}

    return {str(key): str(value) for key, value in data.items()}


def sort_intersections_logically(intersection_counts, clades):
    """Sort intersections by intersection size and requested clade order."""
    clade_order = {clade: index for index, clade in enumerate(clades)}

    def sort_key(intersection):
        members = intersection.split("&")
        return len(members), [clade_order[member] for member in members]

    ordered = intersection_counts.copy()
    ordered["sort_key"] = ordered["Intersection"].apply(sort_key)

    return (
        ordered.sort_values("sort_key")
        .drop(columns="sort_key")
        .reset_index(drop=True)
    )


def write_intersection_gene_files(
    orthogroups, presence_absence, clade_to_species, outdir
):
    """
    Write gene lists for each intersection.

    Output columns are intersection, orthogroup, species and gene. Paralogs are
    written on separate lines.
    """
    output_dir = os.path.join(outdir, "intersection_gene_lists")
    os.makedirs(output_dir, exist_ok=True)

    selected_species = []
    for species_list in clade_to_species.values():
        selected_species.extend(species_list)
    selected_species = list(dict.fromkeys(selected_species))

    intersection_by_orthogroup = dict(
        zip(presence_absence["Orthogroup"], presence_absence["Intersection"])
    )

    all_rows = []
    rows_by_intersection = {}

    for _, row in orthogroups.iterrows():
        orthogroup = row["Orthogroup"]
        intersection = intersection_by_orthogroup.get(orthogroup, "None")

        if intersection == "None":
            continue

        for species in selected_species:
            for gene in split_genes(row[species]):
                output_row = {
                    "intersection": intersection,
                    "orthogroup": orthogroup,
                    "species": species,
                    "gene": gene,
                }
                all_rows.append(output_row)
                rows_by_intersection.setdefault(intersection, []).append(output_row)

    columns = ["intersection", "orthogroup", "species", "gene"]
    pd.DataFrame(all_rows, columns=columns).to_csv(
        os.path.join(outdir, "intersection_genes_all.tsv"),
        sep="\t",
        index=False,
    )

    for intersection, rows in rows_by_intersection.items():
        filename = safe_name(intersection) + ".tsv"
        pd.DataFrame(rows, columns=columns).to_csv(
            os.path.join(output_dir, filename),
            sep="\t",
            index=False,
        )


def write_eligible_orthogroups(presence_absence, outdir):
    """Write the orthogroups defining the UpSet enrichment universe."""
    eligible_orthogroups = (
        presence_absence.loc[
            presence_absence["Intersection"] != "None",
            ["Orthogroup"],
        ]
        .drop_duplicates()
        .sort_values("Orthogroup")
        .reset_index(drop=True)
    )

    eligible_orthogroups.to_csv(
        os.path.join(outdir, "eligible_orthogroups.txt"),
        index=False,
        header=False,
    )

    return len(eligible_orthogroups)


def plot_venn(presence_absence, clades, colors, outdir):
    """Plot a Venn diagram from the orthogroup presence/absence matrix."""
    try:
        from venn import venn
    except ImportError as error:
        raise ImportError(
            "The 'venn' package is required. Install it with: mamba install -c conda-forge venn"
        ) from error

    datasets = {
        clade: set(
            presence_absence.loc[presence_absence[clade] == 1, "Orthogroup"]
        )
        for clade in clades
    }

    fig, ax = plt.subplots(figsize=(10, 8))
    venn_colors = [colors[clade] for clade in clades if clade in colors]
    cmap = venn_colors if len(venn_colors) == len(clades) else "Set3"

    venn(
        datasets,
        fmt="{size}",
        cmap=cmap,
        fontsize=8,
        legend_loc="upper left",
        ax=ax,
    )

    for extension in ["png", "pdf", "svg"]:
        fig.savefig(
            os.path.join(outdir, f"venn_orthogroups.{extension}"),
            dpi=300,
            bbox_inches="tight",
        )

    plt.close(fig)


def plot_upset(presence_absence, clades, outdir):
    """Plot an UpSet plot ordered by intersection size and clade order."""
    try:
        from upsetplot import UpSet, from_memberships
    except ImportError as error:
        raise ImportError(
            "The 'upsetplot' package is required. Install it with: mamba install -c conda-forge upsetplot"
        ) from error

    intersection_counts = (
        presence_absence.groupby("Intersection", sort=False)
        .size()
        .reset_index(name="count")
    )
    intersection_counts = remove_none_intersection(intersection_counts)
    intersection_counts = sort_intersections_logically(intersection_counts, clades)

    memberships = []
    counts = []

    for _, row in intersection_counts.iterrows():
        members = row["Intersection"].split("&")
        memberships.append([clade for clade in clades if clade in members])
        counts.append(row["count"])

    upset_data = from_memberships(memberships, data=counts)
    upset_data = upset_data.reorder_levels(list(reversed(clades)))

    fig = plt.figure(figsize=(12, 7))
    upset = UpSet(
        upset_data,
        subset_size="sum",
        show_counts=True,
        sort_by="input",
        sort_categories_by=None,
    )
    upset.plot(fig=fig)

    for extension in ["png", "pdf", "svg"]:
        fig.savefig(
            os.path.join(outdir, f"upset_orthogroups.{extension}"),
            dpi=300,
            bbox_inches="tight",
        )

    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Create Venn and UpSet plots from OrthoFinder orthogroups."
    )
    parser.add_argument(
        "--gene-counts",
        required=True,
        help="Reformatted Orthogroups.GeneCount.tsv file.",
    )
    parser.add_argument(
        "--orthogroups",
        required=True,
        help="Reformatted Orthogroups.tsv file.",
    )
    parser.add_argument(
        "--taxonomy",
        required=True,
        help="Taxonomy table containing a 'name' column.",
    )
    parser.add_argument(
        "--clades",
        required=True,
        help="Comma-separated clade names.",
    )
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument(
        "--colors",
        default=None,
        help="Optional flat YAML file mapping clade names to colors.",
    )
    parser.add_argument(
        "--keep-single-species",
        action="store_true",
        help="Keep orthogroups present in only one selected species.",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    clades = [clade.strip() for clade in args.clades.split(",") if clade.strip()]
    if len(clades) < 2:
        raise ValueError("--clades must contain at least two comma-separated clades.")

    clade_to_species = read_taxonomy(args.taxonomy, clades)
    gene_counts = read_gene_counts(args.gene_counts)
    orthogroups = read_orthogroups(args.orthogroups)
    validate_species_columns(gene_counts, orthogroups, clade_to_species)

    presence_absence = build_presence_absence_table(
        gene_counts=gene_counts,
        clade_to_species=clade_to_species,
    )

    if not args.keep_single_species:
        presence_absence = remove_single_species_orthogroups(presence_absence)

    presence_absence.to_csv(
        os.path.join(outdir := args.outdir, "orthogroup_presence_absence_by_clade.tsv"),
        sep="\t",
        index=False,
    )

    presence_absence.to_csv(
        os.path.join(outdir := args.outdir, "orthogroup_presence_absence_by_clade.tsv"),
        sep="\t",
        index=False,
    )

    eligible_orthogroup_count = write_eligible_orthogroups(
        presence_absence=presence_absence,
        outdir=outdir,
    )

    intersection_counts = count_intersections(presence_absence)
    intersection_counts = remove_none_intersection(intersection_counts)
    intersection_counts = sort_intersections_logically(
        intersection_counts,
        clades,
    )
    intersection_counts.to_csv(
        os.path.join(outdir, "intersection_counts.tsv"),
        sep="\t",
        index=False,
    )

    write_intersection_gene_files(
        orthogroups=orthogroups,
        presence_absence=presence_absence,
        clade_to_species=clade_to_species,
        outdir=outdir,
    )

    colors = load_colors(args.colors)
    plot_upset(presence_absence, clades, outdir)
    plot_venn(presence_absence, clades, colors, outdir)

    print(f"Done. Results written in: {outdir}")


if __name__ == "__main__":
    main()