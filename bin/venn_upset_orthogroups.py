#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd
import matplotlib.pyplot as plt
import yaml


def read_table(path):
    """Read a tab-separated table as strings and replace missing values by empty strings."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def safe_name(value):
    """Return a filesystem-safe name."""
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def split_genes(cell):
    """Split an OrthoFinder cell containing comma-separated gene IDs."""
    if not isinstance(cell, str) or not cell.strip():
        return []
    return [gene.strip() for gene in cell.split(",") if gene.strip()]


def read_taxonomy(path, groups):
    """
    Read taxonomy table and return:
    group_name -> list of species in this group.

    Also checks that no species belongs to more than one requested group.
    """
    taxonomy = read_table(path)
    taxonomy.columns = [col.strip() for col in taxonomy.columns]

    if "specie" not in taxonomy.columns:
        raise ValueError("The taxonomy table must contain a 'specie' column.")

    group_to_species = {}

    for group in groups:
        mask = (taxonomy.astype(str) == group).any(axis=1)
        species = taxonomy.loc[mask, "specie"].tolist()

        if not species:
            raise ValueError(f"No species found for group: {group}")

        group_to_species[group] = species

    species_to_groups = {}

    for group, species_list in group_to_species.items():
        for species in species_list:
            species_to_groups.setdefault(species, []).append(group)

    duplicated = {
        species: assigned_groups
        for species, assigned_groups in species_to_groups.items()
        if len(assigned_groups) > 1
    }

    if duplicated:
        message = "\n".join(
            f"{species}: {', '.join(assigned_groups)}"
            for species, assigned_groups in duplicated.items()
        )
        raise ValueError(
            "Some species belong to more than one requested group:\n"
            + message
        )

    return group_to_species


def read_gene_counts(path):
    """Read Orthogroups.GeneCount.tsv and convert species columns to integers."""
    gene_counts = read_table(path)
    gene_counts = gene_counts.rename(columns={gene_counts.columns[0]: "Orthogroup"})
    gene_counts = gene_counts.drop(columns=["Total"], errors="ignore")

    for col in gene_counts.columns:
        if col != "Orthogroup":
            gene_counts[col] = pd.to_numeric(
                gene_counts[col],
                errors="coerce"
            ).fillna(0).astype(int)

    return gene_counts


def read_orthogroups(path):
    """Read Orthogroups.tsv as a gene table."""
    orthogroups = read_table(path)
    orthogroups = orthogroups.rename(columns={orthogroups.columns[0]: "Orthogroup"})
    orthogroups = orthogroups.drop(columns=["Total"], errors="ignore")
    return orthogroups


def build_presence_absence_table(gene_counts, group_to_species):
    """
    Build an orthogroup x group presence/absence table.

    A group is coded as 1 if the orthogroup is present in at least
    one species of that group.

    The column n_present_species_in_selected_groups counts the number of
    species, not the number of genes. Paralogs in the same species count as
    one present species.
    """
    rows = []

    for _, row in gene_counts.iterrows():
        out_row = {"Orthogroup": row["Orthogroup"]}
        present_groups = []
        total_present_species = 0

        for group, species_list in group_to_species.items():
            n_present = 0

            for species in species_list:
                if species not in gene_counts.columns:
                    continue

                if row[species] > 0:
                    n_present += 1

            total_present_species += n_present

            is_present = int(n_present >= 1)
            out_row[group] = is_present
            out_row[f"{group}_n_species"] = n_present

            if is_present:
                present_groups.append(group)

        out_row["n_present_species_in_selected_groups"] = total_present_species
        out_row["Intersection"] = "&".join(present_groups) if present_groups else "None"

        rows.append(out_row)

    return pd.DataFrame(rows)


def remove_single_species_orthogroups(presence_absence):
    """Remove orthogroups present in only one species among selected groups."""
    return (
        presence_absence[
            presence_absence["n_present_species_in_selected_groups"] > 1
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
    """Remove orthogroups not assigned to any requested group."""
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


def sort_intersections_logically(intersection_counts, groups):
    """
    Sort intersections by:
    1. number of groups in the intersection;
    2. user-provided group order.
    """
    group_order = {group: i for i, group in enumerate(groups)}

    def sort_key(intersection):
        members = intersection.split("&")
        members_order = [group_order[group] for group in members]
        return (len(members), members_order)

    ordered = intersection_counts.copy()
    ordered["sort_key"] = ordered["Intersection"].apply(sort_key)

    return (
        ordered
        .sort_values("sort_key")
        .drop(columns="sort_key")
        .reset_index(drop=True)
    )


def write_intersection_gene_files(orthogroups, presence_absence, group_to_species, outdir):
    """
    Write gene lists for each UpSet intersection.

    Output columns:
    intersection, orthogroup, species, gene

    If several paralogs are present for the same species and orthogroup,
    one line is written per gene.
    """
    output_dir = os.path.join(outdir, "intersection_gene_lists")
    os.makedirs(output_dir, exist_ok=True)

    selected_species = []
    for species_list in group_to_species.values():
        selected_species.extend(species_list)

    selected_species = list(dict.fromkeys(selected_species))

    intersection_by_og = dict(
        zip(
            presence_absence["Orthogroup"],
            presence_absence["Intersection"],
        )
    )

    all_rows = []
    rows_by_intersection = {}

    for _, row in orthogroups.iterrows():
        orthogroup = row["Orthogroup"]
        intersection = intersection_by_og.get(orthogroup, "None")

        if intersection == "None":
            continue

        for species in selected_species:
            if species not in orthogroups.columns:
                continue

            genes = split_genes(row[species])

            for gene in genes:
                out_row = {
                    "intersection": intersection,
                    "orthogroup": orthogroup,
                    "species": species,
                    "gene": gene,
                }

                all_rows.append(out_row)
                rows_by_intersection.setdefault(intersection, []).append(out_row)

    all_genes = pd.DataFrame(
        all_rows,
        columns=["intersection", "orthogroup", "species", "gene"],
    )

    all_genes.to_csv(
        os.path.join(outdir, "intersection_genes_all.tsv"),
        sep="\t",
        index=False,
    )

    for intersection, rows in rows_by_intersection.items():
        filename = safe_name(intersection) + ".tsv"
        path = os.path.join(output_dir, filename)

        pd.DataFrame(
            rows,
            columns=["intersection", "orthogroup", "species", "gene"],
        ).to_csv(path, sep="\t", index=False)


def plot_venn(presence_absence, groups, colors, outdir):
    """Plot a Venn diagram from the orthogroup presence/absence matrix."""
    try:
        from venn import venn
    except ImportError as e:
        raise ImportError(
            "The 'venn' package is required. Install it with: pip install venn"
        ) from e

    dataset_dict = {}

    for group in groups:
        dataset_dict[group] = set(
            presence_absence.loc[
                presence_absence[group] == 1,
                "Orthogroup"
            ]
        )

    fig, ax = plt.subplots(figsize=(10, 8))

    venn_colors = []
    for group in groups:
        if group in colors:
            venn_colors.append(colors[group])

    cmap = venn_colors if len(venn_colors) == len(groups) else "Set3"

    venn(
        dataset_dict,
        fmt="{size}",
        cmap=cmap,
        fontsize=8,
        legend_loc="upper left",
        ax=ax,
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(
            os.path.join(outdir, f"venn_orthogroups.{ext}"),
            dpi=300,
            bbox_inches="tight",
        )

    plt.close(fig)


def plot_upset(presence_absence, groups, outdir):
    """
    Plot an UpSet plot from the orthogroup presence/absence matrix.

    Intersections are ordered logically:
    single groups -> pairs -> triplets -> quadruplets -> all groups.
    """
    try:
        from upsetplot import UpSet, from_memberships
    except ImportError as e:
        raise ImportError(
            "The 'upsetplot' package is required. Install it with: pip install upsetplot"
        ) from e

    intersection_counts = (
        presence_absence
        .groupby("Intersection", sort=False)
        .size()
        .reset_index(name="count")
    )

    intersection_counts = remove_none_intersection(intersection_counts)

    intersection_counts = sort_intersections_logically(
        intersection_counts=intersection_counts,
        groups=groups,
    )

    memberships = []
    counts = []

    for _, row in intersection_counts.iterrows():
        intersection = row["Intersection"]
        members = intersection.split("&")
        members = [group for group in groups if group in members]

        memberships.append(members)
        counts.append(row["count"])

    upset_data = from_memberships(
        memberships,
        data=counts,
    )

    upset_data = upset_data.reorder_levels(list(reversed(groups)))

    fig = plt.figure(figsize=(12, 7))

    upset = UpSet(
        upset_data,
        subset_size="sum",
        show_counts=True,
        sort_by="input",
        sort_categories_by=None,
    )

    upset.plot(fig=fig)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(
            os.path.join(outdir, f"upset_orthogroups.{ext}"),
            dpi=300,
            bbox_inches="tight",
        )

    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Create Venn and UpSet plots from OrthoFinder orthogroups."
    )

    parser.add_argument("--gene-counts", required=True, help="Orthogroups.GeneCount.tsv")
    parser.add_argument("--orthogroups", required=True, help="Orthogroups.tsv")
    parser.add_argument("--taxonomy", required=True, help="taxonomy_table.tsv")
    parser.add_argument("--groups", required=True, help="Comma-separated group names")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--colors-yaml", default=None, help="Optional colors.yaml file")
    parser.add_argument(
        "--keep-single-species",
        action="store_true",
        help=(
            "Keep orthogroups present in only one species among the selected groups. "
            "By default, orthogroups present in only one selected species are removed."
        ),
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    groups = [group.strip() for group in args.groups.split(",") if group.strip()]

    if len(groups) < 2:
        raise ValueError("--groups must contain at least two groups.")

    group_to_species = read_taxonomy(args.taxonomy, groups)
    gene_counts = read_gene_counts(args.gene_counts)
    orthogroups = read_orthogroups(args.orthogroups)

    presence_absence = build_presence_absence_table(
        gene_counts=gene_counts,
        group_to_species=group_to_species,
    )

    if not args.keep_single_species:
        presence_absence = remove_single_species_orthogroups(presence_absence)

    presence_absence.to_csv(
        os.path.join(args.outdir, "orthogroup_presence_absence_by_group.tsv"),
        sep="\t",
        index=False,
    )

    intersection_counts = count_intersections(presence_absence)
    intersection_counts = remove_none_intersection(intersection_counts)
    intersection_counts = sort_intersections_logically(
        intersection_counts=intersection_counts,
        groups=groups,
    )

    intersection_counts.to_csv(
        os.path.join(args.outdir, "intersection_counts.tsv"),
        sep="\t",
        index=False,
    )

    write_intersection_gene_files(
        orthogroups=orthogroups,
        presence_absence=presence_absence,
        group_to_species=group_to_species,
        outdir=args.outdir,
    )

    colors = load_colors(args.colors_yaml)

    plot_upset(
        presence_absence=presence_absence,
        groups=groups,
        outdir=args.outdir,
    )

    plot_venn(
        presence_absence=presence_absence,
        groups=groups,
        colors=colors,
        outdir=args.outdir,
    )

    print(f"Done. Results written in: {args.outdir}")


if __name__ == "__main__":
    main()