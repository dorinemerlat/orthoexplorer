#!/usr/bin/env python3

import argparse
import math
import re
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable


NAMESPACE_LABELS = {
    "BP": "Biological Process",
    "MF": "Molecular Function",
    "CC": "Cellular Component",
}

NAMESPACE_ORDER = ["BP", "MF", "CC"]


def parse_args():
    parser = argparse.ArgumentParser(description="Create a faceted GO enrichment dotplot from a GOATOOLS result table.")
    parser.add_argument("--input", required=True, type=Path, help="GO enrichment result table.")
    parser.add_argument("--output-prefix", required=True, type=Path, help="Output prefix for PNG, PDF and SVG figures.")
    parser.add_argument("--title", help="Optional figure title.")
    parser.add_argument("--max-terms", type=int, default=10, help="Maximum number of GO terms shown per namespace.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05, help="Maximum adjusted p-value retained.")
    parser.add_argument("--min-fold-enrichment", type=float, default=1.0, help="Minimum fold enrichment retained.")
    parser.add_argument("--label-width", type=int, default=45, help="Maximum number of characters per GO-term label line.")
    parser.add_argument("--width", type=float, default=11, help="Figure width in inches.")
    parser.add_argument("--panel-height", type=float, default=3.5, help="Height of each namespace panel in inches.")
    return parser.parse_args()


def normalise_column_name(column):
    """Return a lowercase column name with normalised separators."""
    return re.sub(r"[^a-z0-9]+", "_", str(column).strip().lower()).strip("_")


def find_column(df, candidates, required=True):
    """Return the first matching column from a list of accepted names."""
    normalised_columns = {normalise_column_name(column): column for column in df.columns}

    for candidate in candidates:
        normalised_candidate = normalise_column_name(candidate)
        if normalised_candidate in normalised_columns:
            return normalised_columns[normalised_candidate]

    if required:
        expected = ", ".join(candidates)
        raise ValueError(f"Missing required column. Expected one of: {expected}")

    return None


def parse_ratio(value):
    """Parse a GOATOOLS ratio represented as 'count/total'."""
    if pd.isna(value):
        return math.nan, math.nan

    value = str(value).strip()
    match = re.fullmatch(r"(\d+)\s*/\s*(\d+)", value)

    if not match:
        raise ValueError(f"Invalid ratio value: {value!r}. Expected the format 'count/total'.")

    return int(match.group(1)), int(match.group(2))


def normalise_namespace(value):
    """Convert GO namespace values to BP, MF or CC."""
    value = str(value).strip().lower()

    namespace_aliases = {
        "bp": "BP",
        "biological_process": "BP",
        "biological process": "BP",
        "mf": "MF",
        "molecular_function": "MF",
        "molecular function": "MF",
        "cc": "CC",
        "cellular_component": "CC",
        "cellular component": "CC",
    }

    return namespace_aliases.get(value, str(value).strip().upper())


def wrap_label(value, width):
    """Wrap a GO-term name over multiple lines."""
    return "\n".join(textwrap.wrap(str(value), width=width, break_long_words=False, break_on_hyphens=False))


def read_enrichment_table(path):
    """Read a tab-separated GO enrichment result table."""
    try:
        df = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError as error:
        raise ValueError(f"The enrichment table is empty: {path}") from error

    # Normalize column names produced by GOATOOLS.
    df.columns = df.columns.str.strip().str.removeprefix("#").str.strip()

    if df.empty:
        raise ValueError(f"The enrichment table contains no results: {path}")

    return df


def prepare_enrichment_results(df, fdr_threshold, min_fold_enrichment):
    """Standardise GOATOOLS columns and calculate plotting statistics."""
    go_column = find_column(df, ["GO", "go_id", "go"])
    namespace_column = find_column(df, ["NS", "namespace", "domain"])
    name_column = find_column(df, ["name", "go_name", "term"])
    study_ratio_column = find_column(df, ["ratio_in_study", "study_ratio"])
    population_ratio_column = find_column(df, ["ratio_in_pop", "ratio_in_population", "background_ratio"])
    fdr_column = find_column(df, ["p_fdr_bh", "fdr", "padj", "p_adjusted"])
    enrichment_column = find_column(df, ["enrichment", "status"], required=False)

    study_ratios = df[study_ratio_column].apply(parse_ratio)
    population_ratios = df[population_ratio_column].apply(parse_ratio)

    results = pd.DataFrame({
        "go_id": df[go_column].astype(str),
        "go_name": df[name_column].astype(str),
        "namespace": df[namespace_column].apply(normalise_namespace),
        "study_count": study_ratios.str[0],
        "study_total": study_ratios.str[1],
        "background_count": population_ratios.str[0],
        "background_total": population_ratios.str[1],
        "fdr": pd.to_numeric(df[fdr_column], errors="coerce"),
    })

    results["study_fraction"] = results["study_count"] / results["study_total"]
    results["background_fraction"] = results["background_count"] / results["background_total"]
    results["fold_enrichment"] = results["study_fraction"] / results["background_fraction"]

    if enrichment_column is not None:
        enrichment = df[enrichment_column].astype(str).str.strip().str.lower()
        enriched_values = {"e", "enriched", "enrichment", "overrepresented", "over-represented"}
        results = results.loc[enrichment.isin(enriched_values)].copy()

    results = results.loc[
        results["namespace"].isin(NAMESPACE_ORDER)
        & results["fdr"].notna()
        & results["fold_enrichment"].notna()
        & results["fold_enrichment"].replace([float("inf"), -float("inf")], pd.NA).notna()
        & (results["fdr"] <= fdr_threshold)
        & (results["fdr"] > 0)
        & (results["fold_enrichment"] >= min_fold_enrichment)
        & (results["study_count"] > 0)
    ].copy()

    results["minus_log10_fdr"] = -results["fdr"].apply(math.log10)

    return results


def select_top_terms(results, max_terms):
    """Keep the most significant GO terms within each namespace."""
    selected = []

    for namespace in NAMESPACE_ORDER:
        subset = results.loc[results["namespace"] == namespace].copy()
        subset = subset.sort_values(
            ["fdr", "fold_enrichment", "study_count"],
            ascending=[True, False, False],
        ).head(max_terms)

        selected.append(subset)

    return pd.concat(selected, ignore_index=True)


def build_size_legend(ax, minimum, maximum):
    """Add a compact legend representing orthogroup counts."""
    if minimum == maximum:
        size_values = [int(round(minimum))]
    else:
        size_values = sorted({
            int(round(minimum)),
            int(round((minimum + maximum) / 2)),
            int(round(maximum)),
        })

    minimum_area = 40
    maximum_area = 260

    def scale_size(value):
        if minimum == maximum:
            return (minimum_area + maximum_area) / 2

        return minimum_area + (value - minimum) / (maximum - minimum) * (maximum_area - minimum_area)

    handles = [
        Line2D(
            [],
            [],
            linestyle="",
            marker="o",
            markersize=math.sqrt(scale_size(value)),
            markerfacecolor="grey",
            markeredgecolor="grey",
            label=str(value),
        )
        for value in size_values
    ]

    return ax.legend(
        handles=handles,
        title="Orthogroups",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )


def plot_enrichment(results, output_prefix, title, max_terms, label_width, width, panel_height):
    """Create one vertically stacked dotplot panel per GO namespace."""
    selected = select_top_terms(results, max_terms)

    if selected.empty:
        for extension in ("png", "pdf", "svg"):
            output_prefix.with_suffix(f".{extension}").touch()
        return

    if selected.empty:
        raise ValueError("No GO terms passed the enrichment and FDR filters.")

    sns.set_theme(style="whitegrid", context="notebook")

    figure, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(width, panel_height * 3),
        sharex=True,
        gridspec_kw={"hspace": 0.28},
    )

    colour_min = selected["minus_log10_fdr"].min()
    colour_max = selected["minus_log10_fdr"].max()

    if colour_min == colour_max:
        colour_max = colour_min + 1e-9

    colour_norm = Normalize(vmin=colour_min, vmax=colour_max)
    size_min = selected["study_count"].min()
    size_max = selected["study_count"].max()

    for axis, namespace in zip(axes, NAMESPACE_ORDER):
        subset = selected.loc[selected["namespace"] == namespace].copy()

        axis.set_title(NAMESPACE_LABELS[namespace], loc="left", fontweight="bold")
        axis.set_ylabel("")
        axis.grid(axis="x", visible=True)
        axis.grid(axis="y", visible=False)

        if subset.empty:
            axis.text(
                0.5,
                0.5,
                "No significant enriched GO terms",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
            axis.set_yticks([])
            continue

        subset = subset.sort_values(
            ["fold_enrichment", "fdr"],
            ascending=[True, False],
        )
        subset["display_name"] = subset["go_name"].apply(lambda value: wrap_label(value, label_width))

        sns.scatterplot(
            data=subset,
            x="fold_enrichment",
            y="display_name",
            hue="minus_log10_fdr",
            size="study_count",
            palette="coolwarm",
            hue_norm=colour_norm,
            sizes=(40, 260),
            legend=False,
            edgecolor="black",
            linewidth=0.3,
            ax=axis,
        )

        axis.tick_params(axis="y", labelsize=9)

    axes[-1].set_xlabel("Fold enrichment")

    for axis in axes[:-1]:
        axis.set_xlabel("")

    colour_mapper = ScalarMappable(norm=colour_norm, cmap=sns.color_palette("coolwarm", as_cmap=True))
    colour_mapper.set_array([])

    colour_bar = figure.colorbar(
        colour_mapper,
        ax=axes,
        fraction=0.025,
        pad=0.12,
        aspect=30,
    )
    colour_bar.set_label("−log10(FDR)")

    size_legend = build_size_legend(axes[0], size_min, size_max)
    axes[0].add_artist(size_legend)

    if title:
        figure.suptitle(title, fontsize=14, fontweight="bold", y=0.995)

    figure.subplots_adjust(left=0.36, right=0.80, top=0.95 if title else 0.98, bottom=0.07)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    for extension in ("png", "pdf", "svg"):
        output_path = output_prefix.with_suffix(f".{extension}")
        figure.savefig(output_path, dpi=300, bbox_inches="tight")

    plt.close(figure)


def main():
    args = parse_args()

    if args.max_terms < 1:
        sys.exit("Error: --max-terms must be greater than zero.")

    if not 0 < args.fdr_threshold <= 1:
        sys.exit("Error: --fdr-threshold must be between 0 and 1.")

    try:
        df = read_enrichment_table(args.input)
        results = prepare_enrichment_results(
            df,
            fdr_threshold=args.fdr_threshold,
            min_fold_enrichment=args.min_fold_enrichment,
        )
        plot_enrichment(
            results=results,
            output_prefix=args.output_prefix,
            title=args.title,
            max_terms=args.max_terms,
            label_width=args.label_width,
            width=args.width,
            panel_height=args.panel_height,
        )
    except (OSError, ValueError) as error:
        sys.exit(f"Error: {error}")


if __name__ == "__main__":
    main()