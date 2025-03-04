#!/usr/bin/env python

import click
import matplotlib.pyplot as plt
import matplotlib_venn
import numpy as np
import pandas as pd
import seaborn as sns
import logomaker as lm

codons_nnn = {
    "TTT": "F",
    "CTT": "L",
    "ATT": "I",
    "GTT": "V",
    "TTC": "F",
    "CTC": "L",
    "ATC": "I",
    "GTC": "V",
    "TTA": "L",
    "CTA": "L",
    "ATA": "I",
    "GTA": "V",
    "TTG": "L",
    "CTG": "L",
    "ATG": "M",
    "GTG": "V",
    "TCT": "S",
    "CCT": "P",
    "ACT": "T",
    "GCT": "A",
    "TCC": "S",
    "CCC": "P",
    "ACC": "T",
    "GCC": "A",
    "TCA": "S",
    "CCA": "P",
    "ACA": "T",
    "GCA": "A",
    "TCG": "S",
    "CCG": "P",
    "ACG": "T",
    "GCG": "A",
    "TAT": "Y",
    "CAT": "H",
    "AAT": "N",
    "GAT": "D",
    "TAC": "Y",
    "CAC": "H",
    "AAC": "N",
    "GAC": "D",
    "TAA": "*",
    "CAA": "Q",
    "AAA": "K",
    "GAA": "E",
    "TAG": "*",
    "CAG": "Q",
    "AAG": "K",
    "GAG": "E",
    "TGT": "C",
    "CGT": "R",
    "AGT": "S",
    "GGT": "G",
    "TGC": "C",
    "CGC": "R",
    "AGC": "S",
    "GGC": "G",
    "TGA": "*",
    "CGA": "R",
    "AGA": "R",
    "GGA": "G",
    "TGG": "W",
    "CGG": "R",
    "AGG": "R",
    "GGG": "G",
}

CODONS_NNK = {
    codon: aa
    for codon, aa in codons_nnn.items()
    if codon.endswith("G") or codon.endswith("T")
}

# Subtract 1 to remove stop codons
# Background amino acid frequencies with NNK library
NNK_AA_FREQ = pd.Series(CODONS_NNK).value_counts() / (len(CODONS_NNK) - 1)
NNK_AA_FREQ = NNK_AA_FREQ.drop("*")

N_ALL_POSSIBLE_7MERS = 20**7


def gini_coefficient(x, w=None):
    # This is the Gini "Coefficient" for e.g. wealth inequality
    # Cribbed From https://stackoverflow.com/a/49571213/1628971
    # The rest of the code requires numpy arrays.
    x = np.asarray(x)
    if w is not None:
        w = np.asarray(w)
        sorted_indices = np.argsort(x)
        sorted_x = x[sorted_indices]
        sorted_w = w[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = np.cumsum(sorted_w, dtype=float)
        cumxw = np.cumsum(sorted_x * sorted_w, dtype=float)
        return np.sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) / (
            cumxw[-1] * cumw[-1]
        )
    else:
        sorted_x = np.sort(x)
        n = len(x)
        cumx = np.cumsum(sorted_x, dtype=float)
        # The above formula, with all weights equal to 1 simplifies to:
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n


def shannon(counts):
    probabilities = counts / counts.sum()
    shannon_index = -(probabilities * np.log2(probabilities)).sum()
    return shannon_index


def compute_global_stats(sequence_counts):
    nucleotide_21mer_counts = sequence_counts.groupby(level=1).sum()
    per_library_nucleotide_diversity_percentage = 100 * nucleotide_21mer_counts.apply(
        lambda x: (x > 0).sum() / x.sum()
    )
    per_library_nucleotide_diversity_percentage.name = "percent_nt_21mer_diversity"

    protein_7mer_counts = sequence_counts.groupby(level=0).sum()
    with_stop_codons = protein_7mer_counts.index.str.contains("\*")

    # Get per library stop codons
    protein_7mer_counts_no_stop_codons = protein_7mer_counts.loc[~with_stop_codons]
    per_library_n_stop_codons = protein_7mer_counts.apply(
        lambda x: x[x > 0].index.str.contains("\*").sum()
    )
    per_library_n_stop_codons.name = "n_stop_codons"
    per_library_stop_codons = 100 * protein_7mer_counts.apply(
        lambda x: x[x > 0].index.str.contains("\*").sum() / (x > 0).sum()
    )
    per_library_stop_codons.name = "percent_stop_codons"

    # Get per libary diversity percentage
    per_libary_diversity_percentage = 100 * protein_7mer_counts_no_stop_codons.apply(
        lambda x: (x > 0).sum() / x.sum()
    )
    per_libary_diversity_percentage.name = "percent_aa_7mer_diversity"

    per_library_shannon = protein_7mer_counts_no_stop_codons.apply(shannon)
    per_library_shannon.name = "shannon_index"

    # Get per library gini coefficients
    per_library_gini_coefficients = 100 * protein_7mer_counts_no_stop_codons.apply(
        gini_coefficient
    )
    per_library_gini_coefficients.name = "gini_coefficient"

    # Get per library number of k-mers
    per_library_n_7mers = (protein_7mer_counts_no_stop_codons > 0).sum()
    per_library_n_7mers.name = "n_7mers"

    # get percentage of total possible 7-mers
    per_libary_percent_7mers = 100 * per_library_n_7mers / N_ALL_POSSIBLE_7MERS
    per_libary_percent_7mers.name = "percent_all_possible_7mers"

    # Concatenate everything for per library stats
    per_library_stats = pd.concat(
        [
            per_libary_diversity_percentage,
            per_library_nucleotide_diversity_percentage,
            per_library_n_stop_codons,
            per_library_stop_codons,
            per_library_gini_coefficients,
            per_library_shannon,
            per_library_n_7mers,
            per_libary_percent_7mers,
        ],
        axis=1,
    ).T
    return per_library_stats, protein_7mer_counts_no_stop_codons


def venn_diagram(
    protein_7mer_counts_no_stop_codons, prefix="", saveas="venn_diagram.png"
):
    if prefix:
        saveas = f"{prefix}__{saveas}"

    n_libraries = len(protein_7mer_counts_no_stop_codons.columns)
    if n_libraries == 2:
        venn = matplotlib_venn.venn2
    elif n_libraries == 3:
        venn = matplotlib_venn.venn3
    else:
        raise ValueError(
            f"Number of libraries is {n_libraries}, can only make venn "
            "diagram for 2 or 3 libraries at a time"
        )
    sets = list(protein_7mer_counts_no_stop_codons.apply(lambda x: set(x.index[x > 0])))
    set_labels = protein_7mer_counts_no_stop_codons.columns.tolist()
    venn(sets, set_labels)
    fig = plt.gcf()
    fig.savefig(saveas)


def _make_position_frequency_matrices(protein_7mer_counts_no_stop_codons):
    observed_7mers_expanded = pd.DataFrame.from_records(
        protein_7mer_counts_no_stop_codons.index.map(list).values, columns=range(1, 8)
    )

    position_freq_matrices = {
        x: pd.DataFrame() for x in protein_7mer_counts_no_stop_codons.columns
    }

    # Count frequency of each amino acid
    for i, col in observed_7mers_expanded.iteritems():
        df = protein_7mer_counts_no_stop_codons.groupby(col.values).sum()

        for library, pwm in position_freq_matrices.items():
            position_freq_matrices[library] = pd.concat([pwm, df[library]], axis=1)

    for library, pwm in position_freq_matrices.items():
        position_freq_matrices[library].columns = observed_7mers_expanded.columns
    return position_freq_matrices


def _make_position_probability_matrices(position_freq_matrices):
    # Make position probability matrices
    position_probability_matrices = {}
    for library, pwm in position_freq_matrices.items():
        ppm = pwm / pwm.sum()

        position_probability_matrices[library] = ppm

        fig, ax = plt.subplots()
        ax = sns.heatmap(ppm.T, center=0.05, cmap="RdBu_r")
        ax.set_title(library)
        fig.savefig(f"{library}__position_probability_matrix.png")

    return position_probability_matrices


def _make_position_weight_matricies(
    position_probability_matrices, background_amino_acid_frequencies=NNK_AA_FREQ
):
    position_weight_matrices = {}

    for library, ppm in position_probability_matrices.items():
        # print(f'--- {library} ---')
        # background_probability = ppm.mean().mean()
        weight_matrix = np.log2(ppm.divide(background_amino_acid_frequencies, axis=0))
        position_weight_matrices[library] = weight_matrix

        fig, ax = plt.subplots()
        ax = sns.heatmap(weight_matrix.T, cmap="RdBu_r", center=0)
        ax.set_title(library)
        fig.savefig(f"{library}__position_weight_matrix.png")

    return position_weight_matrices


def _compute_entropy(position_probability_matrices):
    entropies = {}

    for library, pwm in position_probability_matrices.items():
        # print(f'--- {library} ---')
        library_entropy = pwm.apply(shannon)
        # print(library_entropy)
        entropies[library] = library_entropy

        # fig, ax = plt.subplots()
        # ax = sns.heatmap(library_entropy.T)
        # ax.set_title(library)

    entropies = pd.DataFrame(entropies)
    entropies.index = range(1, 8)
    entropies.index.name = "position"
    return entropies


def barplot_entropies_per_position(entropies, prefix):
    entropies_tidy = entropies.unstack().reset_index()
    entropies_tidy = entropies_tidy.rename(columns={"level_0": "library", 0: "entropy"})

    g = sns.catplot(
        data=entropies_tidy,
        y="entropy",
        x="position",
        hue="library",
        palette="viridis",
        # order=aa_order,
        # hue_order=list(range(1, 8)),
        aspect=2,
        kind="bar",
        height=3,
    )
    g.savefig(f"{prefix}__entropies_barplot.png")


def barplot_position_probability(position_probability_matrices, prefix):
    dfs = []
    for library, matrix in position_probability_matrices.items():
        df = matrix.stack().reset_index()
        df["library"] = library
        dfs.append(df)
    position_prob_df = pd.concat(dfs)
    position_prob_df = position_prob_df.rename(
        columns={"level_0": "aa", "level_1": "position", 0: "probability"}
    )

    # PLot per-position enrichment
    g = sns.catplot(
        data=position_prob_df,
        kind="bar",
        x="aa",
        y="probability",
        col="library",
        row="position",
        height=2,
        aspect=2,
    )
    g.savefig(
        f"{prefix}__position_probability_barplot_separate_libraries_and_positions.pdf"
    )

    # Plot summary across all positions
    g = sns.catplot(
        data=position_prob_df,
        kind="bar",
        x="aa",
        y="probability",
        col="library",
        height=2,
        aspect=2,
    )
    g.savefig(
        f"{prefix}__position_probability_barplot_separate_libraries_mean_positions.pdf"
    )

    # Compare per-position amino acid frequency library-by-library
    g = sns.catplot(
        data=position_prob_df,
        kind="bar",
        x="aa",
        y="probability",
        hue="library",
        height=2,
        aspect=2,
        row="position",
    )
    g.savefig(f"{prefix}__position_probability_barplot_compare_libraries.pdf")

    # Take the mean across positions
    g = sns.catplot(
        data=position_prob_df,
        kind="bar",
        x="aa",
        y="probability",
        hue="library",
        height=2,
        aspect=2,
        errwidth=1,
    )
    g.savefig(
        f"{prefix}__position_probability_barplot_compare_libraries_mean_across_positions.pdf"
    )


def make_position_weight_matrices(protein_7mer_counts_no_stop_codons, prefix):
    position_freq_matrices = _make_position_frequency_matrices(
        protein_7mer_counts_no_stop_codons
    )
    position_probability_matrices = _make_position_probability_matrices(
        position_freq_matrices
    )

    barplot_position_probability(position_probability_matrices, prefix)
    entropies = _compute_entropy(position_probability_matrices)

    fig, ax = plt.subplots()
    sns.heatmap(entropies.T)
    fig.savefig(f"{prefix}__entropies.png")

    barplot_entropies_per_position(entropies, prefix)

    position_weight_matrices = _make_position_weight_matricies(
        position_probability_matrices
    )
    return position_weight_matrices


def make_tidy_position_weights(position_weight_matrices):
    dfs = []

    for lib, data2d in position_weight_matrices.items():
        df = data2d.stack().reset_index()
        df = df.rename(columns={"level_0": "aa", "level_1": "pos", 0: "weight"})
        df["library"] = lib
        dfs.append(df)
    position_weights_tidy = pd.concat(dfs)
    position_weights_tidy.head()
    return position_weights_tidy


def barplot_position_weights(position_weights_tidy, barplot_png, col_wrap=3):
    aa_mean_weights = position_weights_tidy.groupby("aa").weight.mean().sort_values()
    aa_order = aa_mean_weights.index

    g = sns.catplot(
        data=position_weights_tidy,
        col="library",
        y="weight",
        x="aa",
        hue="pos",
        palette="viridis",
        order=aa_order,
        hue_order=list(range(1, 8)),
        aspect=2,
        col_wrap=col_wrap,
        kind="bar",
        height=3,
    )
    g.savefig(barplot_png)

    order = position_weights_tidy.groupby("aa").weight.mean().sort_values().index[::-1]
    g = sns.catplot(
        data=position_weights_tidy,
        x="aa",
        y="weight",
        col="pos",
        kind="bar",
        order=order,
        aspect=1.5,
        height=2,
        col_wrap=4,
    )
    g.savefig("position_weights_across_libraries.pdf")


def make_data_for_bubble_plot(sequence_counts):
    aa = sequence_counts.groupby("aa_7mer").sum()

    # Remove all stop codons
    aa_no_stop = aa.loc[~aa.index.str.contains("*", regex=False)]

    # Create x-axis index/position/coordinate for each amino acid
    aa_index = pd.Series(
        range(len(aa_no_stop)), index=aa_no_stop.index, name="aa_7mer_index"
    )

    # Get tidy series of number of reads
    aa_n_reads_tidy = aa_no_stop.stack()
    aa_n_reads_tidy.name = "n_reads"

    # Get tidy series of percent reads
    aa_no_stop_pct = 100 * aa_no_stop / aa_no_stop.sum()
    aa_no_stop_pct_tidy = aa_no_stop_pct.stack()
    aa_no_stop_pct_tidy.name = "percent_reads"

    # Concatenate percent and number of reads. Join with the x-axis index/coordinate of each amino acid
    aa_combined = pd.concat([aa_no_stop_pct_tidy, aa_n_reads_tidy], axis=1).join(
        aa_index
    )
    aa_combined = aa_combined.reset_index()
    aa_combined = aa_combined.rename(columns={"level_1": "library"})
    aa_combined = aa_combined.query("n_reads > 0")

    return aa_combined


def bubble_plot(
    aa_combined, col_order=None, bubble_plot_png="bubble_plot.png", col_wrap=4
):
    # Create plot
    aa_combined["percent_reads_log10"] = np.log10(aa_combined["percent_reads"])
    col_order = sorted(list(aa_combined["library"].unique()))
    g = sns.relplot(
        x="aa_7mer_index",
        y="percent_reads_log10",
        size="n_reads",
        data=aa_combined,
        hue="library",
        col="library",
        kind="scatter",
        height=3,
        hue_order=col_order,
        col_order=col_order,
        col_wrap=col_wrap,
    )
    # g.set(yscale="log")
    g.savefig(bubble_plot_png)

    return g


def barplot_top_kmers(bubble_plot_data, top_kmers_barplot_png, col_wrap=4):
    top10 = bubble_plot_data.groupby("library").apply(
        lambda x: x.nlargest(10, "n_reads")
    )
    print(top10)
    g = sns.catplot(
        y="aa_7mer",
        x="percent_reads",
        data=top10,
        sharey=False,
        col="library",
        kind="bar",
        height=4,
        hue="library",
        dodge=False,
        col_wrap=col_wrap,
    )
    g.set(xscale="log")

    for (i, j, k), data_ijk in g.facet_data():
        ax = g.facet_axis(i, j)
        for i, (index, row) in enumerate(data_ijk.iterrows()):
            ax.annotate(
                f"{row.n_reads}",
                (row.percent_reads, i),
                ha="left",
                va="center",
                # color="white",
            )
    g.fig.tight_layout()
    g.savefig(top_kmers_barplot_png)


def get_col_wrap(position_weights_tidy):
    """Get number of columns to plot for seaborn catplots"""
    n_libraries = position_weights_tidy["library"].nunique()
    if n_libraries < 3:
        col_wrap = n_libraries
    else:
        col_wrap = 3
    return col_wrap


def make_logos(position_weights_tidy):
    for library, df in position_weights_tidy.groupby("library"):
        matrix = df.pivot(index="pos", columns="aa", values="weight")
        matrix = matrix[matrix > 0].fillna(0)
        lm.Logo(matrix)
        fig = plt.gcf()
        fig.savefig(f"{library}_enriched_sequence_logo.png")


@click.command()
@click.argument("sequence_counts_csv", type=click.Path())
@click.option("--aa-diversity-stats", default="aa-diversity-stats.csv")
@click.option("--bubble-plot-png", default="bubble_plot.png")
@click.option("--position-weights-barplot-png", default="position_weights_barplot.png")
@click.option("--position-weights-csv", default="position_weights.csv")
@click.option("--top-kmers-barplot-png", default="top_kmers_barplot.png")
@click.option("--prefix", default="")
def main(
    sequence_counts_csv,
    aa_diversity_stats,
    bubble_plot_png,
    position_weights_barplot_png,
    position_weights_csv,
    top_kmers_barplot_png,
    prefix,
):
    if prefix:
        bubble_plot_png = f"{prefix}__{bubble_plot_png}"
        aa_diversity_stats = f"{prefix}__{aa_diversity_stats}"
        position_weights_barplot_png = f"{prefix}__{position_weights_barplot_png}"
        top_kmers_barplot_png = f"{prefix}__{top_kmers_barplot_png}"

    sequence_counts = pd.read_csv(sequence_counts_csv, index_col=[0, 1])

    per_library_stats, protein_7mer_counts_no_stop_codons = compute_global_stats(
        sequence_counts
    )
    # Take the transpose so the samples are the row names to make it easy for MultiQC
    per_library_stats.T.to_csv(aa_diversity_stats)

    try:
        venn_diagram(protein_7mer_counts_no_stop_codons, prefix)
    except ValueError:
        print("Not the right number of samples for a Venn diagram --> skipping")
        pass

    # # Make bubble plot
    bubble_plot_data = make_data_for_bubble_plot(sequence_counts)

    col_wrap = get_col_wrap(bubble_plot_data)

    # Bubble plot stopped working and I don't know why
    bubble_plot(bubble_plot_data, bubble_plot_png=bubble_plot_png, col_wrap=col_wrap)
    barplot_top_kmers(bubble_plot_data, top_kmers_barplot_png, col_wrap=col_wrap)

    position_weight_matrices = make_position_weight_matrices(
        protein_7mer_counts_no_stop_codons, prefix
    )

    position_weights_tidy = make_tidy_position_weights(position_weight_matrices)
    position_weights_tidy.to_csv(position_weights_csv)

    # Make barplot of position weights, sorted by amino acid
    barplot_position_weights(
        position_weights_tidy, position_weights_barplot_png, col_wrap=col_wrap
    )

    # Make sequence logos of enriched k-mers
    make_logos(position_weights_tidy)


if __name__ == "__main__":
    main()
