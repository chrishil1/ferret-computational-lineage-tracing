import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# Adjust parameters for output plots
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


def get_scaled_correlation_scores(filename, index_corr, index_low, index_high):
    """
    Load correlation scores and confidence intervals from a CSV file,
    scale them based on extreme quantiles, and clip values to [-1, 1].
    Scaling is performed using the maximum absolute value of the 1st and
    99th percentiles of the correlation scores to reduce the impact of
    outliers before clipping.

    Parameters:
    filename (str): Path to the CSV file containing correlation data.
    index_corr (int): Column index for the correlation scores.
    index_low (int): Column index for the lower confidence interval.
    index_high (int): Column index for the upper confidence interval.

    Returns:
    Scaled and clipped tuple (to [-1, 1]) correlation scores, lower bounds,
    and upper bounds, respectively.
    """

    # Read in csv file containing output scores and confidence intervals
    celltype_df = pd.read_csv(filename, index_col=0)
    scores = celltype_df.iloc[:, index_corr]
    ci_low = celltype_df.iloc[:, index_low]
    ci_high = celltype_df.iloc[:, index_high]

    # Calculate the 1st and 99th quantile and determine scaling factor
    quantile_01 = scores.quantile(0.01)
    quantile_99 = scores.quantile(0.99)
    scaling_factor = max(abs(quantile_01), abs(quantile_99))

    # Scale and clip the scores and confidence intervals
    scores_scaled = np.clip(scores / scaling_factor, -1, 1)
    ci_low_scaled = np.clip(ci_low / scaling_factor, -1, 1)
    ci_high_scaled = np.clip(ci_high / scaling_factor, -1, 1)

    # Return the pandas tuples
    return scores_scaled, ci_low_scaled, ci_high_scaled


def get_cumulative_overlap(scores1, scores2, thresholds):
    """
    Compute cumulative overlap count between two score arrays across thresholds.

    For each threshold, this function counts the number of positions where both
    input score arrays have values greater than or equal to the threshold.

    Overlap is defined element-wise as: (scores1 >= threshold) AND
    (scores2 >= threshold)

    Parameters:
    scores1 (array): First array of scores (e.g., numpy array or pandas Series).
    scores2 (array): Second array of scores, must be the same length as scores1.
    thresholds (array): Iterable of threshold values to compute overlap counts.

    Returns:
    Array of overlap counts corresponding to each threshold.
    """
    # Define list to contain cumulative overlaps
    cumulative_overlaps = []

    # Look for all instances in which the scores exceed the threshold
    for threshold in thresholds:
        overlap = ((scores1 >= threshold) & (scores2 >= threshold)).sum()
        cumulative_overlaps.append(overlap)

    # Convert to an array for the output
    return np.array(cumulative_overlaps)


def get_cumulative_gene_plot(ferret_scores, human_scores, celltype):
    """
    Generate and save a cumulative overlap plot of lineage correlation scores
    for human, ferret, and their intersection across a range of thresholds.

    This function computes cumulative gene overlaps at decreasing correlation
    thresholds and visualizes the results on a log-scaled y-axis. It also
    includes shaded confidence intervals for each curve.

    Parameters:
    ferret_scores: Tuple containing (scaled_scores, ci_low, ci_high) for ferret.
    human_scores: Tuple containing (scaled_scores, ci_low, ci_high) for human.
    celltype (str): Name of the cell type, used for title and output filename.

    Returns
    None - The function saves the generated plot as a PDF file.

    Notes
    -----
    - Thresholds range from 1 to 0 (50 evenly spaced values).
    - Overlaps are computed using `get_cumulative_overlap`.
    - The y-axis is log-scaled and the x-axis is inverted for visualization.
    - Shaded regions represent confidence intervals derived from the
      lower and upper bounds of the input scores.
    - Output file is saved to:
        'output_images/gene_plot_<celltype>.pdf'
    """

    # Define threshold as 50 evenly spaced values
    thresholds = np.linspace(1, 0, 50)

    # Compute cumulative overlap for human and ferret individually and combined
    cumulative_human = get_cumulative_overlap(
        human_scores[0], human_scores[0], thresholds
    )
    cumulative_ferret = get_cumulative_overlap(
        ferret_scores[0], ferret_scores[0], thresholds
    )
    cumulative_combined = get_cumulative_overlap(
        ferret_scores[0], human_scores[0], thresholds
    )

    # Define the figure and plot the colors as ferret, human, or overlap
    fig = plt.figure(figsize=(8, 6))
    plt.plot(thresholds, cumulative_human, label="human", color="red")
    plt.plot(thresholds, cumulative_ferret, label="ferret", color="blue")
    plt.plot(thresholds, cumulative_combined, label="human, ferret", color="green")

    # Get the high and low confidence intervals for the human line
    cumulative_human_ci_low = get_cumulative_overlap(
        human_scores[1], human_scores[1], thresholds
    )
    cumulative_human_ci_high = get_cumulative_overlap(
        human_scores[2], human_scores[2], thresholds
    )
    plt.fill_between(
        thresholds,
        cumulative_human_ci_low,
        cumulative_human_ci_high,
        color="red",
        alpha=0.1,
    )

    # Get the high and low confidence intervals for the ferret line
    cumulative_ferret_ci_low = get_cumulative_overlap(
        ferret_scores[1], ferret_scores[1], thresholds
    )
    cumulative_ferret_ci_high = get_cumulative_overlap(
        ferret_scores[2], ferret_scores[2], thresholds
    )
    plt.fill_between(
        thresholds,
        cumulative_ferret_ci_low,
        cumulative_ferret_ci_high,
        color="blue",
        alpha=0.1,
    )

    # Get the high and low confidence intervals for the combined line
    cumulative_combined_ci_low = get_cumulative_overlap(
        ferret_scores[1], human_scores[1], thresholds
    )
    cumulative_combined_ci_high = get_cumulative_overlap(
        ferret_scores[2], human_scores[2], thresholds
    )
    plt.fill_between(
        thresholds,
        cumulative_combined_ci_low,
        cumulative_combined_ci_high,
        color="green",
        alpha=0.1,
    )

    # Add dashed vertical line at threshold 0.7
    plt.axvline(x=0.7, color="black", linestyle="--")
    plt.yscale("log")
    plt.gca().invert_xaxis()

    # Add labels and legend
    plt.xlabel("Lineage Correlation")
    plt.ylabel("Cumulative Number of Genes")
    plt.title(f"Cumulative Number of Overlapping Lineage Drivers ({celltype})")
    plt.legend(loc="best")

    # Save the output figure
    fig.savefig(f"output_images/gene_plot_{celltype}.pdf", dpi=fig.dpi)


def get_combined_scores(ferret_score, human_score):
    """
    Merge ferret and human gene score data into a single aligned DataFrame.

    This function resets the indices of the input score objects, renames their
    columns for clarity, and performs an inner join on the gene identifiers to
    retain only genes present in both datasets.

    Parameters:
    ferret_score: Ferret gene scores with gene identifiers as the index.
    human_score: Human gene scores with gene identifiers as the index.

    Returns:
        A DataFrame with columns:
        - 'gene': gene identifiers
        - 'Ferret': ferret scores
        - 'Human': human scores
        containing only genes shared between both inputs.
    """

    # Index values converted into a gene column
    ferret_score = ferret_score.reset_index()
    human_score = human_score.reset_index()

    # Rename the columns for clarity
    ferret_score.columns = ["gene", "Ferret"]
    human_score.columns = ["gene", "Human"]

    # Merge the datasets on the gene column, keeping only common genes
    combined_scores = pd.merge(ferret_score, human_score, on="gene")

    # Rows that are NaN after merging are dropped
    combined_scores = combined_scores.dropna()

    return combined_scores


def get_lineage_heatmap(ferret_score, human_score, celltype):
    """
    Generate and save a heatmap of lineage correlation scores for selected genes
    across human and ferret datasets.

    This function merges human and ferret gene scores, filters genes into
    human-specific, ferret-specific, and overlapping categories based on a
    correlation threshold, and visualizes the selected genes as a heatmap.

    Parameters:
    ferret_score: Ferret gene scores with gene identifiers as the index.
    human_score: Human gene scores with gene identifiers as the index.
    celltype (str): Name of the cell type, used for title and output filename.

    Returns:
    None - The function saves the generated heatmap as a PDF file.
    """
    # Get the combined scores from the ferret and human arrays
    df = get_combined_scores(ferret_score, human_score)

    # Filter genes for the three categories (human, ferret, and overlap)
    human_specific = df[(df["Human"] > 0.7) & (df["Ferret"] <= 0.7)]
    ferret_specific = df[(df["Ferret"] > 0.7) & (df["Human"] <= 0.7)]
    human_ferret_overlap = df[(df["Human"] > 0.7) & (df["Ferret"] > 0.7)]

    # Create a combined DataFrame with only those filtered genes
    filtered_df = pd.concat([human_specific, ferret_specific, human_ferret_overlap])
    filtered_df.set_index("gene", inplace=True)
    filtered_df_transposed = filtered_df[["Human", "Ferret"]].T

    # Create a heatmap with lineage correlations
    fig = plt.figure(figsize=(14, 4))
    ax = sns.heatmap(
        filtered_df_transposed,
        cmap="coolwarm",
        annot=False,
        cbar_kws={"label": "Lineage Correlation"},
        vmin=-1,
        vmax=1,
        xticklabels=False,
    )
    ax.set(xlabel="")
    plt.title(f"Lineage Correlation Scores ({celltype})")

    # Save to our output figure
    fig.savefig(f"output_images/gene_plot_heatmap_{celltype}.pdf", dpi=fig.dpi)


if __name__ == "__main__":
    # Define the ferret scores by reading in from files
    ductal_ferret_scores = get_scaled_correlation_scores(
        "csv_files/ductal_proliferating_driver_genes_ferret.csv", 0, 3, 4
    )
    acinar_ferret_scores = get_scaled_correlation_scores(
        "csv_files/acinar_maturing_driver_genes_ferret.csv", 0, 3, 4
    )

    # Define the human scores by reading in from files
    ductal_human_scores = get_scaled_correlation_scores(
        "csv_files/ductal_immature_1_driver_genes_human.csv", 0, 3, 4
    )
    acinar_human_scores = get_scaled_correlation_scores(
        "csv_files/acinar_maturing_driver_genes_human.csv", 0, 3, 4
    )

    # Get the cumulative gene plot for the ductal and acinar celltypes
    get_cumulative_gene_plot(ductal_ferret_scores, ductal_human_scores, "Ductal")
    get_cumulative_gene_plot(acinar_ferret_scores, acinar_human_scores, "Acinar")

    # Get the lineage specific heatmap for the ductal and acinar celltypes
    get_lineage_heatmap(acinar_ferret_scores[0], acinar_human_scores[0], "Acinar")
    get_lineage_heatmap(ductal_ferret_scores[0], ductal_human_scores[0], "Ductal")
