import pandas as pd
import scanpy as sc
import cellrank as cr
import numpy as np
import matplotlib as mpl
import os

# Ensure text is stored as editable fonts, not vector paths
mpl.rcParams["pdf.fonttype"] = 42  # TrueType font
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

# Define palette to be used in plotting
human_palette = {
    "Acinar Immature": "#e64463",
    "Acinar Maturing": "#62be4a",
    "Acinar Progenitor": "#9255cb",
    "Acinar Proliferating": "#b2b336",
    "Acinar Secretory": "#5069d7",
    "Alpha": "#dc8f2f",
    "B lymphocyte": "#d579e1",
    "Beta": "#4c8a35",
    "Delta": "#ca42a9",
    "Dendritic Cell": "#54c081",
    "Ductal ECM": "#d83f86",
    "Ductal Immature": "#52bfad",
    "Ductal Maturing": "#cc4724",
    "Ductal Proliferating": "#4ab4dd",
    "Ducto-Acinar": "#dc6c4d",
    "Endothelial": "#5982c5",
    "Endothelial Lymphatic": "#cca15d",
    "Erythroblast": "#9384d9",
    "Mesenchyme Proliferating": "#727229",
    "Mesenchyme Subtype 1": "#ae57a1",
    "Mesenchyme Subtype 2": "#9cb26a",
    "Mesenchyme Subtype 3": "#755298",
    "Mesenchyme Subtype 4": "#37835d",
    "Mono-Macrophage": "#b83849",
    "Neuron": "#ca96d2",
    "Pancreatic Progenitor": "#9a5e2c",
    "Phagocytic Macrophage": "#e47da1",
    "Schwann": "#e5907a",
    "Stromal": "#9f4970",
    "T lymphocyte": "#aa5659",
}


def ferret_acinar():
    """
    Runs CellRank trajectory analysis on ferret acinar pancreatic cell data.

    Loads pseudotime, lineage weights, and AnnData from disk, computes a
    PseudotimeKernel transition matrix, fits a GPCCA estimator with 5 terminal
    states (Acinar Proliferating, Acinar Maturing, Ductal Immature, Pancreatic
    Progenitor, Endothelial Lymphatic), and computes fate probabilities and
    lineage driver genes for each terminal state. Saves all plots and driver
    gene tables to disk.

    Reads:
        pseudotime_objects/ferret_acinar_s9e10_pseudotime.csv
        pseudotime_objects/ferret_acinar_s9e10_lineage_weights.csv
        pseudotime_objects/ferret_acinar.h5ad

    Writes:
        figures/projection_ferret_acinar.svg
        Several .csv files that detail driver genes for each terminal state
        Several .pdf files that highlight different parts of the fate
        determination
    """

    # Set the output save file directories
    figdir = "figures/12_cellrank_plots/ferret_acinar"
    os.makedirs(figdir, exist_ok=True)
    sc.settings.figdir = figdir
    cr.settings.figdir = figdir

    # Load in pseudotime and weights for this specific acinar lineage
    ferret_acinar_pseudotime = pd.read_csv(
        "pseudotime_objects/ferret_acinar_s9e10_pseudotime.csv", index_col=0
    )
    ferret_acinar_weights = pd.read_csv(
        "pseudotime_objects/ferret_acinar_s9e10_lineage_weights.csv", index_col=0
    )

    # Load in the h5ad file from the ferret Seurat object
    adata_ferret_acinar = sc.read_h5ad("pseudotime_objects/ferret_acinar.h5ad")

    # Specify that we are looking at lineage 2 for this slingshot pseudotime
    adata_ferret_acinar.obs["pseudotime"] = ferret_acinar_pseudotime["Lineage2"]
    adata_ferret_acinar.obs["weights"] = ferret_acinar_weights["Lineage2"]

    # Remove any NaNs, recategorize the majority celltype as a category
    adata_ferret_acinar = adata_ferret_acinar[
        ~np.isnan(adata_ferret_acinar.obs["pseudotime"])
    ]
    adata_ferret_acinar.obs["final_renamed_celltype"] = adata_ferret_acinar.obs[
        "final_renamed_celltype"
    ].astype("category")

    # Calculate the nearest neighbors
    sc.pp.neighbors(adata_ferret_acinar, n_pcs=15, n_neighbors=15)

    # Instantiate a pseudotime kernel and calculate transition matrix
    pk = cr.kernels.PseudotimeKernel(adata_ferret_acinar, time_key="pseudotime")
    pk.compute_transition_matrix()

    # Visualize transition matrix
    pk.plot_projection(
        color="final_renamed_celltype",
        save=os.path.join(figdir, "projection_ferret_acinar.svg"),
        palette=human_palette,
        legend_loc="right",
        figsize=(8, 6),
        add_margin=0.02,
    )

    # Instantiate the GPCCA esimator and fit to the clusters using celltype
    g = cr.estimators.GPCCA(pk)
    g.fit(cluster_key="final_renamed_celltype", n_states=5)

    # Manually set the terminal states after fit, calculate fate probabilities
    g.set_terminal_states(
        states=[
            "Acinar Proliferating",
            "Acinar Maturing",
            "Ductal Immature",
            "Pancreatic Progenitor",
            "Endothelial Lymphatic",
        ]
    )
    g.compute_fate_probabilities()

    # Series of plots to visualize macrostate and fate probabilities
    g.plot_macrostates(
        which="terminal",
        legend_loc="right",
        s=100,
        save=os.path.join(figdir, "ferret_acinar_macrostates_terminal.pdf"),
    )
    g.plot_macrostates(
        which="terminal",
        discrete=False,
        save=os.path.join(figdir, "ferret_acinar_macrostates_terminal_discrete.pdf"),
    )
    g.plot_fate_probabilities(
        same_plot=False,
        save=os.path.join(figdir, "ferret_acinar_fate_probabilities_separate.pdf"),
    )
    cr.pl.circular_projection(
        adata_ferret_acinar,
        keys=["final_renamed_celltype"],
        legend_loc="right",
        save="ferret_acinar_circular_projection.pdf",
    )

    # Name all celltypes to be driving factors when calculating lineage drivers
    driver_clusters = [
        "Acinar Proliferating",
        "Acinar Immature",
        "Acinar Maturing",
        "Ductal Immature",
        "Endothelial Lymphatic",
        "Ductal Proliferating",
        "Pancreatic Progenitor",
        "Acinar Progenitor",
        "Alpha",
        "Mesenchyme Subtype 1",
        "B lymphocyte",
    ]

    # Calculate lineage drivers with ductal immature as the terminal state
    ductal_immature_df = g.compute_lineage_drivers(
        lineages=["Ductal Immature"],
        cluster_key="final_renamed_celltype",
        clusters=driver_clusters,
    )
    adata_ferret_acinar.obs["fate_probabilities_ductal_immature"] = (
        g.fate_probabilities["Ductal Immature"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_ferret_acinar,
        basis="umap",
        color=["fate_probabilities_ductal_immature"]
        + list(ductal_immature_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_ductal_immature_driver_genes_ferret_acinar.pdf",
    )

    # Calculate lineage drivers with acinar maturing as the terminal state
    acinar_mat_df = g.compute_lineage_drivers(
        lineages=["Acinar Maturing"],
        cluster_key="final_renamed_celltype",
        clusters=[
            "Acinar Proliferating",
            "Acinar Immature",
            "Acinar Maturing",
            "Endothelial Lymphatic",
            "Pancreatic Progenitor",
            "Acinar Progenitor",
        ],
    )
    adata_ferret_acinar.obs["fate_probabilities_acinar_maturing"] = (
        g.fate_probabilities["Acinar Maturing"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_ferret_acinar,
        basis="umap",
        color=["fate_probabilities_acinar_maturing"] + list(acinar_mat_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_acinar_maturing_driver_genes_ferret_acinar.pdf",
    )
    acinar_mat_df.to_csv(
        "csv_files/acinar_maturing_driver_genes_ferret.csv", index=True
    )

    # Calculate lineage drivers with acinar proliferating as the terminal state
    acinar_prolif_df = g.compute_lineage_drivers(
        lineages=["Acinar Proliferating"],
        cluster_key="final_renamed_celltype",
        clusters=driver_clusters,
    )
    adata_ferret_acinar.obs["fate_probabilities_acinar_proliferating"] = (
        g.fate_probabilities["Acinar Proliferating"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_ferret_acinar,
        basis="umap",
        color=["fate_probabilities_acinar_proliferating"]
        + list(acinar_prolif_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_acinar_proliferating_driver_genes_ferret_acinar.pdf",
    )

    # Calculate lineage drivers with pancreatic progenitor as the terminal state
    pancreatic_progenitor_df = g.compute_lineage_drivers(
        lineages=["Pancreatic Progenitor"],
        cluster_key="final_renamed_celltype",
        clusters=driver_clusters,
    )
    adata_ferret_acinar.obs["fate_probabilities_pancreatic_progenitor"] = (
        g.fate_probabilities["Pancreatic Progenitor"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_ferret_acinar,
        basis="umap",
        color=["fate_probabilities_pancreatic_progenitor"]
        + list(pancreatic_progenitor_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_pancreatic_progenitor_driver_genes_ferret_acinar.pdf",
    )

    # Calculate lineage drivers with endothelial lympathic as the terminal state
    endothelial_lympathic_df = g.compute_lineage_drivers(
        lineages=["Endothelial Lymphatic"],
        cluster_key="final_renamed_celltype",
        clusters=driver_clusters,
    )
    adata_ferret_acinar.obs["fate_probabilities_endothelial_lympathic"] = (
        g.fate_probabilities["Endothelial Lymphatic"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_ferret_acinar,
        basis="umap",
        color=["fate_probabilities_endothelial_lympathic"]
        + list(endothelial_lympathic_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_endothelial_lympathic_driver_genes_ferret_acinar.pdf",
    )

    # Plot the well-associated gene sets and their distribution
    sc.pl.embedding(
        adata_ferret_acinar,
        basis="umap",
        color=[
            "ATF3",
            "BTG2",
            "CEL",
            "CELA3A",
            "CPA1",
            "CPA2",
            "GATA4",
            "NFKBIZ",
            "PTPN3",
            "RBPJL",
            "SERPINI2",
            "XBP1",
        ],
        color_map="viridis",
        s=50,
        ncols=4,
        vmax="p96",
        save="_WA_ferret_acinar.pdf",
    )
    # Plot the ferret-specific gene sets and their distribution
    sc.pl.embedding(
        adata_ferret_acinar,
        basis="umap",
        color=[
            "ARID3A",
            "CELA1",
            "CPA4",
            "EGFR",
            "EGR1",
            "FOXK2",
            "GATA6",
            "KLF15",
            "KLF6",
            "MYC",
            "NFIB",
            "NR4A1",
            "PRSS3",
            "STAT3",
            "TCF7",
            "VEGFC",
            "ZBTB16",
        ],
        color_map="viridis",
        s=50,
        ncols=4,
        vmax="p96",
        save="_FS_ferret_acinar.pdf",
    )
    # Plot the pseudotime gradient next to the fate probability
    sc.pl.embedding(
        adata_ferret_acinar,
        basis="umap",
        color=["pseudotime"] + ["fate_probabilities_acinar_maturing"],
        color_map="viridis",
        s=50,
        ncols=2,
        vmax="p96",
        save="_pseudotime_ferret_acinar.pdf",
    )


def ferret_ductal():
    """
    Runs CellRank trajectory analysis on ferret ductal pancreatic cell data.

    Loads pseudotime, lineage weights, and AnnData from disk, computes a
    PseudotimeKernel transition matrix, fits a GPCCA estimator with 3 terminal
    states (Ductal Immature_1, Ductal Proliferating, Ductal Immature_2), and
    compute fate probabilities and lineage driver genes for each terminal state.
    Saves all plots and driver gene tables to disk.

    Reads:
        pseudotime_objects/ferret_ductal_s9e0_pseudotime.csv
        pseudotime_objects/ferret_ductal_s9e0_lineage_weights.csv
        pseudotime_objects/ferret_ductal.h5ad

    Writes:
        figures/projection_ferret_ductal.svg
        Several .csv file that detail driver genes for each terminal state
        Several .pdf files that highlight different parts of the fate
        determination
    """
    # Set the output save file directories
    figdir = "figures/12_cellrank_plots/ferret_ductal"
    os.makedirs(figdir, exist_ok=True)
    sc.settings.figdir = figdir
    cr.settings.figdir = figdir

    # Load in pseudotime and weights for this specific ductal lineage
    ferret_ductal_pseudotime = pd.read_csv(
        "pseudotime_objects/ferret_ductal_s9e0_pseudotime.csv", index_col=0
    )
    ferret_ductal_weights = pd.read_csv(
        "pseudotime_objects/ferret_ductal_s9e0_lineage_weights.csv", index_col=0
    )

    # Load in the h5ad file from the ferret Seurat object
    adata_ferret_ductal = sc.read_h5ad("pseudotime_objects/ferret_ductal.h5ad")

    # Specify that we are looking at lineage 3 for this slingshot pseudotime
    adata_ferret_ductal.obs["pseudotime"] = ferret_ductal_pseudotime["Lineage3"]
    adata_ferret_ductal.obs["weights"] = ferret_ductal_weights["Lineage3"]

    # Remove any NaNs, recategorize the majority celltype as a category
    adata_ferret_ductal = adata_ferret_ductal[
        ~np.isnan(adata_ferret_ductal.obs["pseudotime"])
    ]
    adata_ferret_ductal.obs["final_renamed_celltype"] = adata_ferret_ductal.obs[
        "final_renamed_celltype"
    ].astype("category")

    # Calculate the nearest neighbors
    sc.pp.neighbors(adata_ferret_ductal, n_pcs=15, n_neighbors=15)

    # Instantiate a pseudotime kernel and calculate transition matrix
    pk = cr.kernels.PseudotimeKernel(adata_ferret_ductal, time_key="pseudotime")
    pk.compute_transition_matrix()

    # Visualize transition matrix
    pk.plot_projection(
        color="final_renamed_celltype",
        save=os.path.join(figdir, "projection_ferret_ductal.svg"),
        palette=human_palette,
        legend_loc="right",
        figsize=(8, 6),
        add_margin=0.02,
    )

    # Instantiate the GPCCA esimator and fit to the clusters using celltype
    g = cr.estimators.GPCCA(pk)
    g.fit(cluster_key="final_renamed_celltype", n_states=4)

    # Manually set the terminal states after fit, calculate fate probabilities
    g.set_terminal_states(
        states=["Ductal Immature_1", "Ductal Proliferating", "Ductal Immature_2"]
    )
    g.compute_fate_probabilities()

    # Series of plots to visualize macrostate and fate probabilities
    g.plot_macrostates(
        which="terminal",
        legend_loc="right",
        s=100,
        save=os.path.join(figdir, "ferret_ductal_macrostates_terminal.pdf"),
    )
    g.plot_macrostates(
        which="terminal",
        discrete=False,
        save=os.path.join(figdir, "ferret_ductal_macrostates_terminal_discrete.pdf"),
    )
    g.plot_fate_probabilities(
        same_plot=False,
        save=os.path.join(figdir, "ferret_ductal_fate_probabilities_separate.pdf"),
    )
    cr.pl.circular_projection(
        adata_ferret_ductal,
        keys=["final_renamed_celltype"],
        legend_loc="right",
        save="ferret_ductal_circular_projection.pdf",
    )

    # Name all celltypes to be driving factors when calculating lineage drivers
    driver_clusters = [
        "Ductal Immature",
        "Ductal Proliferating",
        "Pancreatic Progenitor",
        "Acinar Proliferating",
        "Acinar Immature",
        "Acinar Progenitor",
        "Mesenchyme Subtype 1",
    ]

    # Calculate lineage drivers with ductal immature_1 as the terminal state
    ductal_immature_1_df = g.compute_lineage_drivers(
        lineages=["Ductal Immature_1"],
        cluster_key="final_renamed_celltype",
        clusters=driver_clusters,
    )
    adata_ferret_ductal.obs["fate_probabilities_ductal_immature_1"] = (
        g.fate_probabilities["Ductal Immature_1"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_ferret_ductal,
        basis="umap",
        color=["fate_probabilities_ductal_immature_1"]
        + list(ductal_immature_1_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_ductal_immature_1_driver_genes_ferret_ductal.pdf",
    )

    # Calculate lineage drivers with ductal immature_2 as the terminal state
    ductal_immature_2_df = g.compute_lineage_drivers(
        lineages=["Ductal Immature_2"],
        cluster_key="final_renamed_celltype",
        clusters=driver_clusters,
    )
    adata_ferret_ductal.obs["fate_probabilities_ductal_immature_2"] = (
        g.fate_probabilities["Ductal Immature_2"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_ferret_ductal,
        basis="umap",
        color=["fate_probabilities_ductal_immature_2"]
        + list(ductal_immature_2_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_ductal_immature_2_driver_genes_ferret_ductal.pdf",
    )

    # Calculate lineage drivers with acinar proliferating as the terminal state
    ductal_proliferating_df = g.compute_lineage_drivers(
        lineages=["Ductal Proliferating"],
        cluster_key="final_renamed_celltype",
        clusters=driver_clusters,
    )
    adata_ferret_ductal.obs["fate_probabilities_ductal_proliferating"] = (
        g.fate_probabilities["Ductal Proliferating"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_ferret_ductal,
        basis="umap",
        color=["fate_probabilities_ductal_proliferating"]
        + list(ductal_proliferating_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_ductal_proliferating_driver_genes_ferret_ductal.pdf",
    )
    ductal_proliferating_df.to_csv(
        "csv_files/ductal_proliferating_driver_genes_ferret.csv", index=True
    )

    # Plot the well-associated gene sets and their distribution
    sc.pl.embedding(
        adata_ferret_ductal,
        basis="umap",
        color=[
            "ANXA4",
            "CFTR",
            "EHF",
            "FGFR2",
            "FGFR3",
            "GLIS3",
            "HNF1B",
            "NOTCH2",
            "SLC4A4",
            "WLS",
        ],
        color_map="viridis",
        s=50,
        ncols=4,
        vmax="p96",
        save="_WA_ferret_ductal.pdf",
    )
    # Plot the ferret-specific gene sets and their distribution
    sc.pl.embedding(
        adata_ferret_ductal,
        basis="umap",
        color=[
            "ARID4A",
            "ATF3",
            "EGFR",
            "ITGB4",
            "NFIB",
            "NFKB1",
            "STAT3",
            "TBP",
            "TSPAN8",
            "VEGFA",
            "VEGFC",
            "ZBTB4",
        ],
        color_map="viridis",
        s=50,
        ncols=4,
        vmax="p96",
        save="_FS_ferret_ductal.pdf",
    )
    # Plot the pseudotime gradient next to the fate probability
    sc.pl.embedding(
        adata_ferret_ductal,
        basis="umap",
        color=["pseudotime"] + ["fate_probabilities_ductal_proliferating"],
        color_map="viridis",
        s=50,
        ncols=2,
        vmax="p96",
        save="_pseudotime_ferret_ductal.pdf",
    )


def human_acinar():
    """
    Runs CellRank trajectory analysis on human acinar pancreatic cell data.

    Loads pseudotime, lineage weights, and AnnData from disk, computes a
    PseudotimeKernel transition matrix, fits a GPCCA estimator with 4 terminal
    states (Acinar Immature, Acinar Maturing, Pancreatic Progenitor_1,
    Pancreatic Progenitor_2), and computes fate probabilities and
    lineage driver genes for each terminal state. Saves all plots and driver
    gene tables to disk.

    Reads:
        pseudotime_objects/human_acinar_s10e2_pseudotime.csv
        pseudotime_objects/human_acinar_s10e2_lineage_weights.csv
        pseudotime_objects/human_acinar.h5ad

    Writes:
        figures/projection_human_acinar.svg
        Several .csv file that detail driver genes for each terminal state
        Several .pdf files that highlight different parts of the fate
        determination
    """

    # Set the output save file directories
    figdir = "figures/12_cellrank_plots/human_acinar"
    os.makedirs(figdir, exist_ok=True)
    sc.settings.figdir = figdir
    cr.settings.figdir = figdir

    # Load in pseudotime and weights for this specific acinar lineage
    human_acinar_pseudotime = pd.read_csv(
        "pseudotime_objects/human_acinar_s10e2_pseudotime.csv", index_col=0
    )
    human_acinar_weights = pd.read_csv(
        "pseudotime_objects/human_acinar_s10e2_lineage_weights.csv", index_col=0
    )

    # Restrict the pseudotime to match the degree of differentiation in ferret
    pseudotime_threshold = 80
    human_acinar_pseudotime[human_acinar_pseudotime > pseudotime_threshold] = np.nan
    human_acinar_weights[human_acinar_pseudotime.isna()] = 0

    # Load in the h5ad file from the ferret Seurat object
    adata_human_acinar = sc.read_h5ad("pseudotime_objects/human_acinar.h5ad")

    # Specify that we are looking at lineage 1 for this slingshot pseudotime
    adata_human_acinar.obs["pseudotime"] = human_acinar_pseudotime["Lineage1"]
    adata_human_acinar.obs["weights"] = human_acinar_weights["Lineage1"]

    # Remove any NaNs, recategorize the majority celltype as a category
    adata_human_acinar = adata_human_acinar[
        ~np.isnan(adata_human_acinar.obs["pseudotime"])
    ]
    adata_human_acinar.obs["Celltype"] = adata_human_acinar.obs["Celltype"].astype(
        "category"
    )

    # Calculate the nearest neighbors
    sc.pp.neighbors(adata_human_acinar, n_pcs=15, n_neighbors=15)

    # Instantiate a pseudotime kernel and calculate transition matrix
    pk = cr.kernels.PseudotimeKernel(adata_human_acinar, time_key="pseudotime")
    pk.compute_transition_matrix()

    # Visualize transition matrix
    pk.plot_projection(
        color="Celltype",
        save=os.path.join(figdir, "projection_human_acinar.svg"),
        palette=human_palette,
        legend_loc="right",
        figsize=(8, 6),
        add_margin=0.02,
    )

    # Instantiate the GPCCA esimator and fit to the clusters using celltype
    g = cr.estimators.GPCCA(pk)
    g.fit(cluster_key="Celltype", n_states=4)

    # Manually set the terminal states after fit, calculate fate probabilities
    g.set_terminal_states(
        states=[
            "Acinar Immature",
            "Acinar Maturing",
            "Pancreatic Progenitor_1",
            "Pancreatic Progenitor_2",
        ]
    )
    g.compute_fate_probabilities()

    # Series of plots to visualize macrostate and fate probabilities
    g.plot_macrostates(
        which="terminal",
        legend_loc="right",
        s=100,
        save=os.path.join(figdir, "human_acinar_macrostates_terminal.pdf"),
    )
    g.plot_macrostates(
        which="terminal",
        discrete=False,
        save=os.path.join(figdir, "human_acinar_macrostates_terminal_discrete.pdf"),
    )
    g.plot_fate_probabilities(
        same_plot=False,
        save=os.path.join(figdir, "human_acinar_fate_probabilities_separate.pdf"),
    )
    cr.pl.circular_projection(
        adata_human_acinar,
        keys=["Celltype"],
        legend_loc="right",
        save="human_acinar_circular_projection.pdf",
    )

    # Name all celltypes to be driving factors when calculating lineage drivers
    driver_clusters = [
        "Acinar Immature",
        "Acinar Maturing",
        "Acinar Progenitor",
        "Acinar Proliferating",
        "Acinar Secretory",
        "Pancreatic Progenitor",
    ]

    # Calculate lineage drivers with pancreatic progenitor_1 as terminal state
    panc_prog_1_df = g.compute_lineage_drivers(
        lineages=["Pancreatic Progenitor_1"],
        cluster_key="Celltype",
        clusters=driver_clusters,
    )
    adata_human_acinar.obs["fate_probabilities_pancreatic_progenitor_1"] = (
        g.fate_probabilities["Pancreatic Progenitor_1"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_human_acinar,
        basis="umap",
        color=["fate_probabilities_pancreatic_progenitor_1"]
        + list(panc_prog_1_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_pancreatic_progenitor_1_driver_genes_human_acinar.pdf",
    )

    # Calculate lineage drivers with pancreatic progenitor_2 as terminal state
    panc_prog_2_df = g.compute_lineage_drivers(
        lineages=["Pancreatic Progenitor_2"],
        cluster_key="Celltype",
        clusters=driver_clusters,
    )
    adata_human_acinar.obs["fate_probabilities_pancreatic_progenitor_2"] = (
        g.fate_probabilities["Pancreatic Progenitor_2"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_human_acinar,
        basis="umap",
        color=["fate_probabilities_pancreatic_progenitor_2"]
        + list(panc_prog_2_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_pancreatic_progenitor_2_driver_genes_human_acinar.pdf",
    )

    # Calculate lineage drivers with acinar maturing as the terminal state
    acinar_maturing_df = g.compute_lineage_drivers(
        lineages=["Acinar Maturing"], cluster_key="Celltype", clusters=driver_clusters
    )
    adata_human_acinar.obs["fate_probabilities_acinar_maturing"] = g.fate_probabilities[
        "Acinar Maturing"
    ].X.flatten()
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_human_acinar,
        basis="umap",
        color=["fate_probabilities_acinar_maturing"]
        + list(acinar_maturing_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_acinar_maturing_driver_genes_human_acinar.pdf",
    )
    acinar_maturing_df.to_csv(
        "csv_files/acinar_maturing_driver_genes_human.csv", index=True
    )

    # Calculate lineage drivers with acinar immature as the terminal state
    acinar_immature_df = g.compute_lineage_drivers(
        lineages=["Acinar Immature"], cluster_key="Celltype", clusters=driver_clusters
    )
    adata_human_acinar.obs["fate_probabilities_acinar_immature"] = g.fate_probabilities[
        "Acinar Immature"
    ].X.flatten()
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_human_acinar,
        basis="umap",
        color=["fate_probabilities_acinar_immature"]
        + list(acinar_immature_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_acinar_immature_driver_genes_human_acinar.pdf",
    )

    # Plot the well-associated gene sets and their distribution
    sc.pl.embedding(
        adata_human_acinar,
        basis="umap",
        color=[
            "ATF3",
            "BTG2",
            "CEL",
            "CELA3A",
            "CPA1",
            "CPA2",
            "GATA4",
            "NFKBIZ",
            "PTPN3",
            "RBPJL",
            "SERPINI2",
            "XBP1",
        ],
        color_map="viridis",
        s=50,
        ncols=4,
        vmax="p96",
        save="_WA_human_acinar.pdf",
    )
    # Plot the ferret-specific gene sets and their distribution
    sc.pl.embedding(
        adata_human_acinar,
        basis="umap",
        color=[
            "ARID3A",
            "CELA1",
            "CPA4",
            "EGFR",
            "EGR1",
            "FOXK2",
            "GATA6",
            "KLF15",
            "KLF6",
            "MYC",
            "NFIB",
            "NR4A1",
            "PRSS3",
            "STAT3",
            "TCF7",
            "VEGFC",
            "ZBTB16",
        ],
        color_map="viridis",
        s=50,
        ncols=4,
        vmax="p96",
        save="_FS_human_acinar.pdf",
    )
    # Plot the pseudotime gradient next to the fate probability
    sc.pl.embedding(
        adata_human_acinar,
        basis="umap",
        color=["pseudotime"] + ["fate_probabilities_acinar_maturing"],
        color_map="viridis",
        s=50,
        ncols=2,
        vmax="p96",
        save="_pseudotime_human_acinar.pdf",
    )


def human_ductal():
    """
    Runs CellRank trajectory analysis on human ductal pancreatic cell data.

    Loads pseudotime, lineage weights, and AnnData from disk, computes a
    PseudotimeKernel transition matrix, fits a GPCCA estimator with 4 terminal
    states (Pancreatic Progenitor_1, Pancreatic Progenitor_2, Ductal Immature_1,
    Ductal Immature_2), and computes fate probabilities and lineage driver genes
    for each terminal state. Saves all plots and driver gene tables to disk.

    Reads:
        pseudotime_objects/human_ductal_s10e19_pseudotime.csv
        pseudotime_objects/human_ductal_s10e19_lineage_weights.csv
        pseudotime_objects/human_ductal.h5ad

    Writes:
        figures/projection_human_ductal.svg
        Several .csv file that detail driver genes for each terminal state
        Several .pdf files that highlight different parts of the fate
        determination
    """

    # Set the output save file directories
    figdir = "figures/12_cellrank_plots/human_ductal"
    os.makedirs(figdir, exist_ok=True)
    sc.settings.figdir = figdir
    cr.settings.figdir = figdir

    # Load in pseudotime and weights for this specific ductal lineage
    human_ductal_pseudotime = pd.read_csv(
        "pseudotime_objects/human_ductal_s10e19_pseudotime.csv", index_col=0
    )
    human_ductal_weights = pd.read_csv(
        "pseudotime_objects/human_ductal_s10e19_lineage_weights.csv", index_col=0
    )

    # Load in the h5ad file from the ferret Seurat object
    adata_human_ductal = sc.read_h5ad("pseudotime_objects/human_ductal.h5ad")

    # Restrict the pseudotime to match the degree of differentiation in ferret
    pseudotime_threshold = 70
    human_ductal_pseudotime[human_ductal_pseudotime > pseudotime_threshold] = np.nan
    human_ductal_weights[human_ductal_weights.isna()] = 0

    # Specify that we are looking at lineage 2 for this slingshot pseudotime
    adata_human_ductal.obs["pseudotime"] = human_ductal_pseudotime["Lineage2"]
    adata_human_ductal.obs["weights"] = human_ductal_weights["Lineage2"]

    # Remove any NaNs, recategorize the majority celltype as a category
    adata_human_ductal = adata_human_ductal[
        ~np.isnan(adata_human_ductal.obs["pseudotime"])
    ]
    adata_human_ductal.obs["Celltype"] = adata_human_ductal.obs["Celltype"].astype(
        "category"
    )

    # Calculate the nearest neighbors
    sc.pp.neighbors(adata_human_ductal, n_pcs=15, n_neighbors=15)

    # Instantiate a pseudotime kernel and calculate transition matrix
    pk = cr.kernels.PseudotimeKernel(adata_human_ductal, time_key="pseudotime")
    pk.compute_transition_matrix()

    # Visualize transition matrix
    pk.plot_projection(
        color="Celltype",
        save=os.path.join(figdir, "projection_human_ductal.svg"),
        palette=human_palette,
        legend_loc="right",
        figsize=(8, 6),
        add_margin=0.02,
    )

    # Instantiate the GPCCA esimator and fit to the clusters using celltype
    g = cr.estimators.GPCCA(pk)
    g.fit(cluster_key="Celltype", n_states=4)

    # Manually set the terminal states after fit, calculate fate probabilities
    g.set_terminal_states(
        states=[
            "Pancreatic Progenitor_1",
            "Pancreatic Progenitor_2",
            "Ductal Immature_1",
            "Ductal Immature_2",
        ]
    )
    g.compute_fate_probabilities()

    # Series of plots to visualize macrostate and fate probabilities
    g.plot_macrostates(
        which="terminal",
        legend_loc="right",
        s=100,
        save=os.path.join(figdir, "human_ductal_macrostates_terminal.pdf"),
    )
    g.plot_macrostates(
        which="terminal",
        discrete=False,
        save=os.path.join(figdir, "human_ductal_macrostates_terminal_discrete.pdf"),
    )
    g.plot_fate_probabilities(
        same_plot=False,
        save=os.path.join(figdir, "human_ductal_fate_probabilities_separate.pdf"),
    )
    cr.pl.circular_projection(
        adata_human_ductal,
        keys=["Celltype"],
        legend_loc="right",
        save="human_ductal_circular_projection.pdf",
    )

    # Name all celltypes to be driving factors when calculating lineage drivers
    driver_clusters = [
        "Ductal Immature",
        "Ductal Maturing",
        "Ductal Proliferating",
        "Ducto-Acinar",
        "Pancreatic Progenitor",
    ]

    # Calculate lineage drivers with ductal immature_1 as the terminal state
    ducto_immature_1_df = g.compute_lineage_drivers(
        lineages=["Ductal Immature_1"], cluster_key="Celltype", clusters=driver_clusters
    )
    adata_human_ductal.obs["fate_probabilities_ductal_immature_1"] = (
        g.fate_probabilities["Ductal Immature_1"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_human_ductal,
        basis="umap",
        color=["fate_probabilities_ductal_immature_1"]
        + list(ducto_immature_1_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_ductal_immature_1_driver_genes_human_ductal.pdf",
    )
    ducto_immature_1_df.to_csv(
        "csv_files/ductal_immature_1_driver_genes_human.csv", index=True
    )

    # Calculate lineage drivers with ductal immature_2 as the terminal state
    ducto_immature_2_df = g.compute_lineage_drivers(
        lineages=["Ductal Immature_2"], cluster_key="Celltype", clusters=driver_clusters
    )
    adata_human_ductal.obs["fate_probabilities_ductal_immature_2"] = (
        g.fate_probabilities["Ductal Immature_2"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_human_ductal,
        basis="umap",
        color=["fate_probabilities_ductal_immature_2"]
        + list(ducto_immature_2_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_ductal_immature_2_driver_genes_human_ductal.pdf",
    )

    # Calculate lineage drivers with pancreatic progenitor_1 as terminal state
    pancreatic_prog_1_df = g.compute_lineage_drivers(
        lineages=["Pancreatic Progenitor_1"],
        cluster_key="Celltype",
        clusters=driver_clusters,
    )
    adata_human_ductal.obs["fate_probabilities_pancreatic_progenitor_1"] = (
        g.fate_probabilities["Pancreatic Progenitor_1"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_human_ductal,
        basis="umap",
        color=["fate_probabilities_pancreatic_progenitor_1"]
        + list(pancreatic_prog_1_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_pancreatic_progenitor_1_driver_genes_human_ductal.pdf",
    )

    # Calculate lineage drivers with pancreatic progenitor_2 as terminal state
    pancreatic_prog_2_df = g.compute_lineage_drivers(
        lineages=["Pancreatic Progenitor_2"],
        cluster_key="Celltype",
        clusters=driver_clusters,
    )
    adata_human_ductal.obs["fate_probabilities_pancreatic_progenitor_2"] = (
        g.fate_probabilities["Pancreatic Progenitor_2"].X.flatten()
    )
    # Plot the top 8 genes associated with the pathway next to the fate plot
    sc.pl.embedding(
        adata_human_ductal,
        basis="umap",
        color=["fate_probabilities_pancreatic_progenitor_2"]
        + list(pancreatic_prog_2_df.index[:8]),
        color_map="viridis",
        s=50,
        ncols=3,
        vmax="p96",
        save="_pancreatic_progenitor_2_driver_genes_human_ductal.pdf",
    )

    # Plot the well-associated gene sets and their distribution
    sc.pl.embedding(
        adata_human_ductal,
        basis="umap",
        color=[
            "ANXA4",
            "CFTR",
            "EHF",
            "FGFR2",
            "FGFR3",
            "GLIS3",
            "HNF1B",
            "NOTCH2",
            "SLC4A4",
            "WLS",
        ],
        color_map="viridis",
        s=50,
        ncols=4,
        vmax="p96",
        save="_WA_human_ductal.pdf",
    )
    # Plot the ferret-specific gene sets and their distribution
    sc.pl.embedding(
        adata_human_ductal,
        basis="umap",
        color=[
            "ARID4A",
            "ATF3",
            "EGFR",
            "ITGB4",
            "NFIB",
            "NFKB1",
            "STAT3",
            "TBP",
            "TSPAN8",
            "VEGFA",
            "VEGFC",
            "ZBTB4",
        ],
        color_map="viridis",
        s=50,
        ncols=4,
        vmax="p96",
        save="_FS_human_ductal.pdf",
    )
    # Plot the pseudotime gradient next to the fate probability
    sc.pl.embedding(
        adata_human_ductal,
        basis="umap",
        color=["pseudotime"] + ["fate_probabilities_ductal_immature_1"],
        color_map="viridis",
        s=50,
        ncols=2,
        vmax="p96",
        save="_pseudotime_human_ductal.pdf",
    )


if __name__ == "__main__":
    # Ferret lineages
    ferret_acinar()
    ferret_ductal()

    # Human lineages
    human_acinar()
    human_ductal()
