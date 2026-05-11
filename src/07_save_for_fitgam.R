options(Seurat.object.assay.version = 'v3')
library(Seurat)
# library(dplyr)
# library(cowplot)
# library(ggbeeswarm)
# library(ggrastr)
library(here)
# library(magrittr)
# library(scAnalysisR)
library(slingshot)
# library(tradeSeq)
# library(scCustomize)

source(here("src/common_setup.R"))

#################
# Ferret saving #
#################
# Load in ferret object
ferret_subset <- readRDS(file.path(data_dir, "05_final_ferret.rds"))

# Subset to ductal clusters only
ductal_clusters <- c(9, 0, 5, 12, 17)
cluster_subset <- subset(ferret_subset,
                         subset = seurat_clusters %in% ductal_clusters
)

# Convert to single cell experiment
sce <- as.SingleCellExperiment(cluster_subset)

# Extract only the top 15 PCs for slingshot object
pca_subset <- reducedDim(sce, "PCA")[, 1:15]

# Calculate the slingshot trajectory for our chosen start/end clusters
slingshot_sce <- slingshot(
  sce,
  reducedDim = pca_subset,
  clusterLabels = sce$seurat_clusters,
  start.clus = 9,
  end.clus = 0
)

# Save the slingshot SCE object and the subset of cells involved
saveRDS(slingshot_sce, file.path(fitgam_dir, "ferret_ductal_S9E0_sce.rds"))
saveRDS(cluster_subset, file.path(fitgam_dir, "ferret_ductal_S9E0.rds"))

# Subset to acinar clusters only
acinar_clusters <- c(9, 4, 14, 1, 16, 8, 11, 10)
cluster_subset <- subset(ferret_subset,
                         subset = seurat_clusters %in% acinar_clusters
)

# Convert to single cell experiment
sce <- as.SingleCellExperiment(cluster_subset)

# Extract only the top 15 PCs for slingshot object
pca_subset <- reducedDim(sce, "PCA")[, 1:15]

# Calculate the slingshot trajectory for our chosen start/end clusters
slingshot_sce <- slingshot(
  sce,
  reducedDim = pca_subset,
  clusterLabels = sce$seurat_clusters,
  start.clus = 9,
  end.clus = 10
)

# Save the slingshot SCE object and the subset of cells involved
saveRDS(slingshot_sce, file.path(fitgam_dir, "ferret_acinar_S9E10_sce.rds"))
saveRDS(cluster_subset, file.path(fitgam_dir, "ferret_acinar_S9E10.rds"))

#################
# Human saving #
#################
# Load in ferret object
human_object <- readRDS(file.path(data_dir, "04_human_remade_seurat.rds"))
Idents(human_object) <- "Leiden_cluster"

# Subset to ductal clusters only
ductal_clusters <- c(5, 8, 10, 19, 21, 27)
cluster_subset <- subset(human_object,
                         subset = Leiden_cluster %in% ductal_clusters
)

# Convert to single cell experiment
sce <- as.SingleCellExperiment(cluster_subset)

# Extract only the top 15 PCs for slingshot object
pca_subset <- reducedDim(sce, "PCA")[, 1:15]

# Calculate the slingshot trajectory for our chosen start/end clusters
slingshot_sce <- slingshot(
  sce,
  reducedDim = pca_subset,
  clusterLabels = sce$Leiden_cluster,
  start.clus = 10,
  end.clus = 19
)

# Save the slingshot SCE object and the subset of cells involved
saveRDS(slingshot_sce, file.path(fitgam_dir, "human_ductal_S10E19_sce.rds"))
saveRDS(cluster_subset, file.path(fitgam_dir, "human_ductal_S10E19.rds"))

# Subset to acinar clusters only
acinar_clusters <- c(2, 3, 4, 6, 10, 11)
cluster_subset <- subset(human_object,
                         subset = Leiden_cluster %in% acinar_clusters
)

# Convert to single cell experiment
sce <- as.SingleCellExperiment(cluster_subset)

# Extract only the top 15 PCs for slingshot object
pca_subset <- reducedDim(sce, "PCA")[, 1:15]

# Calculate the slingshot trajectory for our chosen start/end clusters
slingshot_sce <- slingshot(
  sce,
  reducedDim = pca_subset,
  clusterLabels = sce$Leiden_cluster,
  start.clus = 10,
  end.clus = 2
)

# Save the slingshot SCE object and the subset of cells involved
saveRDS(slingshot_sce, file.path(fitgam_dir, "human_acinar_S10E2_sce.rds"))
saveRDS(cluster_subset, file.path(fitgam_dir, "human_acinar_S10E2.rds"))
