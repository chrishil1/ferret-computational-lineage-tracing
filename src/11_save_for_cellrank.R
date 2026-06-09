library(dplyr)
library(here)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(slingshot)
options(Seurat.object.assay.version = "v3")

source(here("src/common_setup.R"))

# Read in our final species objects
final_ferret <- readRDS(file.path(data_dir, "05_final_ferret.rds"))
final_human <- readRDS(file.path(data_dir, "05_final_human.rds"))

#######################
# Human object saving #
#######################
# Subset to ductal clusters only
ductal_clusters <- c(5, 8, 10, 19, 21, 27)
ductal_subset <- subset(final_human,
  subset = Leiden_cluster %in% ductal_clusters
)

# Convert to single cell experiment
sce <- as.SingleCellExperiment(ductal_subset)

# Extract only the top 15 PCs for slingshot object
pca_subset <- reducedDim(sce, "PCA")

# Run slingshot object with our decided start/end clusters
slingshot_sce <- slingshot(
  sce,
  reducedDim = pca_subset,
  clusterLabels = sce$Leiden_cluster,
  start.clus = 10,
  end.clus = 19
)

# Convert the Celltype entries into characters
ductal_subset$Celltype <- as.character(ductal_subset$Celltype)

# Extract the pseudotime and lineage weights from the slingshot object
ductal_human_pseudotime <- slingPseudotime(slingshot_sce)
ductal_human_lineage_weights <- slingCurveWeights(slingshot_sce)

# Write pseudotime and lineage weights to output files
write.csv(ductal_human_pseudotime,
  file = "pseudotime_objects/human_ductal_s10e19_pseudotime.csv"
)
write.csv(ductal_human_lineage_weights,
  file = "pseudotime_objects/human_ductal_s10e19_lineage_weights.csv"
)

# Save the ductal subset for the human dataset to a python compatible file
SaveH5Seurat(ductal_subset,
  filename = "pseudotime_objects/human_ductal.h5Seurat"
)
Convert("pseudotime_objects/human_ductal.h5Seurat", dest = "h5ad")

# Subset to acinar clusters only
acinar_clusters <- c(2, 3, 4, 6, 10, 11)
acinar_subset <- subset(final_human,
  subset = Leiden_cluster %in% acinar_clusters
)

# Convert to single cell experiment
sce <- as.SingleCellExperiment(acinar_subset)

# Extract only the top 15 PCs for slingshot object
pca_subset <- reducedDim(sce, "PCA")

# Run slingshot object with our decided start/end clusters
slingshot_sce <- slingshot(
  sce,
  reducedDim = pca_subset,
  clusterLabels = sce$Leiden_cluster,
  start.clus = 10,
  end.clus = 2
)

# Convert the Celltype entries into characters
acinar_subset$Celltype <- as.character(acinar_subset$Celltype)

# Extract the pseudotime and lineage weights from the slingshot object
acinar_human_pseudotime <- slingPseudotime(slingshot_sce)
acinar_human_lineage_weights <- slingCurveWeights(slingshot_sce)

# Write pseudotime and lineage weights to output files
write.csv(acinar_human_pseudotime,
  file = "pseudotime_objects/human_acinar_s10e2_pseudotime.csv"
)
write.csv(acinar_human_lineage_weights,
  file = "pseudotime_objects/human_acinar_s10e2_lineage_weights.csv"
)

# Save the acinar subset for the human dataset to a python compatible file
SaveH5Seurat(acinar_subset,
  filename = "pseudotime_objects/human_acinar.h5Seurat"
)
Convert("pseudotime_objects/human_acinar.h5Seurat", dest = "h5ad")

########################
# Ferret object saving #
########################
# Subset to ductal clusters only
ductal_clusters <- c(9, 0, 5, 12, 17)
ductal_subset <- subset(final_ferret,
  subset = seurat_clusters %in% ductal_clusters
)

# Convert to single cell experiment
sce <- as.SingleCellExperiment(ductal_subset)

# Extract only the top 15 PCs for slingshot object
pca_subset <- reducedDim(sce, "PCA")[, 1:15]

# Run slingshot object with our decided start/end clusters
slingshot_sce <- slingshot(
  sce,
  reducedDim = pca_subset,
  clusterLabels = sce$seurat_clusters,
  start.clus = 9,
  end.clus = 0
)

# Extract the pseudotime and lineage weights from the slingshot object
ductal_ferret_pseudotime <- slingPseudotime(slingshot_sce)
ductal_ferret_lineage_weights <- slingCurveWeights(slingshot_sce)

# Write pseudotime and lineage weights to output files
write.csv(ductal_ferret_pseudotime,
  file = "pseudotime_objects/ferret_ductal_s9e0_pseudotime.csv"
)
write.csv(ductal_ferret_lineage_weights,
  file = "pseudotime_objects/ferret_ductal_s9e0_lineage_weights.csv"
)

# Convert the ferret object into RNA3 based naming scheme for CellRank
ductal_subset[["RNA3"]] <- as(object = ductal_subset[["RNA"]], Class = "Assay")
DefaultAssay(ductal_subset) <- "RNA3"
ductal_subset[["RNA"]] <- NULL
ductal_subset <- RenameAssays(object = ductal_subset, RNA3 = "RNA")
ductal_subset[["RNA"]]$data <- NULL
ductal_subset[["RNA"]]$scale.data <- NULL

# Save the ductal subset for the ferret dataset to a python compatible file
SaveH5Seurat(ductal_subset,
  filename = "pseudotime_objects/ferret_ductal.h5Seurat"
)
Convert("pseudotime_objects/ferret_ductal.h5Seurat", dest = "h5ad")

# Subset to acinar clusters only
acinar_clusters <- c(9, 4, 14, 1, 16, 8, 11, 10)
acinar_subset <- subset(final_ferret,
  subset = seurat_clusters %in% acinar_clusters
)

# Convert to single cell experiment
sce <- as.SingleCellExperiment(acinar_subset)

# Extract only the top 15 PCs for slingshot object
pca_subset <- reducedDim(sce, "PCA")[, 1:15]

# Run slingshot object with our decided start/end clusters
slingshot_sce <- slingshot(
  sce,
  reducedDim = pca_subset,
  clusterLabels = sce$seurat_clusters,
  start.clus = 9,
  end.clus = 10
)

# Extract the pseudotime and lineage weights from the slingshot object
acinar_ferret_pseudotime <- slingPseudotime(slingshot_sce)
acinar_ferret_lineage_weights <- slingCurveWeights(slingshot_sce)

# Write pseudotime and lineage weights to output files
write.csv(acinar_ferret_pseudotime,
  file = "pseudotime_objects/ferret_acinar_s9e10_pseudotime.csv"
)
write.csv(acinar_ferret_lineage_weights,
  file = "pseudotime_objects/ferret_acinar_s9e10_lineage_weights.csv"
)

# Convert the ferret object into RNA3 based naming scheme for CellRank
acinar_subset[["RNA3"]] <- as(object = acinar_subset[["RNA"]], Class = "Assay")
DefaultAssay(acinar_subset) <- "RNA3"
acinar_subset[["RNA"]] <- NULL
acinar_subset <- RenameAssays(object = acinar_subset, RNA3 = "RNA")

acinar_subset[["RNA"]]$data <- NULL
acinar_subset[["RNA"]]$scale.data <- NULL

# Save the acinar subset for the ferret dataset to a python compatible file
SaveH5Seurat(acinar_subset,
  filename = "pseudotime_objects/ferret_acinar.h5Seurat"
)
Convert("pseudotime_objects/ferret_acinar.h5Seurat", dest = "h5ad")
