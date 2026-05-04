library(here)
library(Matrix)
library(Seurat)
options(future.globals.maxSize = 8000 * 1024^2)

# Load in the raw samples
seurat_list <- list()
samples <- c(
  "GSM7221604_sample1", "GSM7221605_sample2", "GSM7221606_sample3",
  "GSM7221607_sample4", "GSM7221608_sample5", "GSM7221609_sample6"
)

# Create a Seurat object for each of the six samples individually
sample_names <- c("HFP35", "HFP29", "HFP22", "HFP19", "HFP31", "HFP32")
for (i in 1:6) {
  barcodes <- read.delim(here("raw_data", paste0(samples[i], "_barcodes.tsv")),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  features <- read.delim(here("raw_data", paste0(samples[i], "_features.tsv")),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  matrix <- readMM(here("raw_data", paste0(samples[i], "_matrix.mtx")))

  colnames(matrix) <- paste0(sample_names[i], "_", barcodes$V1)
  rownames(matrix) <- features$V1

  # Add the current Seurat object to a list
  seurat_list[[i]] <- CreateSeuratObject(counts = matrix)
}

# Merged the Seurat objects into a single object
merged_so <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])

# First round of SCTransform and Scale data
merged_so <- SCTransform(merged_so)
merged_so <- ScaleData(merged_so)

# Read in published seurat object
published_object <- readRDS(here("human-fetal-panc-allcells.rds"))

# Subset merged object by cells that are included in final published object
pub_cells <- Cells(published_object)
subset_seurat <- subset(merged_so, cells = pub_cells)

# Match metadata rows to the exact cell order in subset_seurat
matched_cells <- Cells(subset_seurat)
pub_meta <- pub_seurat@meta.data[matched_cells, , drop = FALSE]

# Direct assignment of metadata
subset_seurat@meta.data <- pub_meta

# Second round of SCTransform
subset_seurat <- SCTransform(subset_seurat)
subset_seurat <- RunPCA(subset_seurat,
  features = VariableFeatures(subset_seurat),
  npcs = 50
)

# Get projection in UMAP space
subset_seurat <- ProjectDim(subset_seurat)
subset_seurat <- RunUMAP(subset_seurat,
  dims = 1:40,
  seed.use = 1,
  n.neighbors = 50
)

# Save to output RDS file
saveRDS(subset_seurat, here("recreated_published_object.rds"))
