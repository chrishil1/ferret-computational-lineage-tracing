library(slingshot)
library(Seurat)
library(tradeSeq)
library(here)

harmony_object <- readRDS(here("ferret_ductal_9_26.rds"))
harmony_slingshot <- readRDS(here("ferret_ductal_9_26_sce.rds"))

slingshot_cells <- cellnames(harmony_slingshot$slingshot)
harmony_object <- subset(harmony_object, cells = slingshot_cells)

harmony_object <- FindVariableFeatures(harmony_object, nfeatures = 3000)
counts <- GetAssayData(harmony_object, slot = "data")

counts <- counts[rownames(counts) %in% VariableFeatures(harmony_object),]

#set.seed(0)
#pdf(here("evaluateK_plots_ferret_ductal_9_26.pdf"))
#icMat <- evaluateK(counts = counts, sds = harmony_slingshot$slingshot, k = 3:10, 
#                   nGenes = 200, verbose = T)

#dev.off()
#saveRDS(icMat, here("icmat_ferret_ductal_9_26.rds"))

# Run fitGAM - note for future implementations (ie when accounting for batch,
# covariates can be included here.)

print("Running fitGAM!")
set.seed(0)
nknots <- 7
pseudotime <- slingPseudotime(harmony_slingshot$slingshot, na = FALSE)
cellWeights <- slingCurveWeights(harmony_slingshot$slingshot)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = nknots, verbose = TRUE)

saveRDS(sce, here("fitgam_ferret_ductal_9_26.rda"))
