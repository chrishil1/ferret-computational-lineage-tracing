library(dplyr)
library(here)
library(JASPAR2020)
library(pheatmap)
library(Seurat)
library(TFBSTools)

# Load in common packages
source(here("src/common_setup.R"))
source(here("src/04_celltype_validation_helpers.R"))
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Define the cell types that are found in both datasets
shared_celltypes <- c(
  "Alpha", "Beta", "Delta", "Pancreatic Progenitor", "Acinar Immature",
  "Acinar Maturing", "Acinar Progenitor", "Acinar Proliferating",
  "Ductal Immature", "Ductal Proliferating"
)

# Load in Seurat objects for the two species
ferret_object <- readRDS(file.path(data_dir, "04_final_ferret_w_genes.rds"))
human_object <- readRDS(file.path(data_dir, "04_final_human_w_genes.rds"))
human_subset <- subset(human_object, subset = Celltype %in% shared_celltypes)

# Predefined gene lists
alpha_genes <- c("NLGN1", "ARX", "TNR", "GCG", "SLC15A1", "PYY", "IRX2", "IRX1", "EPB41L3", "FGF14", "KCNH7", "DNAH14", "ISL1", "ACVR1C", "KLHL1", "SLC38A11", "SCG2", "KIF19", "CALB1", "MDGA2", "PIRT", "GPR119", "CDH10", "ACSBG2", "SYNPR", "FILIP1", "TOGARAM2", "OPRM1", "JAKMIP3", "ANKFN1", "ABCC11", "TMEM196", "XKR4", "SLITRK5", "GUCY2C", "KCNK16", "WNT4", "VDR", "CLSTN2", "RSPH9", "ARC", "NFE2L2", "TSPAN1", "SCIN", "RGS4", "SYT1", "RASSF6", "P2RY14", "CDH9", "SLC30A10", "UCN3", "LRRC2", "SCN3B", "LRFN5", "IGFBPL1", "NPTX2", "PRKAG3", "ANKRD55", "SLC22A2", "UBXN10", "WDR17", "ADGRG4", "SYNDIG1L", "MGLL", "OXTR", "KCNJ3", "PLEKHG6", "CDH4", "GPA33", "SLC35G2", "PCSK1N", "HRH1", "KCNA5", "KL", "GPR179", "PXDNL", "CACNA1E", "KLHL29", "F11", "KCNC1", "FAM237A", "GCNT3", "ASB18", "WDR49", "CEMIP", "PRKG2", "KCNK10", "ECT2L", "ACSL1", "CACNG8", "CTBS", "GUCA2A", "DIRAS3", "ATP6V1C2", "SLC17A6", "MGST2", "URAD", "C17orf99", "CDH13", "LIN7A")
beta_genes <- c("MYT1L", "KLHL32", "ASTN2", "SYT13", "PPP2R2C", "PLPP4", "GALNT14", "GDAP1", "SLC24A3", "GRIN2B", "IAPP", "DAB1", "SYN3", "KCNMB2", "SIDT1", "GRIK4", "PAX4", "PAPSS2", "KCNH6", "ABCA5", "FKBP5", "DYNC1I1", "RASGEF1C", "ERO1B", "SYNM", "NKX2-2", "KCNG3", "SLC7A8", "FHL1", "ADAMTSL1", "MFGE8", "VSNL1", "LMO1", "MAP3K15", "DNAH9", "STMN2", "INS", "AEBP1", "SPOCK3", "PACRG", "GRIK3", "SPTBN4", "MEGF11", "SYT16", "BEST3", "ADCYAP1R1", "GAD2", "LIMA1", "MAN1C1", "KCNT1", "CCDC85A", "ADD2", "ADAMTS20", "ZNF804A", "CD40", "MCUB", "B3GNT3", "DUSP26", "SEMA3A", "ADH1B", "EYS", "MCTP1", "GRIN1", "GPER1", "ZNF831", "GRM4", "SCRT1", "ELFN1", "TLL2", "OGDHL", "GABRB3", "MGAM", "MDFIC", "RCOR2", "RBFOX3", "RPS6KL1", "CAMK2A", "PLXNA4", "ERN2", "FAM163B", "PRUNE2", "LOXHD1", "STAC2", "FAM135B", "ECEL1", "FFAR1", "NYAP2", "SLC35D3", "REEP1", "FCER2", "GMPR", "BTBD17", "CNGB1", "NHLH1", "GALNTL6", "VSTM2B", "ZFYVE28", "AGTR1", "ASIC2", "KCNU1")
delta_genes <- c("SST", "ATP2B2", "MGAM", "ANK2", "CCDC85A", "WSCD2", "TMEM200A", "SLC4A8", "ZNF804B", "SCGN", "FYB2", "DRD2", "PDE4B", "SCG3", "SPP1", "ANO5", "NDST3", "PDE7B", "ELAVL4", "CELF3", "SYT9", "EYS", "SLC35F3", "GABRB1", "BHLHE22", "ANKEF1", "SNAP91", "BAIAP3", "THSD7B", "PTPRN", "NSG2", "SRRM4", "LRRTM2", "GABRA5", "VSTM2L", "KCNK3", "ST8SIA5", "ANKRD33B", "SKOR1", "SLC22A15", "GNAZ", "OTOF", "LZTS1", "MAPK8IP2", "PDE2A", "PRPH", "IGSF21", "CDK5R2", "AP3B2", "ABCD2", "AEBP1", "MMP24", "KCNJ5", "CDHR5", "PLCH2", "CD55", "PAK5", "TMEM179", "XKR7", "SP5", "SEZ6", "PLEKHG1", "PPFIA4", "SHANK1", "CELSR2", "RTP4", "CNIH2", "LGI3", "FFAR4", "GNGT1", "LRRC10B", "GRAMD2A", "CALY", "NAP1L2", "PDE6C", "CAMK1G", "GJD2", "NOVA2", "ABCG4", "RHAG", "SLC30A3", "ST8SIA3", "CACNA1I", "FAM181B", "HEBP2", "IGFBPL1", "OPRM1", "TMEM151B", "PARD6A", "CPNE7", "PLA2G4F", "GAD1", "AMZ1", "PLPP7", "C1QL1", "DIRAS3", "GPR27", "SOWAHA", "DNAI1", "BRD3OS")
panc_genes <- c("ZSCAN10", "NOP10", "TUBB1", "PPBP", "ID1", "RPS17", "PPIA", "SLC1A6", "TMSB10", "APELA", "NECAB2", "CA3", "INAVA", "KCTD15", "NKAIN1", "CHRNA4", "CYP11A1", "NEDD8", "NSA2", "MATK", "ATP5F1E", "RPS27", "MIF4GD", "DRAXIN", "FAM180B", "ID3", "NEU1", "RANGRF", "SH2D4B", "ACER1", "RASSF2", "FEM1A", "HMGB3", "FZD2", "SAMD14", "ZNF48", "SLC34A2", "SLC22A12", "CD101", "RPA3", "PPIF", "FJX1", "CLCN1", "SLC2A6", "COA3", "CUEDC2", "MAGEF1", "CALCA", "COX6B2", "PTGER1", "ZACN", "TMEM200C", "MMP11", "C14orf39", "ALDOB", "ZCCHC3", "KRT19", "EMX1", "F10", "CTXN1", "NMB", "ABCG8", "TBKBP1", "HLA-DMA", "GRK1", "CLDN6", "NDUFA2", "FABP12", "KLHL35", "TSKS", "GHDC", "TIMM13", "ZNHIT2", "ISG20L2", "CDK5", "BEX1", "BHLHE41", "IMP3", "WSCD1", "SPNS2", "ATP4B", "ENHO", "MEIOC", "GCHFR", "LIN28A", "SELENOV", "NAA38", "HRC", "CCDC116", "CHAC2", "SHISA4", "FAM124B", "TBXT", "CRB2", "OTUB1", "FZD4", "DKKL1", "PABPC5", "BRICD5", "MARCOL")
acin_im_genes <- c("HHIP", "COL22A1", "ADAMTS5", "OCLN", "LRRTM4", "RPH3A", "FHL2", "KLHL29", "NFASC", "CAMK4", "FAM171B", "ATP8A2", "HHEX", "EYA2", "UNC5D", "EXPH5", "CPNE3", "PMP22", "PCNX4", "ELOC", "ITGA9", "CELSR3", "CFAP299", "PHKG2", "KCNK16", "THUMPD2", "CMTM6", "TRPA1", "COPS3", "NHSL2", "MFAP1", "ARPC3", "OXCT1", "OIP5", "ITGA8", "CAP1", "GNPDA2", "FER1L6", "ARMCX4", "FBXL5", "PNKD", "TCP11L1", "TPRKB", "LHFPL4", "ARSG", "ATP8B2", "DDX20", "NARF", "GRM3", "DECR1", "ABAT", "SF3A3", "ATP5PB", "ALG8", "CDH8", "FAIM", "EVI5L", "POLR3C", "ELOA", "PIP5K1C", "SPOCK1", "MFSD11", "CAMKK1", "GPD1", "RFC2", "COLQ", "SDHC", "CCDC43", "PSMD7", "CACNB1", "TIMM22", "PSMB5", "MAPKAPK5", "SDF4", "ATP6V0D1", "AK6", "ZNF25", "ARMC10", "KLHDC2", "CRIPT", "DNAJA3", "SCD", "ARFGAP1", "ZNF48", "MTERF4", "IGHMBP2", "TMEM98", "TDRD9", "SEMA3E", "BRSK2", "WDR12", "FUT8", "GAS8", "PEX3", "SNRPE", "TMEM168", "ME2", "SMIM15", "PRKRIP1", "CDC37")
acin_mat_genes <- c("SYT15", "SLC22A4", "VIPR2", "SLC1A2", "MYO7A", "COL24A1", "SREBF1", "ZBTB16", "SV2B", "AZGP1", "RETREG1", "LARGE2", "PYGL", "TMX4", "IL22RA1", "PPP1R16B", "SLC1A5", "CPB1", "PTP4A3", "IER3", "KLF15", "BBIP1", "AQP12B", "DRGX", "RNF212", "CLEC2L", "CLPS", "ADCY1", "CDC42EP4", "PDK4", "DDIT4", "RGS2", "CTRL", "MAT1A", "OSGIN1", "GATM", "RARRES1", "FKBP14", "SERP1", "MKNK1", "TXNDC16", "KLB", "BSPRY", "ALDH3B1", "STS", "FAM3B", "GCAT", "SLC17A9", "ABHD12B", "LGALS2", "CD82", "FLVCR2", "GMPPB", "NPC1L1", "KCNJ8", "XKR9", "ENPEP", "PON1", "FAS", "PLA2G10", "BICDL2", "CRELD2", "HS1BP3", "RASSF6", "BCL2L14", "REXO2", "CSF2RA", "MST1", "HDDC2", "SLC30A2", "MANF", "PTGS1", "SLC27A2", "PRPS2", "AFM", "DNAJB9", "EHD3", "F7", "PKD2L2", "EXTL1", "SUSD1", "TCIRG1", "PRR5L", "BNIP3", "EGF", "DNLZ", "MYO1G", "KCNA3", "SNUPN", "PRIMPOL", "C1RL", "IMPDH1", "WDR12", "CDC42EP2", "GPR63", "ZNF473", "GFI1", "BOK", "CISH", "ARHGDIG")
acin_prog_genes <- c("TBX10", "SPOCD1", "CHRNA4", "SMIM28", "CCER2", "PKIB", "SDCBP2", "DEPP1", "OPCML", "ASCL1", "DLL3", "BTK", "LRP1B", "DIAPH2", "TUBA8", "ALDOB", "CDK6", "CDH17", "LIPI", "TRDN", "ITGA9", "FRMPD1", "PSENEN", "TMEM100", "RNF113A", "FAM162B", "NCF2", "TMEM238L", "VENTX", "MAGI1", "COL22A1", "PCDH7", "MAP3K1", "KIF21A", "VSTM2A", "FBXL15", "ZDHHC14", "DIAPH1", "MAP3K13", "NFIB", "FLNB", "BDNF", "RBMS1", "SIK2", "NAV3", "SPON2", "TP73", "SMAD3", "CFTR", "DPYSL5", "CHCHD3", "PADI2", "APOBEC1", "DDX1", "QSER1", "CPA1", "VEPH1", "SH3RF1", "SMIM3", "HNRNPH1", "LMBR1", "CRIM1", "PAQR6", "GABRG2", "SLC2A13", "LAMC1", "CTDSPL", "DCDC2", "ECPAS", "GRIA3", "SEC22A", "CHD3", "SOX11", "GJC2", "NPY1R", "INO80", "EPB41L2", "EYA3", "DNMT1", "GLIS3", "CACHD1", "ATP1B1", "DENND1A", "NOTCH2", "SLC4A4", "ADGRB1", "EEA1", "NTRK2", "MAPK1", "SLC8A1", "MKRN3", "AKAP14", "ARHGAP9", "SLC16A10", "PKP4", "PUS7", "PTPRT")
acin_prol_genes <- c("CPA2", "CPA1", "CEP164", "GP2", "FANCA", "SHCBP1", "LRP1B", "ALMS1", "IGF2BP1", "DYM", "LGR4", "STXBP5", "DHX15", "RBM25", "TRDN", "PAWR", "GATA6", "SNTB1", "MAPK1", "CERS6", "PCDH7", "ARHGAP12", "SLIT2", "FBXW11", "INTS8", "PTN", "PNN", "PRKD3", "FMNL2", "NAV3", "RAB10", "LRRC28", "NR6A1", "NTN1", "PAK3", "ADGRF5", "CNOT2", "DCBLD1", "KIAA1958", "SENP6", "HERC4", "MTF2", "NRXN1", "PAPOLA", "SOS1", "DCP2", "KLF7", "DHX36", "REEP3", "LAMA1", "GTF2E2", "CUL2", "NIPA1", "PPP2R5E", "MYO16", "FAT4", "BCL11A", "SCAI", "SLAIN2", "PEG10", "RARB", "P3H2", "HFM1", "TNRC6B", "SIM1", "PDE1C", "DST", "CCNI", "HECTD4", "GALNT13", "STK32B", "PTBP3", "UST", "CDH11", "TET1", "CHD7", "BMPR1B", "DDX17", "ST6GAL2", "RFX7", "GPC6", "RCAN2", "FLRT2", "COL4A5", "BTBD3", "KLHDC10", "SORCS2", "FBXL7", "ZNF521", "NTRK2", "SENP7", "MAP3K1", "KIDINS220", "KHDRBS1", "DZIP3", "UNC13C", "FBXO42", "CACHD1", "DPP6", "ADAMTS6")
duct_im_genes <- c("PKHD1", "DSCAM", "DCDC2", "FGFR2", "WWTR1", "RASGRF2", "BICC1", "LYPD6B", "RBPMS", "MPC1", "CLU", "CDH6", "MID1", "MAP3K13", "EPB41L4A", "LMCD1", "MAML2", "SHANK2", "NFIB", "FSTL4", "ANXA4", "CADM2", "SHROOM3", "CREB5", "KCNQ1", "TNIK", "NCK2", "PTPRK", "CTNND1", "NPNT", "LAMC1", "MTCL1", "MAN1A1", "MAGI2", "SH3RF1", "PTPN14", "YAP1", "ANXA3", "PLCB4", "FAT3", "FRMPD4", "COBL", "SLC16A7", "TCF7L2", "CDH1", "ACTN4", "AHNAK", "RBMS1", "SORBS1", "KLF5", "F11R", "APLP2", "ZNF385B", "FLNB", "CNN3", "SMAD3", "NRP1", "SDK2", "FBXO32", "PKP4", "RNF24", "TNFRSF21", "EGR1", "EPB41L1", "FMN2", "PFKFB3", "ATF3", "STAT3", "SEMA5A", "EZR", "RNF19A", "HOMER1", "ILDR2", "NRG3", "SPTAN1", "SMOC1", "SHB", "TIMP3", "RYR2", "PDGFD", "SNTG1", "LIPC", "TENM4", "ZSWIM6", "MYO1E", "FRMD4B", "PAH", "PPM1H", "SLC2A13", "SRGAP1", "PARVA", "NECTIN2", "SIK2", "XYLT1", "FAM91A1", "TBC1D9", "GRIP1", "MOB3B", "VCL", "TJP1")
duct_prol_genes <- c("TOP2A", "CENPF", "KIF15", "TPX2", "ECT2", "BUB1B", "MKI67", "DLGAP5", "KIF14", "HMMR", "RACGAP1", "DEPDC1", "ARHGAP11B", "PRC1", "CFTR", "INCENP", "KIFC1", "DNA2", "FOXM1", "E2F8", "HJURP", "FANCM", "BIRC5", "SDK1", "SCTR", "CDC20", "FANCD2", "TTF2", "ARL6IP1", "AURKA", "CCNB3", "PDE3A", "SAMD12", "ADAMTS16", "PPARGC1A", "UBE2C", "LTF", "RCC2", "PLK1", "NUP205", "GMNN", "POLH", "NDE1", "NEMP1", "SALL1", "FBN1", "GDNF", "HS3ST4", "NUAK1", "CPNE4", "SSRP1", "RYR3", "DTX1", "PTPDC1", "NRTN", "SLCO1C1", "SCARA3", "CFP", "FGF1", "FAM110A", "ILF2", "ARHGEF37", "ST6GALNAC3", "GANAB", "RBPJL", "ETV4", "TRIM36", "SKP2", "SMU1", "PPP2R5A", "IGFBP4", "SLC25A33", "GRB10", "PPP2R3B", "SRSF2", "KCNC2", "CDKN2C", "SRSF1", "LMF2", "ABCC3", "EIF5A", "DAZAP1", "MCM4", "KHSRP", "MED1", "MEOX2", "NUDCD3", "DSTN", "CTDSP2", "JAG1", "KANK2", "ZCCHC4", "KLHDC4", "CALR", "KCNT2", "YBX3", "FREM2", "ARPP19", "METTL5", "FOXP3")

# Combine all of the gene lists into a single object
combined_gene_list <- c(
  alpha_genes, beta_genes, delta_genes, panc_genes,
  acin_im_genes, acin_mat_genes, acin_prog_genes,
  acin_prol_genes, duct_im_genes, duct_prol_genes
)

# Calculate the initial ferret heatmap including all genes extracted
initial_ferret_heatmap <- get_validation_heatmap(
  ferret_object,
  custom_order = shared_celltypes,
  idents = "majority_celltype",
  genes_of_interest = combined_gene_list,
  cluster_rows = FALSE,
  show_rownames = FALSE
)

# Calculate the initial human heatmap including all genes extracted
initial_human_heatmap <- get_validation_heatmap(
  human_subset,
  custom_order = shared_celltypes,
  idents = "Celltype",
  genes_of_interest = combined_gene_list,
  cluster_rows = FALSE,
  show_rownames = FALSE
)

# Get expression data from the first heatmap expression
expr_ferret <- initial_ferret_heatmap$expression_data
expr_human <- initial_human_heatmap$expression_data

genes_per_cluster <- list()
all_genes_ordered <- c()

# For each celltype, only extract the genes that had high expression z-scores
# (> 2) in both species datasets
for (cluster in shared_celltypes) {
  if (cluster %in% colnames(expr_ferret) && cluster %in% colnames(expr_human)) {
    genes_ferret_high <- rownames(expr_ferret)[expr_ferret[, cluster] > 2]
    genes_human_high <- rownames(expr_human)[expr_human[, cluster] > 2]
    common_genes <- intersect(
      intersect(genes_ferret_high, genes_human_high),
      combined_gene_list
    )
    genes_per_cluster[[cluster]] <- common_genes
  }
}

# Extract all of these genes per cluster into a final list
for (cluster in shared_celltypes) {
  if (cluster %in% names(genes_per_cluster)) {
    all_genes_ordered <- c(all_genes_ordered, genes_per_cluster[[cluster]])
  }
}

# Calculate a 2nd iteration of the ferret heatmap using the identified genes
final_ferret_heatmap <- get_validation_heatmap(
  ferret_object,
  custom_order = shared_celltypes,
  idents = "majority_celltype",
  genes_of_interest = all_genes_ordered,
  cluster_rows = FALSE,
  show_rownames = FALSE
)
pdf(file.path(image_dir, "04_final_ferret_heatmap_old_method.pdf"), height = 8, width = 3)
print(final_ferret_heatmap$heatmap)
dev.off()

# Calculate a 2nd iteration of the human heatmap using the identified genes
final_human_heatmap <- get_validation_heatmap(
  human_subset,
  custom_order = shared_celltypes,
  idents = "Celltype",
  genes_of_interest = all_genes_ordered,
  cluster_rows = FALSE,
  show_rownames = FALSE
)
pdf(file.path(image_dir, "04_final_human_heatmap_old_method.pdf"), height = 8, width = 3)
print(final_human_heatmap$heatmap)
dev.off()


# Save these genes and their expresion data to output csv files
write.csv(
  final_ferret_heatmap$expression_data,
  file = file.path(csv_dir, "heatmap_ferret_expression_data_old_method.csv"),
  row.names = TRUE
)

write.csv(
  final_human_heatmap$expression_data,
  file = file.path(csv_dir, "heatmap_human_expression_data_old_method.csv"),
  row.names = TRUE
)
