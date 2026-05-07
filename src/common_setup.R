library(here)

# Define color palettes
ferret_palette <- c(
  "#0083bc", "#00bef1", "#019538", "#ff98c0", "#efa700",
  "#ff594b", "#ffa993", "#ed087d", "#404aa2", "#c90049",
  "#4fad15", "#21e198", "#00906b", "#80421e", "#c4b5ff",
  "#9b1767", "#0d6113", "#afd098", "#b77eff", "#f646bc",
  "#c09e00", "#626699", "#55a1ff", "#91a800"
)

human_palette <- c(
  "Acinar Immature" = "#e64463",
  "Acinar Maturing" = "#62be4a",
  "Acinar Progenitor" = "#9255cb",
  "Acinar Proliferating" = "#b2b336",
  "Acinar Secretory" = "#5069d7",
  "Alpha" = "#dc8f2f",
  "B lymphocyte" = "#d579e1",
  "Beta" = "#4c8a35",
  "Delta" = "#ca42a9",
  "Dendritic Cell" = "#54c081",
  "Ductal ECM" = "#d83f86",
  "Ductal Immature" = "#52bfad",
  "Ductal Maturing" = "#cc4724",
  "Ductal Proliferating" = "#4ab4dd",
  "Ducto-Acinar" = "#dc6c4d",
  "Endothelial" = "#5982c5",
  "Endothelial Lymphatic" = "#cca15d",
  "Erythroblast" = "#9384d9",
  "Mesenchyme Proliferating" = "#727229",
  "Mesenchyme Subtype 1" = "#ae57a1",
  "Mesenchyme Subtype 2" = "#9cb26a",
  "Mesenchyme Subtype 3" = "#755298",
  "Mesenchyme Subtype 4" = "#37835d",
  "Mono-Macrophage" = "#b83849",
  "Neuron" = "#ca96d2",
  "Pancreatic Progenitor" = "#9a5e2c",
  "Phagocytic Macrophage" = "#e47da1",
  "Schwann" = "#e5907a",
  "Stromal" = "#9f4970",
  "T lymphocyte" = "#aa5659"
)

plasma <- c(
  "#001889", "#5F008C", "#91008D", "#B62485",
  "#D24E71", "#E4774D", "#F8BA43", "#F6DA23",
  "#F8FA0D"
)

blue_yellow <- c(
  "#352A86", "#343DAE", "#0262E0", "#1389D2",
  "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"
)

# Point to useful directories
data_dir <- here("rds_objects")
image_dir <- here("figures")
csv_dir <- here("csv_files")
mapping_file <- file.path(csv_dir, "species_mapping.csv")
