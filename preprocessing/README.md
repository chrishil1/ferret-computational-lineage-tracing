# Parse Biosciences Snakemake Pipeline

A Snakemake workflow for processing Parse Biosciences split-seq data. The pipeline runs Parse Biosciences Split-Pipe for initial processing, then executes a configurable set of R scripts for per-sample analysis, merged analysis, and optional subset analysis.

Written by Kristen Wells.

## Features

* Runs Parse Biosciences Split-Pipe for initial processing
* Processes samples independently with R scripts
* Supports merged analysis across selected samples
* Supports optional subset-specific downstream analysis
* Uses Snakemake for reproducible workflow management
* Runs R scripts inside a Singularity/Apptainer container

## Pipeline philosophy

This workflow is designed to balance automation with interactive analysis. Snakemake manages dependencies, parallelization, and reproducibility, while the biological analyses are performed by modular R scripts whose execution is controlled through config.yaml.

Many analysis decisions in single-cell RNA-seq—such as the number of principal components to retain or the clustering resolution—benefit from manual inspection rather than being fixed in advance. To support this, analysis parameters are stored in sample_info.tsv, allowing the workflow to pause at key decision points. After inspecting the outputs, users can update the metadata file and simply rerun Snakemake to continue the analysis from the appropriate stage.

This design makes the workflow reproducible while preserving the flexibility needed for exploratory single-cell analysis.## 

## Workflow overview

```text
FASTQs
  ↓
Parse Biosciences Split-Pipe
  ↓
Per-sample R analysis
  ↓
Merged R analysis
  ↓
Optional subset analysis
```

The exact scripts run at each stage are controlled by `config.yaml`.

## Repository layout

```text
.
├── Snakefile
├── config.yaml
├── files/
│   ├── samples.tsv
│   └── sample_info.tsv
├── docker/
│   └── single_cell.sif
├── src/
│   ├── rules/
│   │   ├── parse_pipeline.snake
│   │   └── run_rscripts.snake
│   └── scripts/
│       ├── individual_analysis/
│       │   ├── 01_Initial_processing.R
│       │   ├── 02_remove_doublets.R
│       │   ├── 03_PCA.R
│       │   ├── 04a_UMAP_find_resolution.R
│       │   ├── 04b_UMAP.R
│       │   ├── 05_name_clusters.R
│       │   └── 06_find_markers.R
│       └── integrated_analysis/
│           ├── 01_Initial_processing.R
│           ├── 02_PCA.R
│           ├── 03_batch_correction.R
│           ├── 04a_UMAP_find_resolution.R
│           ├── 04b_UMAP.R
│           ├── 05_name_clusters.R
│           ├── 06_find_markers.R
│           ├── 07_subset_cells.R
│           └── 08_subset_markers.R
└── results/
```

## Requirements

* Snakemake
* Parse Biosciences Split-Pipe
* R
* Required R packages for the analysis scripts
* Singularity or Apptainer, if running the R scripts through the provided container image

## Quick start

1. Edit `config.yaml` so that all paths are correct for your environment.
2. Prepare `files/samples.tsv` and `files/sample_info.tsv`.
3. Run the workflow:

```bash
snakemake --cores 24 --use-singularity
```

If you run on a cluster, use the appropriate Snakemake profile and resource settings for your environment.

## Configuration

All pipeline settings are defined in `config.yaml`.

### Core settings

| Setting        | Description                                                          |
| -------------- | -------------------------------------------------------------------- |
| `SAMPLE_SHEET` | Path to the sample sheet (`files/samples.tsv`)                       |
| `RESULTS`      | Output directory                                                     |
| `GENOME`       | Cell Ranger genome reference used by Split-Pipe                      |
| `CHEMISTRY`    | Parse chemistry version (`v1` or `v2`)                               |
| `MAX_JOBS`     | Maximum number of jobs per sample submitted by Cell Ranger at a time |

### R script settings

| Setting          | Description                                      |
| ---------------- | ------------------------------------------------ |
| `SCRIPT_PATH`    | Path to the directory containing R scripts       |
| `SCRIPTS_RUN`    | Scripts to run globally, if any                  |
| `SAMPLE_SCRIPTS` | Scripts to run on individual samples             |
| `MERGE_SCRIPTS`  | Scripts to run on merged samples                 |
| `SUBSET_SCRIPTS` | Scripts to run on subsets of the merged analysis |

### Metadata and container

| Setting             | Description                                 |
| ------------------- | ------------------------------------------- |
| `SAMPLE_INFO`       | Path to `files/sample_info.tsv`             |
| `RSCRIPT_CONTAINER` | Singularity image used to run the R scripts |

## Input files

### `samples.tsv`

This file defines the samples and libraries to process.

Required columns:

```text
sample_name	fastq_files	sample_wells	data_dir	library_name
```

Column descriptions:

| Column         | Description                                                                                                      |
| -------------- | ---------------------------------------------------------------------------------------------------------------- |
| `sample_name`  | Unique sample identifier                                                                                         |
| `fastq_files`  | Comma-separated list of FASTQ filenames for the sample                                                           |
| `sample_wells` | Sample well identifier(s)                                                                                        |
| `data_dir`     | Directory containing the FASTQ files; may be a single directory or a comma-separated list matching `fastq_files` |
| `library_name` | Library name(s) associated with the FASTQ files; must match the length of `fastq_files`                          |

### `sample_info.tsv`

This file stores sample-level analysis settings and controls which samples are included in the individual analysis stage.

Required columns:

```text
sample	merged	HTO	ADT	hash_ident	adt_pca	adt_umap	PCs	resolution	batch_correction
```

Column descriptions:

| Column             | Description                                                                   |
| ------------------ | ----------------------------------------------------------------------------- |
| `sample`           | Sample identifier                                                             |
| `merged`           | Whether the sample should be included in the individual sample analysis stage |
| `HTO`              | Whether hashtag oligos were used                                              |
| `ADT`              | Whether antibody-derived tags were used                                       |
| `hash_ident`       | Hashtag identity assignment                                                   |
| `adt_pca`          | Whether ADT PCA should be run                                                 |
| `adt_umap`         | Whether ADT UMAP should be run                                                |
| `PCs`              | Number of principal components to use                                         |
| `resolution`       | Clustering resolution                                                         |
| `batch_correction` | Batch correction setting used in merged analysis                              |

Only rows where `merged == False` are used for the individual sample analysis stage.

## Analysis stages

**Note** R scripts for the analysis are found in `src/scripts`. There are two directories here. One called `individual_analysis` includes scripts to process an individual sample. `intigrated_analysis` contains scripts to run analysis on multiple samples.

To run these scripts, you will need to install my custom analysis package `scAnalysisR`:

```R
install.packages("devtools")
library(devtools)
install_github("CUAnschutzBDC/scAnalysisR")
```

For many of these steps, I've included different methods of performing them (for example batch correction can be done with `harmony` or `fastMNN`) and ways of comparing the different approaches so you can make an informed decision about what method is best for your samples.

For all scripts in this analysis, please make sure you understand the processing steps so you know what are the appropriate parameters for each step. These scripts are not intended to be run blindly with no understanding of the underlying algorithims.

### 1. Initial processing with Split-Pipe

For each library, the workflow runs Parse Biosciences Split-Pipe using:

* the FASTQ files listed in `samples.tsv`
* the reference genome specified by `GENOME`
* the chemistry specified by `CHEMISTRY`

This step produces the initial count outputs used by the downstream R analysis.

### 2. Per-sample analysis

The following R scripts are run in order for each sample listed in `sample_info.tsv` with `merged == False`:

1. `01_Initial_processing.R`
2. `02_remove_doublets.R`
3. `03_PCA.R`
4. `04a_UMAP_find_resolution.R`
5. `04b_UMAP.R`
6. `05_name_clusters.R`
7. `06_find_markers.R`

Scripts are run sequentially within each sample. Different samples may run in parallel.

### 3. Merged analysis

Merged analysis is controlled by `MERGE_SCRIPTS` in `config.yaml`. In the current configuration, the scripts run in order are:

1. `01_Initial_processing.R`
2. `02_PCA.R`
3. `03_batch_correction.R`
4. `04a_UMAP_find_resolution.R`
5. `04b_UMAP.R`
6. `05_name_clusters.R`
7. `06_find_markers.R`
8. `07_subset_cells.R`

### 4. Subset analysis

Optional subset analyses can be defined under `SUBSET_SCRIPTS`. In the current configuration, the example subset is:

* `pancreas`

The subset workflow runs:

1. `02_PCA.R`
2. `03_batch_correction.R`
3. `04a_UMAP_find_resolution.R`
4. `04b_UMAP.R`
5. `05_name_clusters.R`
6. `08_subset_markers.R`

## Outputs

All results are written to the directory specified by `RESULTS`.

Expected outputs include:

* Split-Pipe outputs for each library
* Completion marker files for each library
* Per-sample analysis outputs
* Merged analysis outputs
* Subset analysis outputs
* Script completion marker files

The workflow creates the output directory if it does not already exist.

## Customization

The workflow can be customized by editing `config.yaml`:

* Add or remove scripts in `SAMPLE_SCRIPTS`
* Add or remove scripts in `MERGE_SCRIPTS`
* Define new subset analyses under `SUBSET_SCRIPTS`
* Update `SAMPLE_INFO` to match the metadata required for each sample

## Notes

* Scripts listed in the config run in the order provided.
* Different samples may be processed in parallel.
* Merged and subset analyses are fully controlled through the config file.
* The pipeline uses strict shell settings so that errors fail fast.

## Example input files

See the example files in `files/` for the expected input format:

* `files/samples.tsv`
* `files/sample_info.tsv`
