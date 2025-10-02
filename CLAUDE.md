# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`stepprofiler` is an R package for temporal gene expression analysis that identifies step-wise expression patterns in RNA-seq data using DESeq2. The package specializes in detecting four key temporal expression patterns: one-step-up, one-step-down, impulsed-up, and impulsed-down gene expression patterns across experimental time points or conditions.

## Development Commands

### Package Building and Checking
```r
# Install development version
devtools::install()

# Build package
devtools::build()

# Check package (CRAN standards)
devtools::check()

# Run R CMD check locally
R CMD check .

# Generate documentation
devtools::document()
```

### Documentation and Website
```r
# Build pkgdown site
pkgdown::build_site()

# Preview documentation locally
pkgdown::preview_site()

# Update function documentation
devtools::document()
```

### CI/CD
- GitHub Actions automatically runs R-CMD-check on push/PR to master branch
- Tests on Windows, macOS, and Ubuntu with multiple R versions
- Package must pass `R CMD check` without warnings

## Code Architecture

### Core Workflow Classes
The package uses S4 classes and specialized list structures:

- **`DE_Results`**: S4 class extending DataFrame for differential expression results
- **`pairewise_compare`**: List class containing significance matrices, fold changes, and metadata
- **`stepgenes`**: List class with subclasses (`onestepup`, `onestepdown`, `impulsedup`, `impulseddown`)

### Main Analysis Pipeline
1. **Data Import**: `ezrnaseq.R` - Functions for importing count, sample, and annotation data
2. **Pairwise Comparison**: `stepProfiler.R::pairewise_compare()` - DESeq2-based differential expression across groups
3. **Pattern Detection**: `stepProfiler.R` - Functions to identify temporal patterns:
   - `one_step_up()` / `one_step_down()` - Sustained expression changes
   - `impulsed_up()` / `impulsed_down()` - Transient expression pulses
4. **Visualization**: `plot.R` - ggplot2-based plotting functions with multiple transformation options
5. **Pathway Analysis**: `stepProfiler.R::pathways()` - Integration with Reactome pathway enrichment

### Core Algorithms (`utilities_stepProfiler.R`)
- **Stability Detection**: Functions to assess expression stability across transitions
- **Impulse Detection**: Algorithms to identify transient up/down expression pulses
- **Pattern Classification**: Logic for categorizing genes into temporal expression patterns

### Data Processing Modules
- **`limma.R`**: Alternative differential expression using limma/voom pipeline
- **`summarize.R`**: Expression summarization and grouping functions
- **`get_diff.R`**: Differential expression utilities and result formatting
- **`expressed.R`**: Active expression threshold determination and detection calls

## Function Naming Conventions

- **Core workflow**: `pairewise_compare()`, `one_step_up()`, `impulsed_down()`, etc.
- **Utilities**: Prefixed with `.` for internal functions (`.subset_index()`, `.is_stable_in_nextsteps()`)
- **Data processing**: `summarizeby()`, `get_diff()`, `expressed()`
- **Import/export**: `import_rnaseq()`, `compare_two_groups()`
- **Analysis**: `pathways()` for pathway enrichment integration

## Key Dependencies

### Bioconductor Core
- **DESeq2**: Primary differential expression analysis
- **limma**: Alternative differential expression method
- **Biobase**: Core Bioconductor data structures
- **S4Vectors**: S4 class definitions and DataFrame extensions

### Analysis and Visualization
- **ggplot2** + **ggpubr**: Plotting and visualization
- **tidyverse** (dplyr, tidyr, readr): Data manipulation and import
- **reshape2**: Data reshaping for plotting

### Optional Pathway Analysis
- **clusterProfiler**: Gene ID conversion and pathway analysis
- **ReactomePA**: Reactome pathway enrichment
- **DOSE**: Disease ontology and pathway visualization

## Data Standards

### Input Data Format
- **Count matrix**: Genes as rows, samples as columns, tab-delimited
- **Sample metadata**: Must contain `group` column for experimental conditions
- **Gene annotation**: Gene IDs and symbols for downstream analysis

### Expression Thresholds
- **Active expression**: Configurable threshold (default 64 normalized counts) for presence/absence calls
- **Fold change**: Minimum fold change threshold (default 1.5-2x) for significance
- **Stability**: Maximum fold change for "stable" expression (default 1.4x)

## Temporal Pattern Definitions

### One-Step Patterns
- **One-step-up**: Low expression before step i, then sustained upregulation
- **One-step-down**: High expression before step i, then sustained downregulation

### Impulse Patterns
- **Impulsed-up**: Gradual increase to peak at step i, then gradual decrease
- **Impulsed-down**: Gradual decrease to trough at step i, then gradual increase

### Detection Logic
- Uses significance matrices (-1/0/1) combined with fold change thresholds
- Incorporates stability analysis for sustained vs. transient changes
- Applies detection calls based on active expression thresholds

## Analysis Workflow Example

```r
# 1. Import data
data_dir <- system.file("rnaseq", "multiclass", package = "stepprofiler")
raw_count <- read.delim(file.path(data_dir, "raw.count.txt"), row.names = 1)
samples <- read.delim(file.path(data_dir, "samples.txt"), row.names = 1)

# 2. Create DESeq dataset and run analysis
dds <- DESeqDataSetFromMatrix(countData = raw_count, colData = samples, design = ~ group)
dds <- DESeq(dds, parallel = TRUE)

# 3. Pairwise comparisons across groups
pwc <- pairewise_compare(dds, "group", c("BM", "prePB", "PB", "PC"),
                        active_exprs = 64, fc = 2)

# 4. Identify temporal patterns
step_up <- one_step_up(pwc, fc = 2)
step_down <- one_step_down(pwc, fc = 2)
impul_up <- impulsed_up(pwc, fc = 2)

# 5. Visualize patterns
plot(step_up, transformby = "firststep", getRegFunc = get_human_regulators)

# 6. Pathway analysis
pathways(step_up, gn_id_type = "ENSEMBL")
```

## Integration Notes

This package is part of the larger OmniCompass AI bioinformatics ecosystem and focuses on temporal RNA-seq analysis patterns. It provides foundation algorithms for identifying biologically meaningful expression transitions that can be used in downstream multi-omics integration and machine learning applications.