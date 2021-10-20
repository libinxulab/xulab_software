## xulab_software/rnaseq/analysis

### `two_group_deseq2_analysis.R`
Performs differential expression analysis on RNAseq data between two groups of samples

#### Requires:
* `<input>.csv` (gene count data)
  * input data shoulbe be in .csv format with the headers: `Geneid, group_A_sample_1, ..., group_A_sample_N, 
  group_B_sample_1, ..., group_B_sample_N`
* `individual_genes/` (a directory to store plots of individual gene counts)

#### Produces:
* plots of mean counts vs. LFC with and without LFC shrinkage applied (a.k.a. MA plots):
  * `MA_plot_no_shrinkage.png`
  * `MA_plot_shrinkage.png`
* plots of individual gene counts (stored in individual_genes directory)":
  * `individual_genes/<gene_name>_counts.png`
* LFC and statistical significance data (.csv format, with and without LFC shrinkage applied):
  * `deseq2_results.csv`
  * `deseq2_results_shrinkage.csv`
* variance stabilizing transformed gene count data (.csv format): 
  * `deseq2_counts_vst.csv`

#### Usage:
1. Copy `two_group_deseq2_analysis.R` script to working directory with count data and `individual_genes/` directory
2. Edit `main` function in `two_group_deseq2_analysis.R` to reflect your input data: 
```R
# ! SET RUN PARAMETERS HERE !
# number of CPU cores to use
cpu_cores <- 4
# raw gene counts file (.csv format)
raw_counts_csv <- "CTRL_BC16_counts_raw.csv"
# group A label and number of samples in group
group_A <- "CTRL"
n_A <- 4
# group A label and number of samples in group
group_B <- "BC16" 
n_B <- 5
# subset of genes to use in VST (must be < total number of genes)
vst_nsub <- 10000
```
3. Run the analysis: `Rscript two_group_deseq2_analysis.R`

### `two_group_gene_counts_pca.R`
Computes PCA on gene count data (raw or VST transformed from DESeq2) and plots the projections for the first three 
principal components.

#### Requires:
* `<input>.csv` (gene count data, raw or VST transformed from DESeq2)

#### Produces:
* plots of PCA projections for the first 3 principal components:
  * `PCA_proj_PC1_PC2.png`
  * `PCA_proj_PC2_PC3.png`

#### Usage: 
1. Copy `two_group_gene_counts_pca.R` script to working directory with count data
2. Edit `main` function in `two_group_gene_counts_pca.R` to reflect your input data: 
```R
# ! SET RUN PARAMETERS HERE !
# gene counts file
gene_counts_file <- "CTRL_BC16_counts_raw.csv"
# number of samples in group A
n_A <- 4
# number of samples in group B
n_B <- 5
```
3. Run the analysis: `Rscript two_group_gene_counts_pca.R`
