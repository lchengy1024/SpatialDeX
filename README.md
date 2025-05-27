# 🧬 SpatialDeX: A Reference-Free Spatial Transcriptomics Deconvolution Tool

**SpatialDeX** (Spatial Deconvolution Explorer) is a regression model-based method for estimating **cell type proportions** in tumor spatial transcriptomics (ST) data. It enables **reference-free** deconvolution of spatial spots, allowing exploration of the **spatial cellular organization** of tissues without requiring single-cell RNA-seq references.

---

## 📖 Background

While recent advances have enabled single-cell resolution ST, many current ST platforms remain at **multi-cellular resolution**. This makes **deconvolution** critical to infer the identity and proportion of cell types within each spot. SpatialDeX addresses this by integrating **CNV smoothing**, **ridge regression**, and **statistical filtering** to estimate major cell types in tumor tissue sections.

---

### 📦 Prerequisites

- R (>= 3.5.0)
- R packages:
  - `GSVA (1.50.5)`
  - `copykat` [v1.1.0](https://github.com/navinlabcode/copykat)
  - `dplyr`
  - `devtools`

### 📥 Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install SpatialDeX from GitHub
devtools::install_github("wang-lab/SpatialDeX", build_vignettes = TRUE)
