# 🛠️ SpatialDeX: Development Repository

This repository documents the **development process** of **SpatialDeX** (Spatial Deconvolution Explorer), a regression-based method for estimating **cell type proportions** in tumor spatial transcriptomics (ST) data.

> ⚠️ This is **not a finalized R package**, but a development workspace. Scripts and functions here support model building, feature extraction, and testing of core algorithms used in SpatialDeX.

---

## 📖 Overview

SpatialDeX aims to address the challenge of **cell type deconvolution** in spatial transcriptomics data, especially for **multi-cell resolution** platforms. It incorporates:

- **CNV smoothing and segmentation**
- **Ridge regression**
- **Feature coefficient inference**
- **Cell-type identity modeling**

---



## ⚙️ Dependencies

This repository was developed using:

- R (>= 3.5.0)
- R packages:
  - `GSVA (1.50.5)`
  - `copykat` ([v1.1.0](https://github.com/navinlabcode/copykat))
  - `dplyr`
  - `devtools`

---

## For final version installation

```r
# Install devtools if needed
install.packages("devtools")

# Install SpatialDeX from GitHub
devtools::install_github("wang-lab/SpatialDeX", build_vignettes = TRUE)
