
# 🧬 SpatialDeX: A Reference-Free Spatial Transcriptomics Deconvolution Tool

**SpatialDeX** is a regression model-based tool designed for estimating cell-type proportions in multicellular-resolution spatial transcriptomics (ST) spots. It enables the characterization of tumor microenvironmental heterogeneity without requiring external single-cell RNA-seq references.

---

## 📌 Highlights

- ✅ **Reference-free** cell type deconvolution for ST data  
- 🎯 **Accurate estimation** of cell-type composition in tumor regions  
- 🔬 Supports **pan-cancer analyses** to explore tumor progression mechanisms  
- 📊 Performs comparably or better than existing reference-based and reference-free methods  
- 💡 Enables deep insights into spatial cellular architecture across tissue sections  

---

## 🚀 Getting Started

### Requirements

- R (>= 4.2.0)
- Suggested packages:
  - `Seurat`
  - `spdep`
  - `Matrix`
  - `ggplot2`
  - `data.table`

Install dependencies in R:

```r
install.packages(c("Seurat", "spdep", "Matrix", "ggplot2", "data.table"))
