# Using Seurat for Flow Cytometry Data Analysis
> Yes, you read that right — single-cell analysis tools like **Seurat** can be repurposed for **flow cytometry data**!
While Seurat was originally designed for single-cell RNA-seq, its powerful toolkit for dimensionality reduction, clustering, and visualization is also highly applicable to flow cytometry data — *provided the data is properly preprocessed*.

## 📋 Prerequisites
The built-in preprocessing functions in Seurat assume count-based gene expression matrices and are **not suitable** for raw flow cytometry FCS files. Therefore, **preprocessing must be performed externally** using conventional flow cytometry software (e.g., FlowJo, FCS Express, or Cytobank) before importing the data into Seurat.

### Required preprocessing steps (before export):
1. **Debris and doublet/multiplet removal**  
2. **Dead cell exclusion**  
3. **Fluorescence compensation** for spectral overlap  
4. *(Optional)* **Gating** to retain only the cell populations of interest  

### Data format for Seurat
After preprocessing, export your data as a **numeric matrix** (e.g., a CSV file) where:
- **Rows** =  individual cells 
- **Columns** = markers (e.g., CD3, CD4, CD8, etc.) 
- Values = **compensated but un-scaled** intensities readout.

> ⚠️ **Note**: This repository **does not include data loading code**. You will need to write your own script to read the preprocessed matrix (e.g., from a CSV file). The provided analysis code starts with creating Seurat object using `CreateSeuratObject()`.

## 🧪 Inspiration & Acknowledgements
This workflow is inspired by the [Seurat multimodal vignette](https://satijalab.org/seurat/articles/multimodal_vignette.html). Big thanks to the Satija Lab for developing such a flexible and powerful toolkit!

---

Happy analyzing! 🧬🔬 
