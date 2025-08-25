# rna-seq-deseq2-pipeline

# 🧬 Differential Gene Expression Analysis of Lung Cancer Endothelial Cells using DESeq2

This repository provides a complete RNA-Seq analysis pipeline to identify **differentially expressed genes (DEGs)** in **lung cancer endothelial cells**, using data from the **E-MTAB-8031** study. The analysis was performed using the R package **DESeq2**, along with tools for visualization and functional enrichment analysis.

---

## 📚 Background

### 🎯 Study Objective

The dataset (E-MTAB-8031) includes RNA-Seq data from **endothelial cells** sampled from **lung cancer tissue** and other anatomical locations. The aim is to detect **differential expression patterns** based on sampling sites, which may reveal **tumor-associated transcriptional changes** in endothelial cells.

### 🔬 What is Differential Expression?

Differential gene expression (DGE) analysis identifies genes whose expression levels significantly differ between experimental groups (e.g. tumor vs normal). These genes may play roles in:

- Disease progression
- Therapeutic resistance
- Cellular differentiation or transformation

---

## 🧪 Theoretical Overview

### 🧬 RNA-Seq and Count Data

RNA sequencing quantifies mRNA levels in biological samples. After mapping reads to the genome, we obtain a **count matrix**: genes as rows and samples as columns.

### 📊 Why DESeq2?

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) models count data using the **negative binomial distribution** and estimates:

- Size factors (for library normalization)
- Dispersion (biological variance)
- Fold change (log2)
- Adjusted p-values (FDR)

It uses shrinkage estimators to reduce noise and improve effect size estimates—particularly for genes with low counts.


## 🛠️ Tools & Packages

- **DESeq2** – Differential expression analysis
- **apeglm** – Shrinkage of log fold changes
- **EnhancedVolcano** – Publication-ready volcano plots
- **clusterProfiler** – Functional enrichment (GO)
- **org.Hs.eg.db** – Gene annotation
- **pheatmap** – Heatmaps and clustering
- **tidyverse** – Data wrangling and cleaning

---

## 🚀 Analysis Workflow

### 1. 📥 Download Raw Data

- Count matrix and metadata are downloaded from the following links:
  Count matrix - https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-8031/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts
  Metadata - https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8031/E-MTAB-8031.sdrf.txt

### 2. 🧹 Data Wrangling

- Metadata is filtered and formatted (factors: disease, sampling site, sex)
- Count matrix is aligned with metadata
- Low-expressed genes are filtered out

### 3. 📊 Differential Expression with DESeq2

- Experimental design: `~ sampling_site`
- DEGs are identified using thresholds:
  - Adjusted p-value (FDR) < 0.05
  - |log2FoldChange| > 1
- Shrinkage applied using `apeglm` for more stable estimates

### 4. 🖼️ Visualization

- **MA Plot**: Mean expression vs log2 fold change
- **PCA Plot**: Visualizes sample clustering by sampling site
- **Volcano Plot**: Highlights significant DEGs

### 5. 🧠 GO Enrichment Analysis

- Significant DEGs mapped to **ENTREZ IDs**
- GO analysis (Biological Process) performed via `enrichGO`
- Top terms visualized in a **bar plot**
---

## 📈 Key Results

- Significant DEGs: 126 genes

Top Enriched Biological Processes:

Regulation of body fluid levels

Mucus secretion

Body fluid secretion

Negative regulation of viral genome replication

Interpretation:

These enriched processes suggest that tumor-associated endothelial cells are involved in regulation of body fluid levels, mucus secretion, body fluid secretion and negative regulation of viral genome replication.

Notably, genes like OAS2, EGFR, and VAMP8 appear in multiple enriched pathways, indicating a central role in the tumor microenvironment.


## 🧪 How to Run the Pipeline

### 📦 Prerequisites

Install R packages:

install.packages("tidyverse")
BiocManager::install(c("DESeq2", "EnhancedVolcano", "apeglm",
                       "AnnotationDbi", "org.Hs.eg.db",
                       "clusterProfiler", "pheatmap"))

---
▶️ Run the Script

source("DESEQ2_analysis.R")
All output files (CSV tables and plots) will be saved in the /results directory.

---
👤 Author
Vasavi K 
🔗 Email: kurravasavi2@gmail.com

---
