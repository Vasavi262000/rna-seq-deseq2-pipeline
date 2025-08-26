# ==============================================================================
# RNA-Seq differential expression analysis using DESeq2
# Project: E-MTAB-8031 - Lung cancer endothelial cells
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Load required libraries

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(EnhancedVolcano)
library(apeglm)

# ------------------------------------------------------------------------------
# Download data
# ------------------------------------------------------------------------------
# Raw counts
counts <- read.delim(
  "https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-8031/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts",
  header = TRUE, check.names = FALSE
)
head(counts)

# Metadata
metadata <- read.delim(
  "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8031/E-MTAB-8031.sdrf.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)
head(metadata)

# ------------------------------------------------------------------------------
# Wrangle metadata
# ------------------------------------------------------------------------------
clean_metadata <- metadata[, c(
  "Comment[ENA_RUN]",
  "Characteristics[disease]",
  "Characteristics[individual]",
  "Characteristics[sex]",
  "Characteristics[sampling site]",
  "Characteristics[organism part]",
  "Characteristics[cell type]"
)]
colnames(clean_metadata) <- c(
  "sample_id", "disease", "individual", "sex", "sampling_site", "organ", "cell_type"
)
rownames(clean_metadata) <- clean_metadata$sample_id
clean_metadata$disease <- factor(clean_metadata$disease)
clean_metadata$sex <- factor(clean_metadata$sex)
clean_metadata$sampling_site <- factor(clean_metadata$sampling_site)

# ------------------------------------------------------------------------------
# Wrangle count matrix
# ------------------------------------------------------------------------------
# Save gene annotation info
genes <- counts[, c("Gene ID", "Gene Name")]

# Clean counts
count_data <- counts
rownames(count_data) <- count_data$`Gene ID`
count_data <- count_data[, -c(1, 2)]

# Align metadata
clean_metadata <- clean_metadata[colnames(count_data), ]
all(colnames(count_data) == rownames(clean_metadata))

# ------------------------------------------------------------------------------
# Run DESeq2
# ------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = clean_metadata,
  design = ~ sampling_site
)

# Filter low-expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq
dds <- DESeq(dds)
res=results(dds)
res

#Shrink log2 fold changes for accurate effect size
resLFC <- lfcShrink(dds,
                    coef = resultsNames(dds)[2],
                    type = "apeglm")

#Summary & Export
summary(resLFC)

#Filter significant results
resSig=resLFC[which(resLFC$padj<0.05 & 
                      abs(resLFC$log2FoldChange)>1),]
resSig=resSig[order(resSig$padj),]
nrow(resSig)


#Annotate gene symbols
resSig$symbol=mapIds(org.Hs.eg.db,
                     keys=rownames(resSig),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multivals="first")

#Create results folder
if(!dir.exists("results"))dir.create("results")

#Save results tables
write.csv(as.data.frame(resLFC),
          "results/all_results_lfcShrink.csv")
write.csv(as.data.frame(resSig),
          "results/significant_DEGs.csv")


# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------
#MA Plot
png("results/MAplot_lfcShrink.png")
plotMA(resLFC,ylim=c(-5,5))
dev.off()

#Variance Stabilizing Transformation(VST)
vsd=vst(dds,blind=FALSE)


colnames(colData(vsd))
head(vsd$sampling_site)


#PCA Plot by sampling site
png("results/PCA_plot.png")
plotPCA(vsd,intgroup="sampling_site")
dev.off()

#Volcano plot
png("results/Volcano_plot.png",width=800,height=600)
EnhancedVolcano(resLFC,
                lab=rownames(resLFC),
                x='log2FoldChange',
                y='padj',
                title='Volcano Plot',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize=2.0,
                labSize=3.0)
dev.off()


# Prepare significant genes
sig_genes <- rownames(resSig)

# Remove version numbers from ENSEMBL IDs if present
sig_genes <- sub("\\..*$", "", sig_genes)

# Convert ENSEMBL to ENTREZ
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = sig_genes,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Remove NAs
entrez_ids <- na.omit(entrez_ids)

# Perform ORA using enrichGO
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",           
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE)

# Bar plot
png("results/GO_barplot.png", width = 1000, height = 800)
barplot(ego, showCategory = 20, title = "GO Enrichment (BP)")
dev.off()

# Save the enrichment results
write.csv(as.data.frame(ego), file = "results/GO_ORA_results.csv")




