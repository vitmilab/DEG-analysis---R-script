# =========================================================
# Differential Expression Analysis of GSE140928 RNA-Seq Dataset
# Platform: Illumina HiSeq 2000
# Organism: Sorghum bicolor
# Analysis Method: DESeq2
# =========================================================

# =========================================================
# 1. Set Working Directory
# =========================================================
setwd("")

# =========================================================
# 2. Load Required Libraries
# =========================================================
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

# =========================================================
# 3. Import RNA-Seq Count Data
# =========================================================
counts <- read.csv("saline.csv", header = TRUE, row.names =1, sep = ",")
counts

# =========================================================
# 4. Filter Low-Expression Genes
# =========================================================

# Remove genes with low read counts
# Threshold: row sum > 50
counts1 <- counts[which(rowSums(counts) > 50),]
counts1

# =========================================================
# 5. Define Experimental Conditions
# =========================================================
condition1<-factor(c("stress","stress","stress","control","control","control","stress","stress","stress"))
coldata1 <- data.frame(row.names= colnames(counts1), condition1)

# =========================================================
# 6. Construct DESeq2 Dataset
# =========================================================
dds <- DESeqDataSetFromMatrix(countData =counts1, colData = coldata1,design = ~condition1)

# =========================================================
# 7. Differential Expression Analysis Using DESeq2
# =========================================================
dds<- DESeq(dds) 
res <- results(dds, contrast = c("condition1", "stress", "control")) 
res

# =========================================================
# 8. Extract Significant Differentially Expressed Genes
# =========================================================

# Thresholds:
# Adjusted p-value < 0.05
# |log2FoldChange| > 1
sigs <- na.omit(res) 
sigs <-sigs[sigs$padj < 0.05,]
sigs
write.csv(res, file = "GSE140928_DEGs.csv")

# Export DESeq2 results
write.csv(as.data.frame(res),
          file = "GSE140928_DEGs.csv")

# Import DEG results for visualization
de_genes <- read.csv("GSE140928_DEGs.csv")colnames(de_genes)

# =========================================================
# 9. Classify Differential Expression Status
# ========================================================= 
# Define differential expression categories
de_genes$diffexpressed <- "Not Significant"

de_genes$diffexpressed[
  de_genes$log2FoldChange > 1 &
  de_genes$padj <= 0.05
] <- "Upregulated"

de_genes$diffexpressed[
  de_genes$log2FoldChange < -1 &
  de_genes$padj <= 0.05
] <- "Downregulated"

# =========================================================
# 10. Volcano Plot Visualization
# =========================================================
# Generate volcano plot
ggplot(data = de_genes,
       aes(x = log2FoldChange,
           y = -log10(pvalue),
           color = diffexpressed)) +

  geom_point(size = 2,
             alpha = 0.7) +

  theme_minimal() +

  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +

  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +

  xlim(-5, 5) +

  scale_color_manual(values = c(
    "Downregulated" = "darkblue",
    "Not Significant" = "pink",
    "Upregulated" = "olivedrab"
  )) +

  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +

  theme(text = element_text(size = 14),
        legend.position = "right")