# =========================================================
# Differential Expression Analysis of GSE48205 Microarray Dataset
# Platform: GPL17335
# Organism: Sorghum bicolor
# Analysis Method: limma (Linear Models for Microarray Data)
# =========================================================
# =========================================================
# 1. Load Required Libraries
# =========================================================
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# =========================================================
# 2. Download and Load GEO Dataset
# =========================================================


gset <- getGEO("GSE48205", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17335", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


# =========================================================
# 3. Prepare Feature Annotation
# =========================================================

fvarLabels(gset) <- make.names(fvarLabels(gset))


## =========================================================
# 4. Define Experimental Groups
# =========================================================

# 1 = Stress-treated samples
# 0 = Control samples
# X = Excluded samples
gsms <- "111XXXXXX000"
sml <- strsplit(gsms, split="")[[1]]


# =========================================================
# 5. Filter Selected Samples
# =========================================================
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]


# =========================================================
# 6. Data Preprocessing and Log2 Transformation
# =========================================================
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }


# =========================================================
# 7. Design Matrix Construction
# =========================================================
gs <- factor(sml)
groups <- make.names(c("Test","Control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values


# =========================================================
# 8. Differential Expression Analysis Using limma
# =========================================================
fit <- lmFit(gset, design)  # fit linear model


# =========================================================
# 9. Contrast Matrix and Statistical Testing
# =========================================================
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)


# =========================================================
# 10. Extract Significant Differentially Expressed Genes
# =========================================================

# Thresholds:
# Adjusted p-value < 0.05
# |log2FC| > 1
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","ORF","SPOT_ID"))
write.csv(tT, "GSE48205_DEGs.csv", row.names = FALSE)


## =========================================================
# 11. Volcano Plot Visualization
# =========================================================
# Define significance thresholds
pvalue_cutoff <- 0.05
logFC_cutoff <- 1

# Extract differential expression values
logFC <- fit2$coefficients[, ct]
P.Value <- fit2$p.value[, ct]
logP <- -log10(P.Value)

# Classify genes based on significance thresholds
Status <- ifelse(P.Value < pvalue_cutoff & logFC > logFC_cutoff,
                 "Upregulated",
                 ifelse(P.Value < pvalue_cutoff & logFC < -logFC_cutoff,
                        "Downregulated",
                        "Not Significant"))

# Combine results into a dataframe
volcano_data <- data.frame(logFC, logP, Status)

# Generate volcano plot
volcano_plot <- ggplot(volcano_data,
                       aes(x = logFC,
                           y = logP,
                           color = Status)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "olivedrab",
                                "Downregulated" = "darkblue",
                                "Not Significant" = "pink")) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "right")

# Display volcano plot
volcano_plot