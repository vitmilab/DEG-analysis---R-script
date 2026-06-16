# ============================================================
# External Transcriptomic Validation of Candidate Hub Genes
# Dataset: GSE80699
# ============================================================

library(ggplot2)
library(tidyr)

# Read GSE80699 files
control1 <- read.table(
  "GSE80699_Genotype1_control.Gene.rpkm.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

stress1 <- read.table(
  "GSE80699_Genotype1-treatment.Gene.rpkm.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

control2 <- read.table(
  "GSE80699_Genotype2-control.Gene.rpkm.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

stress2 <- read.table(
  "GSE80699_Genotype2-treatment.Gene.rpkm.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Candidate hub genes
genes <- c(
  "Sb01g016780.1",  # Sobic.001G191000
  "Sb01g042610.1",  # Sobic.001G453300
  "Sb04g027550.1",  # Sobic.004G231800
  "Sb10g027590.1"   # Sobic.010G233800
)

# Extract RPKM values
c1 <- control1[control1$GeneID %in% genes, c("GeneID", "RPKM")]
s1 <- stress1[stress1$GeneID %in% genes, c("GeneID", "RPKM")]
c2 <- control2[control2$GeneID %in% genes, c("GeneID", "RPKM")]
s2 <- stress2[stress2$GeneID %in% genes, c("GeneID", "RPKM")]

colnames(c1)[2] <- "IS20351_Control"
colnames(s1)[2] <- "IS20351_Drought"
colnames(c2)[2] <- "IS22330_Control"
colnames(s2)[2] <- "IS22330_Drought"

# Merge tables
validation_table <- merge(c1, s1, by = "GeneID")
validation_table <- merge(validation_table, c2, by = "GeneID")
validation_table <- merge(validation_table, s2, by = "GeneID")

# Calculate log2 fold changes
validation_table$log2FC_IS20351 <- log2(
  validation_table$IS20351_Drought /
    validation_table$IS20351_Control
)

validation_table$log2FC_IS22330 <- log2(
  validation_table$IS22330_Drought /
    validation_table$IS22330_Control
)

# Add candidate gene names
validation_table$Candidate_Gene <- c(
  "Sobic.001G191000",
  "Sobic.001G453300",
  "Sobic.004G231800",
  "Sobic.010G233800"
)

validation_table <- validation_table[, c(
  "Candidate_Gene",
  "GeneID",
  "IS20351_Control",
  "IS20351_Drought",
  "log2FC_IS20351",
  "IS22330_Control",
  "IS22330_Drought",
  "log2FC_IS22330"
)]

# Export Table S4
write.csv(
  validation_table,
  "Table_S4_GSE80699_Validation.csv",
  row.names = FALSE
)

print(validation_table)



plot_data <- data.frame(
  Gene = c(
    "Sobic.001G191000",
    "Sobic.001G453300",
    "Sobic.004G231800",
    "Sobic.010G233800",
    "Sobic.001G191000",
    "Sobic.001G453300",
    "Sobic.004G231800",
    "Sobic.010G233800"
  ),
  Genotype = c(
    "IS20351","IS20351","IS20351","IS20351",
    "IS22330","IS22330","IS22330","IS22330"
  ),
  Control = c(
    154.86,224.20,256.62,901.91,
    122.11,129.21,165.79,587.79
  ),
  Drought = c(
    65.71,74.05,62.95,310.65,
    82.78,97.07,128.23,425.89
  )
)

plot_data <- pivot_longer(
  plot_data,
  cols = c(Control, Drought),
  names_to = "Condition",
  values_to = "RPKM"
)

p <- ggplot(
  plot_data,
  aes(
    x = Gene,
    y = RPKM,
    fill = Condition
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~Genotype, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    x = "Candidate Hub Genes",
    y = "RPKM",
    fill = "Condition"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p)

ggsave(
  "Supplementary_Figure_GSE80699_Validation.tiff",
  plot = p,
  width = 8,
  height = 5,
  dpi = 600,
  compression = "lzw"
)