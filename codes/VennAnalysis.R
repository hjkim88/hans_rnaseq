### This is a simple Venn analysis
###
### Two comparisons:
### A. Adult High vs Low Ca2+
### B. Fetal High vs Low Ca2+ 
###
### We would like to  know which genes are in:
### 1. A (FDR < 0.05) & B (P > 0.1)
### 2. A (P > 0.1) & B (FDR < 0.05)
### 3. A (FDR < 0.05) & B (FDR < 0.05)

### load library
if(!require(VennDiagram, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("VennDiagram")
  require(VennDiagram, quietly = TRUE)
}
if(!require(gridExtra, quietly = TRUE)) {
  install.packages("gridExtra")
  library(gridExtra, quietly = TRUE)
}

### the path of the DE results
path_A <- "C:/Research/CUMC/Hans_RNASeq/results/differential_expression/ABM_High_Ca2_vs_ABM_Low_Ca2/DE_result_ABM_High_Ca2_vs_ABM_Low_Ca2.txt"
path_B <- "C:/Research/CUMC/Hans_RNASeq/results/differential_expression/FL_High_Ca2_vs_FL_Low_Ca2/DE_result_FL_High_Ca2_vs_FL_Low_Ca2.txt"

### the result path
result_path <- "C:/Research/CUMC/Hans_RNASeq/results/venn/"

### load the data
A <- read.table(path_A, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
B <- read.table(path_B, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

### 1. A (FDR < 0.05) & B (P > 0.1)
genes <- intersect(A$Gene_Symbol[which(A$adj.P.Val < 0.05)], B$Gene_Symbol[which(B$P.Value > 0.1)])
write.table(genes, file = paste0(result_path, "AD_FDR_0.05_FT_P_0.1.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
v1 <- venn.diagram(list(A$Gene_Symbol[which(A$adj.P.Val < 0.05)], B$Gene_Symbol[which(B$P.Value > 0.1)]),
                   category.names = c("Adult\nHigh vs Low Ca2+\nFDR < 0.05",
                                      "Fetal\nHigh vs Low Ca2+\nP > 0.1"),
                   cat.cex = 1.5, cat.pos = 1, cat.dist = 0.035, cex = 1.5,
                   filename = NULL)
### save the diagram as png
png(paste0(result_path, "AD_FDR_0.05_FT_P_0.1.png"),
    width = 2000, height = 1500, res = 150)
grid.arrange(gTree(children=v1),
             top=paste0("Adult High vs Low FDR < 0.05 & Fetal High vs Low P > 0.1"),
             bottom="")
dev.off()

### 2. A (P > 0.1) & B (FDR < 0.05)
genes <- intersect(A$Gene_Symbol[which(A$P.Value > 0.1)], B$Gene_Symbol[which(B$adj.P.Val < 0.05)])
write.table(genes, file = paste0(result_path, "AD_P_0.1_FT_FDR_0.05.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
v2 <- venn.diagram(list(A$Gene_Symbol[which(A$P.Value > 0.1)], B$Gene_Symbol[which(B$adj.P.Val < 0.05)]),
                   category.names = c("Adult\nHigh vs Low Ca2+\nP > 0.1",
                                      "Fetal\nHigh vs Low Ca2+\nFDR < 0.05"),
                   cat.cex = 1.5, cat.pos = 1, cat.dist = 0.035, cex = 1.5,
                   filename = NULL)
### save the diagram as png
png(paste0(result_path, "AD_P_0.1_FT_FDR_0.05.png"),
    width = 2000, height = 1500, res = 150)
grid.arrange(gTree(children=v2),
             top=paste0("Adult High vs Low P > 0.1 & Fetal High vs Low FDR < 0.05"),
             bottom="")
dev.off()

### 3. A (FDR < 0.05) & B (FDR < 0.05)
genes <- intersect(A$Gene_Symbol[which(A$adj.P.Val < 0.05)], B$Gene_Symbol[which(B$adj.P.Val < 0.05)])
write.table(genes, file = paste0(result_path, "AD_FDR_0.05_FT_FDR_0.05.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
v3 <- venn.diagram(list(A$Gene_Symbol[which(A$adj.P.Val < 0.05)], B$Gene_Symbol[which(B$adj.P.Val < 0.05)]),
                   category.names = c("Adult\nHigh vs Low Ca2+\nFDR < 0.05",
                                      "Fetal\nHigh vs Low Ca2+\nFDR < 0.05"),
                   cat.cex = 1.5, cat.pos = 1, cat.dist = 0.035, cex = 1.5,
                   filename = NULL)
### save the diagram as png
png(paste0(result_path, "AD_FDR_0.05_FT_FDR_0.05.png"),
    width = 2000, height = 1500, res = 150)
grid.arrange(gTree(children=v3),
             top=paste0("Adult High vs Low FDR < 0.05 & Fetal High vs Low FDR < 0.05"),
             bottom="")
dev.off()
