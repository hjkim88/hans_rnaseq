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
if(!require(org.Mm.eg.db, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Mm.eg.db")
  require(org.Mm.eg.db, quietly = TRUE)
}
if(!require(gridExtra, quietly = TRUE)) {
  install.packages("gridExtra")
  library(gridExtra, quietly = TRUE)
}
### load library
if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}

### the path of the DE results
path_A <- "C:/Research/CUMC/Hans_RNASeq/results/differential_expression/ABM_High_Ca2_vs_ABM_Low_Ca2/DE_result_ABM_High_Ca2_vs_ABM_Low_Ca2.txt"
path_B <- "C:/Research/CUMC/Hans_RNASeq/results/differential_expression/FL_High_Ca2_vs_FL_Low_Ca2/DE_result_FL_High_Ca2_vs_FL_Low_Ca2.txt"
path_C <- "C:/Research/CUMC/Hans_RNASeq/results/differential_expression/(ABM_High_Ca2 - ABM_Low_Ca2)_vs_(FL_High_Ca2 - FL_Low_Ca2)/DE_result_(ABM_High_Ca2 - ABM_Low_Ca2)_vs_(FL_High_Ca2 - FL_Low_Ca2).txt"

### the result path
result_path <- "C:/Research/CUMC/Hans_RNASeq/results/venn/"

### load the data
A <- read.table(path_A, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
B <- read.table(path_B, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
C <- read.table(path_C, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

### 1. A (FDR < 0.05) & B (P > 0.1)
a <- A$Gene_Symbol[which(A$adj.P.Val < 0.05)]
b <- B$Gene_Symbol[which(B$P.Value > 0.1)]
remove_idx <- unique(c(which(is.na(a)), which(is.null(a)), which(a == "NULL"), which(a == "NA")))
if(length(remove_idx) > 0) {
  a <- a[-remove_idx]
}
remove_idx <- unique(c(which(is.na(b)), which(is.null(b)), which(b == "NULL"), which(b == "NA")))
if(length(remove_idx) > 0) {
  b <- b[-remove_idx]
}
genes <- intersect(a, b)
result_table <- data.frame(Gene_Symbol=genes,
                           Gene_Name=mapIds(org.Mm.eg.db, genes, c("GENENAME"), "ALIAS"),
                           logFC_ABM_High_Ca2_vs_Low_Ca2BM_High_Ca2_vs_Low_Ca2=A$logFC[which(A$Gene_Symbol %in% genes)],
                           logFC_FL_High_Ca2_vs_Low_Ca2=B$logFC[which(B$Gene_Symbol %in% genes)],
                           PVal_ABM_High_Ca2_vs_Low_Ca2=A$P.Value[which(A$Gene_Symbol %in% genes)],
                           PVal_FL_High_Ca2_vs_Low_Ca2=B$P.Value[which(B$Gene_Symbol %in% genes)],
                           FDR_ABM_High_Ca2_vs_Low_Ca2=A$adj.P.Val[which(A$Gene_Symbol %in% genes)],
                           FDR_FL_High_Ca2_vs_Low_Ca2=B$adj.P.Val[which(B$Gene_Symbol %in% genes)],
                           PVal_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=C$P.Value[which(C$Gene_Symbol %in% genes)],
                           FDR_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=C$adj.P.Val[which(C$Gene_Symbol %in% genes)],
                           FDR2_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=p.adjust(C$P.Value[which(C$Gene_Symbol %in% genes)], method = "BH"))
write.xlsx2(result_table, file = paste0(result_path, "AD_FDR_0.05_FT_P_0.1.xlsx"),
            sheetName = "Venn_Result", row.names = FALSE)
v1 <- venn.diagram(list(a, b),
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
a <- A$Gene_Symbol[which(A$P.Value > 0.1)]
b <- B$Gene_Symbol[which(B$adj.P.Val < 0.05)]
remove_idx <- unique(c(which(is.na(a)), which(is.null(a)), which(a == "NULL"), which(a == "NA")))
if(length(remove_idx) > 0) {
  a <- a[-remove_idx]
}
remove_idx <- unique(c(which(is.na(b)), which(is.null(b)), which(b == "NULL"), which(b == "NA")))
if(length(remove_idx) > 0) {
  b <- b[-remove_idx]
}
genes <- intersect(a, b)
result_table <- data.frame(Gene_Symbol=genes,
                           Gene_Name=mapIds(org.Mm.eg.db, genes, c("GENENAME"), "ALIAS"),
                           logFC_ABM_High_Ca2_vs_Low_Ca2=A$logFC[which(A$Gene_Symbol %in% genes)],
                           logFC_FL_High_Ca2_vs_Low_Ca2=B$logFC[which(B$Gene_Symbol %in% genes)],
                           PVal_ABM_High_Ca2_vs_Low_Ca2=A$P.Value[which(A$Gene_Symbol %in% genes)],
                           PVal_FL_High_Ca2_vs_Low_Ca2=B$P.Value[which(B$Gene_Symbol %in% genes)],
                           FDR_ABM_High_Ca2_vs_Low_Ca2=A$adj.P.Val[which(A$Gene_Symbol %in% genes)],
                           FDR_FL_High_Ca2_vs_Low_Ca2=B$adj.P.Val[which(B$Gene_Symbol %in% genes)],
                           PVal_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=C$P.Value[which(C$Gene_Symbol %in% genes)],
                           FDR_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=C$adj.P.Val[which(C$Gene_Symbol %in% genes)],
                           FDR2_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=p.adjust(C$P.Value[which(C$Gene_Symbol %in% genes)], method = "BH"))
write.xlsx2(result_table, file = paste0(result_path, "AD_P_0.1_FT_FDR_0.05.xlsx"),
            sheetName = "Venn_Result", row.names = FALSE)
v2 <- venn.diagram(list(a, b),
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
a <- A$Gene_Symbol[which(A$adj.P.Val < 0.05)]
b <- B$Gene_Symbol[which(B$adj.P.Val < 0.05)]
remove_idx <- unique(c(which(is.na(a)), which(is.null(a)), which(a == "NULL"), which(a == "NA")))
if(length(remove_idx) > 0) {
  a <- a[-remove_idx]
}
remove_idx <- unique(c(which(is.na(b)), which(is.null(b)), which(b == "NULL"), which(b == "NA")))
if(length(remove_idx) > 0) {
  b <- b[-remove_idx]
}
genes <- intersect(a, b)
result_table <- data.frame(Gene_Symbol=genes,
                           Gene_Name=mapIds(org.Mm.eg.db, genes, c("GENENAME"), "ALIAS"),
                           logFC_ABM_High_Ca2_vs_Low_Ca2=A$logFC[which(A$Gene_Symbol %in% genes)],
                           logFC_FL_High_Ca2_vs_Low_Ca2=B$logFC[which(B$Gene_Symbol %in% genes)],
                           PVal_ABM_High_Ca2_vs_Low_Ca2=A$P.Value[which(A$Gene_Symbol %in% genes)],
                           PVal_FL_High_Ca2_vs_Low_Ca2=B$P.Value[which(B$Gene_Symbol %in% genes)],
                           FDR_ABM_High_Ca2_vs_Low_Ca2=A$adj.P.Val[which(A$Gene_Symbol %in% genes)],
                           FDR_FL_High_Ca2_vs_Low_Ca2=B$adj.P.Val[which(B$Gene_Symbol %in% genes)],
                           PVal_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=C$P.Value[which(C$Gene_Symbol %in% genes)],
                           FDR_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=C$adj.P.Val[which(C$Gene_Symbol %in% genes)],
                           FDR2_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=p.adjust(C$P.Value[which(C$Gene_Symbol %in% genes)], method = "BH"))
write.xlsx2(result_table, file = paste0(result_path, "AD_FDR_0.05_FT_FDR_0.05.xlsx"),
            sheetName = "Venn_Result", row.names = FALSE)
### ++
write.xlsx2(result_table[intersect(which(result_table$logFC_ABM_High_Ca2_vs_Low_Ca2 > 0), which(result_table$logFC_FL_High_Ca2_vs_Low_Ca2 > 0)),],
            file = paste0(result_path, "AD_FDR_0.05_FT_FDR_0.05_logFC++.xlsx"),
            sheetName = "Venn_Result", row.names = FALSE)
### --
write.xlsx2(result_table[intersect(which(result_table$logFC_ABM_High_Ca2_vs_Low_Ca2 < 0), which(result_table$logFC_FL_High_Ca2_vs_Low_Ca2 < 0)),],
            file = paste0(result_path, "AD_FDR_0.05_FT_FDR_0.05_logFC--.xlsx"),
            sheetName = "Venn_Result", row.names = FALSE)
### +-
write.xlsx2(result_table[intersect(which(result_table$logFC_ABM_High_Ca2_vs_Low_Ca2 > 0), which(result_table$logFC_FL_High_Ca2_vs_Low_Ca2 < 0)),],
            file = paste0(result_path, "AD_FDR_0.05_FT_FDR_0.05_logFC+-.xlsx"),
            sheetName = "Venn_Result", row.names = FALSE)
### -+
write.xlsx2(result_table[intersect(which(result_table$logFC_ABM_High_Ca2_vs_Low_Ca2 < 0), which(result_table$logFC_FL_High_Ca2_vs_Low_Ca2 > 0)),],
            file = paste0(result_path, "AD_FDR_0.05_FT_FDR_0.05_logFC-+.xlsx"),
            sheetName = "Venn_Result", row.names = FALSE)
v3 <- venn.diagram(list(a, b),
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


### second Venn analysis
### Venn analysis with the signature genes from Wilson et al.

### load the signature genes
sig_genes1 <- read.table(file = "C:/Research/CUMC/Hans_RNASeq/data/genes_of_interest/MolO_vs_NoMO.txt",
                         sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)[,1]
sig_genes2 <- read.table(file = "C:/Research/CUMC/Hans_RNASeq/data/genes_of_interest/DE_signature_genes.txt",
                         sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)[,1]
sig_genes3 <- read.table(file = "C:/Research/CUMC/Hans_RNASeq/data/genes_of_interest/DE_signature_genes2.txt",
                         sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)[,1]
sig_genes <- unique(union(union(sig_genes1, sig_genes2), sig_genes3))

### annotate our exp's DE statistics to the genes
sig_result <- data.frame(Gene_Symbol=sig_genes,
                         Gene_Name=mapIds(org.Mm.eg.db, sig_genes, c("GENENAME"), "ALIAS"),
                         Appeared=sapply(sig_genes, function(x) {
                           y <- ""
                           if(length(which(sig_genes1 == x)) > 0) {
                             y <- paste(y, "Fig2", sep = "/")
                           }
                           if(length(which(sig_genes2 == x)) > 0) {
                             y <- paste(y, "Fig3", sep = "/")
                           }
                           if(length(which(sig_genes3 == x)) > 0) {
                             y <- paste(y, "Fig5", sep = "/")
                           }
                           return(substring(y, 2))
                         }),
                         stringsAsFactors = FALSE, check.names = FALSE)
sig_result <- merge(sig_result, A[,c("Gene_Symbol", "logFC", "P.Value", "adj.P.Val")],
                    by.x = "Gene_Symbol", by.y = "Gene_Symbol",
                    all.x = TRUE)
sig_result <- merge(sig_result, B[,c("Gene_Symbol", "logFC", "P.Value", "adj.P.Val")],
                    by.x = "Gene_Symbol", by.y = "Gene_Symbol",
                    all.x = TRUE)
sig_result <- merge(sig_result, C[,c("Gene_Symbol", "P.Value", "adj.P.Val")],
                    by.x = "Gene_Symbol", by.y = "Gene_Symbol",
                    all.x = TRUE)

colnames(sig_result) <- c("Gene_Symbol", "Gene_Name", "Source in Wilson paper",
                          "logFC_ABM_High_Ca2_vs_Low_Ca2", "PVal_ABM_High_Ca2_vs_Low_Ca2",
                          "FDR_ABM_High_Ca2_vs_Low_Ca2", "logFC_FL_High_Ca2_vs_Low_Ca2",
                          "PVal_FL_High_Ca2_vs_Low_Ca2", "FDR_FL_High_Ca2_vs_Low_Ca2",
                          "PVal_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2",
                          "FDR_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2")

sig_result <- data.frame(sig_result,
                         FDR2_ABM_High_Ca2_vs_Low_Ca2_VS_FL_High_Ca2_vs_Low_Ca2=p.adjust(sig_result[,10], method = "BH"),
                         stringsAsFactors = FALSE, check.names = FALSE)

### print out the results
write.xlsx2(sig_result, file = paste0(result_path, "Wilson_Signature_Genes.xlsx"),
            sheetName = "Venn_Result", row.names = FALSE)
