### This is an easy version of the "DEAnalysis.R" of Corey's (Snoeck) Project

### load raw count data
load("//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Hans/data/raw_counts.rda")

### make sample info
sampleType <- as.factor(rep(c("ABM_High_Ca2", "ABM_Low_Ca2", "FL_High_Ca2", "FL_Low_Ca2"), 3))
sampleType1 <- as.factor(rep(c(rep("ABM", 2), rep("FL", 2)), 3))
sampleType2 <- as.factor(rep(c("High_Ca2", "Low_Ca2"), 6))

### load library
if(!require(edgeR)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  library(edgeR)
}
if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}

### make a data frame for design matrix
Coldata <- data.frame(sampleType)
rownames(Coldata) <- colnames(rawCnt[,-c(1,2)])

### prepare design matrix
design <- model.matrix(~0+sampleType, data = Coldata)
colnames(design) <- levels(sampleType)

### raw counts to DGEList object
d <- DGEList(rawCnt[,-c(1,2)])

### remove 0 or low counts
keep <- filterByExpr(d, design)
d <- d[keep,,keep.lib.sizes=FALSE]

### calculate normalization factors
d <- calcNormFactors(d)

### voom
d <- voomWithQualityWeights(d, design)

### fit the linear model
fit <- lmFit(d, design)

#   1. ABM_High_Ca2+ vs ABM_Low_Ca2+
### extract specific comparisons of interest
contrastMat <- makeContrasts(contrasts = "ABM_High_Ca2-ABM_Low_Ca2", levels = design)
### fit the contrasts
fit2 <- contrasts.fit(fit, contrastMat)
fit2 <- eBayes(fit2)
### get the differentially expressed genes
result <- topTable(fit2, adjust.method="BH", number=Inf)
### order based on adj.p.val
result <- result[order(result$adj.P.Val),]
### add baseMean for each group
exp_rowMeans <- apply(d$E[,which(Coldata$sampleType == "ABM_High_Ca2"), drop=FALSE], 1, mean)
ctrl_rowMeans <- apply(d$E[,which(Coldata$sampleType == "ABM_Low_Ca2"), drop=FALSE], 1, mean)
result <- data.frame(V1=exp_rowMeans[rownames(result)],
                     V2=ctrl_rowMeans[rownames(result)],
                     result[,1:6],
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(result)[1:2] <- c(paste0("normMean_", "ABM_High_Ca2"), paste0("normMean_", "ABM_Low_Ca2")) 
### write out the DE result tables
write.xlsx2(data.frame(Entrez_ID=rownames(result), Gene_Symbol=rawCnt[rownames(result),"Gene_Symbol"], result,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = "./DE_result_ABM_High_Ca2_vs_ABM_Low_Ca2.xlsx",
            sheetName = "DE_Result", row.names = FALSE)

#   2. FL_High_Ca2+ vs FL_Low_Ca2+
### extract specific comparisons of interest
contrastMat <- makeContrasts(contrasts = "FL_High_Ca2-FL_Low_Ca2", levels = design)
### fit the contrasts
fit2 <- contrasts.fit(fit, contrastMat)
fit2 <- eBayes(fit2)
### get the differentially expressed genes
result <- topTable(fit2, adjust.method="BH", number=Inf)
### order based on adj.p.val
result <- result[order(result$adj.P.Val),]
### add baseMean for each group
exp_rowMeans <- apply(d$E[,which(Coldata$sampleType == "FL_High_Ca2"), drop=FALSE], 1, mean)
ctrl_rowMeans <- apply(d$E[,which(Coldata$sampleType == "FL_Low_Ca2"), drop=FALSE], 1, mean)
result <- data.frame(V1=exp_rowMeans[rownames(result)],
                     V2=ctrl_rowMeans[rownames(result)],
                     result[,1:6],
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(result)[1:2] <- c(paste0("normMean_", "FL_High_Ca2"), paste0("normMean_", "FL_Low_Ca2")) 
### write out the DE result tables
write.xlsx2(data.frame(Entrez_ID=rownames(result), Gene_Symbol=rawCnt[rownames(result),"Gene_Symbol"], result,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = "./DE_result_FL_High_Ca2_vs_FL_Low_Ca2.xlsx",
            sheetName = "DE_Result", row.names = FALSE)

#   3. ABM_High_Ca2+ vs  FL_High_Ca2+
### extract specific comparisons of interest
contrastMat <- makeContrasts(contrasts = "ABM_High_Ca2-FL_High_Ca2", levels = design)
### fit the contrasts
fit2 <- contrasts.fit(fit, contrastMat)
fit2 <- eBayes(fit2)
### get the differentially expressed genes
result <- topTable(fit2, adjust.method="BH", number=Inf)
### order based on adj.p.val
result <- result[order(result$adj.P.Val),]
### add baseMean for each group
exp_rowMeans <- apply(d$E[,which(Coldata$sampleType == "ABM_High_Ca2"), drop=FALSE], 1, mean)
ctrl_rowMeans <- apply(d$E[,which(Coldata$sampleType == "FL_High_Ca2"), drop=FALSE], 1, mean)
result <- data.frame(V1=exp_rowMeans[rownames(result)],
                     V2=ctrl_rowMeans[rownames(result)],
                     result[,1:6],
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(result)[1:2] <- c(paste0("normMean_", "ABM_High_Ca2"), paste0("normMean_", "FL_High_Ca2")) 
### write out the DE result tables
write.xlsx2(data.frame(Entrez_ID=rownames(result), Gene_Symbol=rawCnt[rownames(result),"Gene_Symbol"], result,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = "./DE_result_ABM_High_Ca2_vs_FL_High_Ca2.xlsx",
            sheetName = "DE_Result", row.names = FALSE)

#   4. ABM_Low_Ca2+ vs FL_Low_Ca2+
### extract specific comparisons of interest
contrastMat <- makeContrasts(contrasts = "ABM_Low_Ca2-FL_Low_Ca2", levels = design)
### fit the contrasts
fit2 <- contrasts.fit(fit, contrastMat)
fit2 <- eBayes(fit2)
### get the differentially expressed genes
result <- topTable(fit2, adjust.method="BH", number=Inf)
### order based on adj.p.val
result <- result[order(result$adj.P.Val),]
### add baseMean for each group
exp_rowMeans <- apply(d$E[,which(Coldata$sampleType == "ABM_Low_Ca2"), drop=FALSE], 1, mean)
ctrl_rowMeans <- apply(d$E[,which(Coldata$sampleType == "FL_Low_Ca2"), drop=FALSE], 1, mean)
result <- data.frame(V1=exp_rowMeans[rownames(result)],
                     V2=ctrl_rowMeans[rownames(result)],
                     result[,1:6],
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(result)[1:2] <- c(paste0("normMean_", "ABM_Low_Ca2"), paste0("normMean_", "FL_Low_Ca2"))
### write out the DE result tables
write.xlsx2(data.frame(Entrez_ID=rownames(result), Gene_Symbol=rawCnt[rownames(result),"Gene_Symbol"], result,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = "./DE_result_ABM_Low_Ca2_vs_FL_Low_Ca2.xlsx",
            sheetName = "DE_Result", row.names = FALSE)

#   5. ABM_High_Ca2+ vs FL_Low_Ca2+
### extract specific comparisons of interest
contrastMat <- makeContrasts(contrasts = "ABM_High_Ca2-FL_Low_Ca2", levels = design)
### fit the contrasts
fit2 <- contrasts.fit(fit, contrastMat)
fit2 <- eBayes(fit2)
### get the differentially expressed genes
result <- topTable(fit2, adjust.method="BH", number=Inf)
### order based on adj.p.val
result <- result[order(result$adj.P.Val),]
### add baseMean for each group
exp_rowMeans <- apply(d$E[,which(Coldata$sampleType == "ABM_High_Ca2"), drop=FALSE], 1, mean)
ctrl_rowMeans <- apply(d$E[,which(Coldata$sampleType == "FL_Low_Ca2"), drop=FALSE], 1, mean)
result <- data.frame(V1=exp_rowMeans[rownames(result)],
                     V2=ctrl_rowMeans[rownames(result)],
                     result[,1:6],
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(result)[1:2] <- c(paste0("normMean_", "ABM_High_Ca2"), paste0("normMean_", "FL_Low_Ca2"))
### write out the DE result tables
write.xlsx2(data.frame(Entrez_ID=rownames(result), Gene_Symbol=rawCnt[rownames(result),"Gene_Symbol"], result,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = "./DE_result_ABM_High_Ca2_vs_FL_Low_Ca2.xlsx",
            sheetName = "DE_Result", row.names = FALSE)

#   6. ABM_Low_Ca2+ vs FL_High_Ca2+
### extract specific comparisons of interest
contrastMat <- makeContrasts(contrasts = "ABM_Low_Ca2-FL_High_Ca2", levels = design)
### fit the contrasts
fit2 <- contrasts.fit(fit, contrastMat)
fit2 <- eBayes(fit2)
### get the differentially expressed genes
result <- topTable(fit2, adjust.method="BH", number=Inf)
### order based on adj.p.val
result <- result[order(result$adj.P.Val),]
### add baseMean for each group
exp_rowMeans <- apply(d$E[,which(Coldata$sampleType == "ABM_Low_Ca2"), drop=FALSE], 1, mean)
ctrl_rowMeans <- apply(d$E[,which(Coldata$sampleType == "FL_High_Ca2"), drop=FALSE], 1, mean)
result <- data.frame(V1=exp_rowMeans[rownames(result)],
                     V2=ctrl_rowMeans[rownames(result)],
                     result[,1:6],
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(result)[1:2] <- c(paste0("normMean_", "ABM_Low_Ca2"), paste0("normMean_", "FL_High_Ca2"))
### write out the DE result tables
write.xlsx2(data.frame(Entrez_ID=rownames(result), Gene_Symbol=rawCnt[rownames(result),"Gene_Symbol"], result,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = "./DE_result_ABM_Low_Ca2_vs_FL_High_Ca2.xlsx",
            sheetName = "DE_Result", row.names = FALSE)


#   7. (ABM_High_Ca2+ vs ABM_Low_Ca2+) vs (FL_High_Ca2+ vs FL_Low_Ca2+)
### make a data frame for design matrix
Coldata <- data.frame(sampleType1, sampleType2)
rownames(Coldata) <- colnames(rawCnt[,-c(1,2)])
Coldata$sampleType1 <- relevel(Coldata$sampleType1, ref = "FL")
Coldata$sampleType2 <- relevel(Coldata$sampleType2, ref = "Low_Ca2")

### prepare design matrix
design <- model.matrix(~0+sampleType1*sampleType2, data = Coldata)

### raw counts to DGEList object
d <- DGEList(rawCnt[,-c(1,2)])

### remove 0 or low counts
keep <- filterByExpr(d, design)
d <- d[keep,,keep.lib.sizes=FALSE]

### calculate normalization factors
d <- calcNormFactors(d)

### voom
d <- voomWithQualityWeights(d, design)

### fit the linear model
fit <- lmFit(d, design)
fit2 <- eBayes(fit)

### get the differentially expressed genes
result <- topTable(fit2, coef = "sampleType1ABM:sampleType2High_Ca2",
                   adjust.method="BH", number=Inf)

### order based on adj.p.val
result <- result[order(result$adj.P.Val),]

### add baseMean for each group
exp1_ctrl1_rowMeans <- apply(d$E[,intersect(which(Coldata$sampleType1 == "ABM"),
                                            which(Coldata$sampleType2 == "High_Ca2")),
                                 drop=FALSE], 1, mean)
exp1_ctrl2_rowMeans <- apply(d$E[,intersect(which(Coldata$sampleType1 == "ABM"),
                                            which(Coldata$sampleType2 == "Low_Ca2")),
                                 drop=FALSE], 1, mean)
exp2_ctrl1_rowMeans <- apply(d$E[,intersect(which(Coldata$sampleType1 == "FL"),
                                            which(Coldata$sampleType2 == "High_Ca2")),
                                 drop=FALSE], 1, mean)
exp2_ctrl2_rowMeans <- apply(d$E[,intersect(which(Coldata$sampleType1 == "FL"),
                                            which(Coldata$sampleType2 == "Low_Ca2")),
                                 drop=FALSE], 1, mean)
result <- data.frame(V1=exp1_ctrl1_rowMeans[rownames(result)],
                     V2=exp1_ctrl2_rowMeans[rownames(result)],
                     V3=exp2_ctrl1_rowMeans[rownames(result)],
                     V4=exp2_ctrl2_rowMeans[rownames(result)],
                     result[,1:6],
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(result)[1:4] <- c(paste0("normMean_", "ABM", "_", "High_Ca2"),
                           paste0("normMean_", "ABM", "_", "Low_Ca2"),
                           paste0("normMean_", "FL", "_", "High_Ca2"),
                           paste0("normMean_", "FL", "_", "Low_Ca2"))

### write out the DE result tables
write.xlsx2(data.frame(Entrez_ID=rownames(result), Gene_Symbol=rawCnt[rownames(result),"Gene_Symbol"], result,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = "./DE_result_(ABM_High_Ca2 - ABM_Low_Ca2)_vs_(FL_High_Ca2 - FL_Low_Ca2).xlsx",
            sheetName = "DE_Result", row.names = FALSE)
