rm(list = ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0) {
  for (i in 1:la)
    eval(parse(text = args[[i]]))
}

require(DESeq2)
require(gplots)
require(dplyr)
require(magrittr)
require(stringr)
source(polyA_lib)

## Load Peaks file
peaks <- read.csv(peakfile, header = TRUE, sep = "\t", check.names = FALSE)
rownames(peaks) <- peaks$peak

## Get condition
samples_id <- colnames(peaks)[10:(ncol(peaks) - 1)]

design <- read.csv(input_list, sep = "\t", check.names = FALSE, header = FALSE, colClasses = c(rep("character",3)))

if (ncol(design == 4)) {
  condition <- design[match(samples_id, design[,1]), 3]
  group <- design[match(samples_id, design[,1]), 4]
}else{
  print("the input_list_compare file is false")
}

## If we want to add replicates in the statistical model
#if (ncol(design == 5)) {
#  condition <- design[match(samples_id, design[,1]), 3]
#  replicate <- design[match(samples_id, design[,1]), 4]
#  group <- design[match(samples_id, design[,1]), 5]
#}else{
#  print("the input_list_compare file is false")
#}

cond <- c(unique(condition[group==1]), unique(condition[group==0]))
condition <- as.factor(condition)
print(samples_id)
print(group)
print(condition)
print(cond)

## Filering
### coverage sum by condition > MIN_COUNT_PER_COND
peaks <- peaks[which(rowSums(peaks[,samples_id[group == 0]]) >= as.numeric(min_count_cond) | rowSums(peaks[,samples_id[group == 1]]) >= as.numeric(min_count_cond)),]
print(dim(peaks))

## keep gene with at least 2 peaks
peaks <- peaks[peaks$gene %in% peaks$gene[duplicated(peaks$gene)],]
peaks$gene <- factor(peaks$gene)
print(dim(peaks))

## Estimate size factor
ddsFirst <- DESeqDataSetFromMatrix(countData = as.matrix(peaks[,samples_id]), colData = data.frame(condition), design = ~ condition)
ddsNorm <- estimateSizeFactors(ddsFirst)
print(ddsNorm)

## DESEQ
## peakIntron vs LastExon (LE sum)
Data4DESEQ <- group_by(peaks, gene) %>% do(res=PeakIntron_LastPeakSum(.)) %>% extract2(2) %>% bind_rows
rownames(Data4DESEQ) <- paste0(Data4DESEQ$chr,":",Data4DESEQ$start,"-",Data4DESEQ$end,":",Data4DESEQ$peak)

LE <- c(rep("NLE", length(condition)), rep("LE", length(condition)))
#repli <- rep(replicate, 2)		## If we want to add replicates in the statistical model
#colData <- data.frame(t(rbind("condition" = as.vector(condition), LE, "replicate" = as.vector(repli))))		## If we want to add replicates in the statistical model
colData <- data.frame(t(rbind("condition" = as.vector(condition), LE)))
myData <- as.matrix(Data4DESEQ[,8:(ncol(Data4DESEQ))])
rownames(myData) <- rownames(Data4DESEQ)


## Statistical model
dds <- DESeqDataSetFromMatrix(countData = myData, colData = colData, design = ~ LE + condition + LE:condition)
#dds <- DESeqDataSetFromMatrix(countData = myData, colData = colData, design = ~ LE + condition + LE:condition + replicate)	## Statistical model with replicate 
sizeFactors(dds) <- rep(sizeFactors(ddsNorm), 2)
print(sizeFactors(dds))
mydds <- DESeq(dds)
res <- results(mydds,contrast = c(0,0,0,1))
#res <- results(mydds,contrast = c(0,0,0,0,0,1))		## Statistical model with replicate 


## Histogramm of pvalue & padj
pdf(file = paste(dirname(peakfile), "/hist_pval_padj.pdf", sep=""))
hist(res$pvalue, breaks=100, main="Histogramm of pvalue", xlab="pvalue")
hist(res$padj, breaks=100, main="Histogramm of padj", xlab="padj")
dev.off()

## Hierarchical clustering
cdslog <- log(counts(dds)+1)
dist.cor <- 1-cor(cdslog, use="pairwise.complete.obs", method="spearman")
hist.cor <- hclust(as.dist(dist.cor), method="ward.D")
pdf(file = paste(dirname(peakfile), "/Ascending_hierarchical_classification.pdf", sep=""))
plot(hist.cor, main=paste("samples classification by \n", nrow(cdslog), "peaks", sep=""), sub="Distance= 1-correlation")
dev.off()

NLE_cond1 <- which(colData[,"condition"] == cond[1] & colData[,"LE"] == "NLE")
NLE_cond2 <- which(colData[,"condition"] == cond[2] & colData[,"LE"] == "NLE")
LE_cond1 <- which(colData[,"condition"] == cond[1] & colData[,"LE"] == "LE")
LE_cond2 <- which(colData[,"condition"] == cond[2] & colData[,"LE"] == "LE")
logFC_intro <- log2(rowMeans(counts(dds,normalized = TRUE)[,NLE_cond1])/rowMeans(counts(dds,normalized = TRUE)[,NLE_cond2]))
logFC_LE <- log2(rowMeans(counts(dds,normalized = TRUE)[,LE_cond1])/rowMeans(counts(dds,normalized = TRUE)[,LE_cond2]))

resO <- cbind(as.matrix(res),counts(dds,normalized = TRUE), "logFC_intro" = logFC_intro, "logFC_LE" = logFC_LE)
colnames(resO)[6 + c(1:ncol(myData))] <- colnames(myData)

## Create a complete final file
resOTot <- cbind(resO,"logFC_intro-logFC_LE" = (resO[,"logFC_intro"] - resO[,"logFC_LE"]))
splitNames <- str_split_fixed(rownames(resOTot),"_",7)
resOFinal <- cbind(resOTot,"gene" = splitNames[,5])
outfile <- sub(".bed", "_peakIntron_SumLE.res", peakfile)
write.table(resOFinal, file = outfile, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

## Up file
resOUp <- resOFinal[which(as.numeric(resOFinal[,"padj"]) < 0.05 & as.numeric(resOFinal[,"logFC_intro-logFC_LE"]) > 0),]
outUpfile <- sub(".bed", "_peakIntron_SumLE_UpPeaks.res", peakfile)
write.table(resOUp, file = outUpfile, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)


## Make plots
pdf(file = paste(dirname(peakfile), "/diffan_heatmap.pdf", sep=""))

par(font.lab = 2, cex = .8)
plot(x = logFC_LE, y = logFC_intro, pch = 20, main = paste("Intronic peak vs Last Exon", paste(cond[1], cond[2], sep = "-"), "(FDR 5%)"), xlab = "log2 FC Last Exon", ylab = "log2 FC Intronic peak", frame = FALSE)
points(x = logFC_LE[which(res$padj < 0.05 & logFC_intro - logFC_LE > 0)], y = logFC_intro[which(res$padj < 0.05 & logFC_intro - logFC_LE > 0)], pch = 20, col = "red")
points(x = logFC_LE[which(res$padj < 0.05 & logFC_intro - logFC_LE < 0)],y = logFC_intro[which(res$padj < 0.05 & logFC_intro - logFC_LE < 0)], pch = 20, col = "blue")
up <- length(logFC_LE[which(res$padj < 0.05 & logFC_intro - logFC_LE > 0)])
down <- length(logFC_LE[which(res$padj < 0.05 & logFC_intro - logFC_LE < 0)])
noDiff <- length(logFC_LE[which(res$padj > 0.05)])
legend("topleft", legend = up, text.col = "red", bty = "n", cex = 1.5)
legend("bottomright", legend = down, text.col = "blue", bty = "n", cex = 1.5)
legend("topright", legend = noDiff, text.col = "black", bty = "n", cex = 1.5)

plot(x = logFC_LE, y = logFC_intro, pch = 20, main = paste("Intronic peak vs Last Exon", paste(cond[1], cond[2], sep = "-"), "(FDR 10%)"), xlab = "log2 FC Last Exon", ylab = "log2 FC Intronic peak", frame = FALSE)
points(x = logFC_LE[which(res$padj < 0.1 & logFC_intro - logFC_LE > 0)], y = logFC_intro[which(res$padj < 0.1 & logFC_intro - logFC_LE > 0)], pch = 20, col = "red")
points(x = logFC_LE[which(res$padj < 0.1 & logFC_intro - logFC_LE < 0)], y = logFC_intro[which(res$padj < 0.1 & logFC_intro - logFC_LE < 0)], pch = 20, col = "blue")
up <- length(logFC_LE[which(res$padj < 0.1 & logFC_intro - logFC_LE > 0)])
down <- length(logFC_LE[which(res$padj < 0.1 & logFC_intro - logFC_LE < 0)])
noDiff <- length(logFC_LE[which(res$padj > 0.1)])
legend("topleft", legend = up, text.col = "red", bty = "n", cex = 1.5)
legend("bottomright", legend = down, text.col = "blue", bty = "n", cex = 1.5)
legend("topright", legend = noDiff, text.col = "black", bty = "n", cex = 1.5)

plot(x = logFC_LE, y = logFC_intro, pch = 20, main = paste("Intronic peak vs Last Exon", paste(cond[1], cond[2], sep = "-"), "(pvalue 5%)"), xlab = "log2 FC Last Exon", ylab = "log2 FC Intronic peak", frame = FALSE)
points(x = logFC_LE[which(res$pvalue < 0.05 & logFC_intro - logFC_LE > 0)], y = logFC_intro[which(res$pvalue < 0.05 & logFC_intro - logFC_LE > 0)], pch = 20, col = "red")
points(x = logFC_LE[which(res$pvalue < 0.05 & logFC_intro - logFC_LE < 0)], y = logFC_intro[which(res$pvalue < 0.05 & logFC_intro - logFC_LE < 0)], pch = 20, col = "blue")
up <- length(logFC_LE[which(res$pvalue < 0.05 & logFC_intro - logFC_LE > 0)])
down <-  length(logFC_LE[which(res$pvalue < 0.05 & logFC_intro - logFC_LE < 0)])
noDiff <- length(logFC_LE[which(res$pvalue > 0.05)])
legend("topleft", legend = up, text.col = "red", bty = "n", cex = 1.5)
legend("bottomright", legend = down, text.col = "blue", bty = "n", cex = 1.5)
legend("topright", legend = noDiff, text.col = "black", bty = "n", cex = 1.5)

dev.off()

