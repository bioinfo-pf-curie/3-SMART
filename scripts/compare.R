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
group <- strsplit(group, ",")[[1]]
## reorder 0/1
group <- group[c(which(group == 0), which(group == 1))]
samples_id <- samples_id[c(which(group == 0), which(group == 1))]
stopifnot(length(samples_id) == length(group))

design <- read.csv(input_list, sep = "\t", check.names = FALSE, header = FALSE)
if (ncol(design == 3)) {
  condition <- design[match(samples_id, design[,1]), 3]
  condition <- c(rep(paste0(unique(as.character(condition[which(group == 0)])), collapse = "_"), length(which(group == 0))),
                 rep(paste0(unique(as.character(condition[which(group == 1)])), collapse = "_"), length(which(group == 1))))
}else{
  condition <- group
}

condition <- as.factor(as.character(condition))
cond <- c(paste0(as.character(unique(condition[which(group == 0)])), collapse = "_"), paste0(as.character(unique(condition[which(group == 1)])), collapse = "_"))
print(samples_id)
print(group)
print(condition)
print(cond)

## Filering
### coverage sum by condition > MIN_COUNT_PER_COND
peaks <- peaks[which(rowSums(peaks[,samples_id[which(group == 0)]]) >= as.numeric(min_count_cond) | rowSums(peaks[,samples_id[which(group == 1)]]) >= as.numeric(min_count_cond)),]
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
Data4DESEQ <- group_by(peaks, gene) %>% do(res=PeakIntron_LastPeakSum(.)) %>% extract2(2) %>% rbind_all
rownames(Data4DESEQ) <- paste0(Data4DESEQ$chr,":",Data4DESEQ$start,"-",Data4DESEQ$end,":",Data4DESEQ$peak)

LE <- c(rep("NLE", length(condition)), rep("LE", length(condition)))
colData <- data.frame(t(rbind("condition" = as.vector(condition), LE)))
myData <- as.matrix(Data4DESEQ[,8:(ncol(Data4DESEQ))])
rownames(myData) <- rownames(Data4DESEQ)

## Ratio (ELE/LE)*100 
RatELE_LE <- (rowSums(myData[,(1:length(condition))])/rowSums(myData[,(length(condition) + 1:length(condition))]))*100

## Ratio (LE/ELE)*100 
RatLE_ELE <- (rowSums(myData[,(length(condition) + 1:length(condition))])/rowSums(myData[,(1:length(condition))]))*100

## Add this 2 ratio at myData & make a plot
NewMyData <- cbind(RatELE_LE, RatLE_ELE, myData)
#write.table(file="NewMyData_ELE_LE_ratio_For_Histogramm.txt", NewMyData, col.names=T, row.names=T, quote=F, sep="\t")

## Select only the gene with at least ratioValue
SelectPeak <- NewMyData[which(NewMyData[,1] >= as.numeric(ratioValue) & NewMyData[,2] >= as.numeric(ratioValue)),]
#write.table(file="Peak_ratio.res", SelectPeak, col.names=T, row.names=T, quote=F, sep="\t")

## Remove columns RatELE_LE & RatLE_ELE
myDataB <- SelectPeak[,-c(1:2)]
myData <- myDataB

#
dds <- DESeqDataSetFromMatrix(countData = myData, colData = colData, design = ~ LE + condition + LE:condition)
sizeFactors(dds) <- rep(sizeFactors(ddsNorm), 2)
mydds <- DESeq(dds)
res <- results(mydds,contrast = c(0,0,0,1))

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

