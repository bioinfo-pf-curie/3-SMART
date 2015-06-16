rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

require(gplots)
require(dplyr)
require(magrittr)
require(DESeq2)
source("scripts/polyA_lib.R")

## Load Peaks file
peaks <- read.csv(peakfile, header=TRUE,sep="\t", check.names=FALSE)
rownames(peaks) <- peaks$peak

## Get condition
samples_id <- colnames(peaks)[10:(ncol(peaks)-1)]
group <- strsplit(group,",")[[1]]
## reorder 0/1
group <- group[c(which(group==0), which(group==1))]
samples_id <- samples_id[c(which(group==0), which(group==1))]
stopifnot(length(samples_id)==length(group))

design <- read.csv(input_list, sep="\t", check.names=FALSE, header=FALSE)
if (ncol(design==3)){
  condition <- design[match(samples_id, design[,1]),3]
  condition <- c(rep(paste0(unique(as.character(condition[which(group==0)])), collapse="_"), length(which(group==0))),
                 rep(paste0(unique(as.character(condition[which(group==1)])), collapse="_"), length(which(group==1))))
}else{
  condition <-  group
}

condition <- as.factor(as.character(condition))
cond <- c(paste0(as.character(unique(condition[which(group==0)])), collapse="_"), paste0(as.character(unique(condition[which(group==1)])), collapse="_"))
print (samples_id)
print (group)
print(condition)
print(cond)

## Filering
### coverage sum by condition > MIN_COUNT_PER_COND
peaks <- peaks[which(rowSums(peaks[,samples_id[which(group==0)]]) >= as.numeric(min_count_cond) | rowSums(peaks[,samples_id[which(group==1)]]) >= as.numeric(min_count_cond)),]
print(dim(peaks))

## keep gene with at least 2 peaks
peaks <- peaks[peaks$gene %in% peaks$gene[duplicated(peaks$gene)],]
peaks$gene <- factor(peaks$gene)
print(dim(peaks))

## Estimate size factor
ddsFirst <- DESeqDataSetFromMatrix(countData=as.matrix(peaks[,samples_id]), colData=data.frame(condition),design = ~ condition)
ddsNorm <- estimateSizeFactors(ddsFirst)
print(ddsNorm)

## DESEQ
## peakIntron vs LastExon (LE sum)
Data4DESEQ <- group_by(peaks, gene) %>% do(res=PeakIntron_LastPeakSum(.)) %>% extract2(2) %>% rbind_all
rownames(Data4DESEQ) <- paste0(Data4DESEQ$chr,":",Data4DESEQ$start,"-",Data4DESEQ$end,":",Data4DESEQ$peak)

LE <- c(rep("NLE",length(condition)),rep("LE",length(condition)))
colData <- data.frame(t(rbind("condition"=as.vector(condition),LE)))
myData <- as.matrix(Data4DESEQ[,8:(ncol(Data4DESEQ))])
rownames(myData) <- rownames(Data4DESEQ)
dds <- DESeqDataSetFromMatrix(countData=myData, colData=colData, design = ~ LE + condition + LE:condition)
sizeFactors(dds) <- rep(sizeFactors(ddsNorm),2)
mydds <- DESeq(dds)
res <- results(mydds,contrast=c(0,0,0,1))

NLE_cond1 <- which(colData[,"condition"] == cond[1] & colData[,"LE"] == "NLE")
NLE_cond2 <- which(colData[,"condition"] == cond[2] & colData[,"LE"] == "NLE")
LE_cond1 <- which(colData[,"condition"] == cond[1] & colData[,"LE"] == "LE")
LE_cond2 <- which(colData[,"condition"] == cond[2] & colData[,"LE"] == "LE")
logFC_intro <- log2(rowMeans(counts(dds,normalized=TRUE)[,NLE_cond1])/rowMeans(counts(dds,normalized=TRUE)[,NLE_cond2]))
logFC_LE <- log2(rowMeans(counts(dds,normalized=TRUE)[,LE_cond1])/rowMeans(counts(dds,normalized=TRUE)[,LE_cond2]))

resO <- cbind(as.matrix(res),counts(dds,normalized=T), "logFC_intro"=logFC_intro,"logFC_LE"=logFC_LE)
colnames(resO)[6+c(1:ncol(myData))] <- colnames(myData)
outfile <- sub(".bed", "_peakIntron_SumLE.res", peakfile)
write.table(resO,file=outfile,col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)

scaled_rpkm <- t(scale(t(counts(dds,normalized=TRUE))))
scaled_rpkm_means <- cbind(rowMeans(scaled_rpkm[,NLE_cond1]),rowMeans(scaled_rpkm[,NLE_cond2]),rowMeans(scaled_rpkm[,LE_cond1]),rowMeans(scaled_rpkm[,LE_cond2]))

## Make plots

pdf(file=paste(dirname(peakfile),"/diffan_heatmap.pdf",sep=""))

par(font.lab=2, cex=.8)
plot(x=logFC_LE, y=logFC_intro,pch=20,main=paste("Intronic peak vs Last Exon",paste(cond[1],cond[2],sep="-"), "(FDR 5%)"),xlab="log2 FC Last Exon",ylab="log2 FC Intronic peak", frame=FALSE)
points(x=logFC_LE[which(res$padj<0.05 & logFC_intro-logFC_LE > 0)],y=logFC_intro[which(res$padj<0.05 & logFC_intro-logFC_LE > 0)],pch=20,col="red")
points(x=logFC_LE[which(res$padj<0.05 & logFC_intro-logFC_LE < 0)],y=logFC_intro[which(res$padj<0.05 & logFC_intro-logFC_LE < 0)],pch=20,col="blue")
up <- length(logFC_LE[which(res$padj<0.05 & logFC_intro-logFC_LE > 0)])
down <-  length(logFC_LE[which(res$padj<0.05 & logFC_intro-logFC_LE < 0)])
legend("topleft",legend=up,text.col="red",bty="n",cex=1.5)
legend("bottomright",legend=down,text.col="blue",bty="n",cex=1.5)

plot(x=logFC_LE, y=logFC_intro,pch=20,main=paste("Intronic peak vs Last Exon",paste(cond[1],cond[2],sep="-"), "(FDR 10%)"),xlab="log2 FC Last Exon",ylab="log2 FC Intronic peak", frame=FALSE)
points(x=logFC_LE[which(res$padj<0.1 & logFC_intro-logFC_LE > 0)],y=logFC_intro[which(res$padj<0.1 & logFC_intro-logFC_LE > 0)],pch=20,col="red")
points(x=logFC_LE[which(res$padj<0.1 & logFC_intro-logFC_LE < 0)],y=logFC_intro[which(res$padj<0.1 & logFC_intro-logFC_LE < 0)],pch=20,col="blue")
up <- length(logFC_LE[which(res$padj<0.1 & logFC_intro-logFC_LE > 0)])
down <-  length(logFC_LE[which(res$padj<0.1 & logFC_intro-logFC_LE < 0)])
legend("topleft",legend=up,text.col="red",bty="n",cex=1.5)
legend("bottomright",legend=down,text.col="blue",bty="n",cex=1.5)

plot(x=logFC_LE, y=logFC_intro,pch=20,main=paste("Intronic peak vs Last Exon",paste(cond[1],cond[2],sep="-"), "(pvalue 5%)"),xlab="log2 FC Last Exon",ylab="log2 FC Intronic peak", frame=FALSE)
points(x=logFC_LE[which(res$pvalue<0.05 & logFC_intro-logFC_LE > 0)],y=logFC_intro[which(res$pvalue<0.05 & logFC_intro-logFC_LE > 0)],pch=20,col="red")
points(x=logFC_LE[which(res$pvalue<0.05 & logFC_intro-logFC_LE < 0)],y=logFC_intro[which(res$pvalue<0.05 & logFC_intro-logFC_LE < 0)],pch=20,col="blue")
up <- length(logFC_LE[which(res$pvalue<0.05 & logFC_intro-logFC_LE > 0)])
down <-  length(logFC_LE[which(res$pvalue<0.05 & logFC_intro-logFC_LE < 0)])
legend("topleft",legend=up,text.col="red",bty="n",cex=1.5)
legend("bottomright",legend=down,text.col="blue",bty="n",cex=1.5)

dev.off()




#list1 <- read.table("/data/kdi_prod/.kdi/project_workspace_0/713/acl/03.00/ANALYSES/pip_nicolas/Liste_1.txt")
#list2 <- read.table("/data/kdi_prod/.kdi/project_workspace_0/713/acl/03.00/ANALYSES/pip_nicolas/Liste_2.txt")

#heatmap.2(scaled_rpkm_means, density.info="none",Colv=FALSE, labRow=NA, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Intronic peak vs Last Exon")

#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), list1$V1),]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Intronic peak vs Last Exon",cexRow=0.5,labRow=sapply(strsplit(rownames(extract),"\\|"),FUN=function(x) x[2]))
#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), list2$V1),]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Intronic peak vs Last Exon",cexRow=0.5,labRow=sapply(strsplit(rownames(extract),"\\|"),FUN=function(x) x[2]))
#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), rbind(list1,list2)$V1),]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Intronic peak vs Last Exon",cexRow=0.5,labRow=sapply(strsplit(rownames(extract),"\\|"),FUN=function(x) x[2]))




##### peakIntron vs Last peak 
#Data4DESEQ <- group_by(peaks, gene) %>% do(res=PeakIntron_LastPeak(.)) %>% extract2(2) %>% rbind_all
#rownames(Data4DESEQ) <- paste0(Data4DESEQ$chr,":",Data4DESEQ$start,"-",Data4DESEQ$end,":",Data4DESEQ$peak)
#colData <- data.frame(t(rbind("condition"=as.vector(condition),LE)))
#myData <- as.matrix(Data4DESEQ[,8:(ncol(Data4DESEQ))])
#rownames(myData) <- rownames(Data4DESEQ)
#dds <- DESeqDataSetFromMatrix(countData=myData, colData=colData,design = ~ LE + condition + LE:condition)
#sizeFactors(dds) <- rep(sizeFactors(ddsNorm),2)
#mydds <- DESeq(dds)
#res <- results(mydds,contrast=c(0,0,0,1))

#resO <- cbind(as.matrix(res),counts(dds,normalized=T))
#colnames(resO)[6+c(1:ncol(myData))] <- colnames(myData)
#outfile <- sub(".bed", "_peakIntron_LE.res", peakfile)
#write.table(resO,file=outfile,col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)

#scaled_rpkm <- t(scale(t(counts(dds,normalized=TRUE))))
#scaled_rpkm_means <- cbind(rowMeans(scaled_rpkm[,1:2]),rowMeans(scaled_rpkm[,3:4]),rowMeans(scaled_rpkm[,5:6]),rowMeans(scaled_rpkm[,7:8]))
#heatmap.2(scaled_rpkm_means, density.info="none",Colv=FALSE, labRow=NA, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE, sep="_")[c(1,3,5,7)],cexCol=0.7,main = "Intronic peak vs Last Exonic peak")

#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), list1$V1),3:4]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Intronic peak vs Last Exonic peak",cexRow=0.8,labRow=sapply(strsplit(rownames(extract),"_"),FUN=function(x) x[3]))
#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), list2$V1),3:4]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Intronic peak vs Last Exonic peak",cexRow=0.8,labRow=sapply(strsplit(rownames(extract),"_"),FUN=function(x) x[3]))
#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), rbind(list1,list2)$V1),]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Intronic peak vs Last Exonic peak",cexRow=0.5,labRow=sapply(strsplit(rownames(extract),"\\|"),FUN=function(x) x[2]))


#logFC_intro <- log2(counts(dds,normalized=TRUE)[,1:2]/counts(dds,normalized=TRUE)[,3:4])
#logFC_LE <- log2(counts(dds,normalized=TRUE)[,5:6]/counts(dds,normalized=TRUE)[,7:8])
#plot(x=logFC_LE, y=logFC_intro,pch=20,main="Intronic peak vs Last Exonic peak",xlab="log2 FC Last Exon",ylab="log2 FC Intronic peak")
#points(x=logFC_LE[which(res$padj<0.05),],y=logFC_intro[which(res$padj<0.05),],pch=20,col="red")


### peak Last Exon vs Last peak
#Data4DESEQ <- group_by(peaks, gene) %>% do(res=PeakLastExon_LastPeak(.)) %>% extract2(2) %>% rbind_all
#rownames(Data4DESEQ) <- paste0(Data4DESEQ$chr,":",Data4DESEQ$start,"-",Data4DESEQ$end,":",Data4DESEQ$peak)
#colData <- data.frame(t(rbind("condition"=as.vector(condition),LE)))
#myData <- as.matrix(Data4DESEQ[,8:(ncol(Data4DESEQ))])
#rownames(myData) <- rownames(Data4DESEQ)
#dds <- DESeqDataSetFromMatrix(countData=myData, colData=colData,design = ~ LE + condition + LE:condition)
#sizeFactors(dds) <- rep(sizeFactors(ddsNorm),2)
#mydds <- DESeq(dds)
#res <- results(mydds,contrast=c(0,0,0,1))

#resO <- cbind(as.matrix(res),counts(dds,normalized=T))
#colnames(resO)[6+c(1:ncol(myData))] <- colnames(myData)

#outfile <- sub(".bed", "_peakLE_LP.res", peakfile)
#write.table(resO,file=outfile,col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)

#scaled_rpkm <- t(scale(t(counts(dds,normalized=TRUE))))
#scaled_rpkm_means <- cbind(rowMeans(scaled_rpkm[,1:2]),rowMeans(scaled_rpkm[,3:4]),rowMeans(scaled_rpkm[,5:6]),rowMeans(scaled_rpkm[,7:8]))
#heatmap.2(scaled_rpkm_means, density.info="none",Colv=FALSE, labRow=NA, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE, sep="_")[c(1,3,5,7)],cexCol=0.7,main = "Exonic peak vs Last Exonic peak")

#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), list1$V1),]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Exonic peak vs Last Exonic peak",cexRow=0.5,labRow=sapply(strsplit(rownames(extract),"\\|"),FUN=function(x) x[2]))
#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), list2$V1),]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Exonic peak vs Last Exonic peak",cexRow=0.5,labRow=sapply(strsplit(rownames(extract),"\\|"),FUN=function(x) x[2]))
#extract <- scaled_rpkm_means[sapply(strsplit(rownames(scaled_rpkm_means),"_"),FUN=function(x, listG) any(x%in%listG), rbind(list1,list2)$V1),]
#heatmap.2(extract, density.info="none",Colv=FALSE, symkey=F, trace="none",col=colorpanel(40,low="blue",mid="white",high="red"), hclust=hclustW ,margins=c(6,6),labCol=paste(condition,LE,sep="_")[c(1,3,5,7)],cexCol=0.7, main = "Exonic peak vs Last Exonic peak",cexRow=0.5,labRow=sapply(strsplit(rownames(extract),"\\|"),FUN=function(x) x[2]))


#logFC_intro <- log2(counts(dds,normalized=TRUE)[,1:2]/counts(dds,normalized=TRUE)[,3:4])
#logFC_LE <- log2(counts(dds,normalized=TRUE)[,5:6]/counts(dds,normalized=TRUE)[,7:8])
#plot(x=logFC_LE, y=logFC_intro,pch=20,main="Exonic peak vs Last Exonic peak",xlab="log2 FC Last Exon",ylab="log2 FC Intronic peak")
#points(x=logFC_LE[which(res$padj<0.05),],y=logFC_intro[which(res$padj<0.05),],pch=20,col="red")

