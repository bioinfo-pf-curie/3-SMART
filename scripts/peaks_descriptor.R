rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}


###################################################################
##
## Peak filtering - PolyA-seq pipeline
##
###################################################################
require(RColorBrewer)
require(GenomicRanges)
require(rtracklayer)
require(ggplot2)
source("scripts/polyA_lib.R")



## Load Genome
message("Load genome info ...")
if (org=="mm9"){
  genomePack <- "BSgenome.Mmusculus.UCSC.mm9"
}else if (org=="hg19"){
  genomePack <- "BSgenome.Hsapiens.UCSC.hg19"
}
stopifnot(require(genomePack, character.only = TRUE))
genome <- eval(as.name(genomePack))

## Import peaks as GRanges object
message("Import peaks file ...")
rois <- import(peakfile)
names(rois) <- rois$name
chrname <- seqlevels(rois)
if (is.element("chrMT",  chrname)){
  chrname[which(chrname=="chrMT")] <- "chrM"
  seqlevels(rois) <- chrname
}
seqinfo(rois) <- seqinfo(genome)[seqlevels(rois)]
message(length(rois)," loaded")

## Describe peaks
message("Extract flanking regions [",wsizeup," - ",wsizedown,"] ...")
fseq <- getFlankingRegions(rois, wdwn=as.numeric(wsizedown), wup=as.numeric(wsizeup), genome)

###################################################################
##
## Peak motifs 
##
###################################################################

## PolyA site
rois.fseq=list()
rois.fseq$peak=fseq$peak[names(rois)]
rois.fseq$up=fseq$up[names(rois)]
rois.fseq$dw=fseq$dw[names(rois)]

polyA <- read.csv(polyAfile, header=FALSE)
polyA <- as.list(as.character(polyA[,1]))
names(polyA) <- unlist(polyA)


## PolyA motifs
message("Look for known polyA site ...")
pap <- containsPolyAsignal(rois.fseq$peak, polyA)
rownames(pap) <- names(rois.fseq$peak)
paup <- containsPolyAsignal(rois.fseq$up, polyA)
rownames(paup) <- names(rois.fseq$up)
padw <- containsPolyAsignal(rois.fseq$dw, polyA)
rownames(padw) <- names(rois.fseq$dw)

outfile <- sub(".bed$", "_peak_polyAsites.csv", peakfile)
write.csv(pap, file=outfile, quote=FALSE)
outfile <- sub(".bed$", "_up_polyAsites.csv", peakfile)
write.csv(paup, file=outfile, quote=FALSE)
outfile <- sub(".bed$", "_dw_polyAsites.csv", peakfile)
write.csv(padw, file=outfile, quote=FALSE)                  

###################################################################
##
## Peak plots 
##
###################################################################

## Nb peaks with motifs

np <- table(rowSums(pap))
nu <- table(rowSums(paup))
nd <- table(rowSums(padw))

mnp <- np[as.character(0:10)]
names(mnp) <- as.character(0:10)
mnp["10"] <- mnp["10"] + sum(np[which(as.numeric(rownames(np))>10)])

mnu <- nu[as.character(0:10)]
names(mnu) <- as.character(0:10)
mnu["10"] <- mnu["10"] + sum(nu[which(as.numeric(rownames(nu))>10)])

mnd <- nd[as.character(0:10)]
names(mnd) <- as.character(0:10)
mnd["10"] <- mnd["10"] + sum(nd[which(as.numeric(rownames(nd))>10)])

z <- data.frame(nb=c(0:10,0:10,0:10), counts=c(mnp, mnd, mnu), win=c(rep("Peaks", 11), rep("Upstream", 11), rep("Downstream", 11)))

outfile <- sub(".bed$", "_polyA_ppeak.pdf", peakfile)
pdf(file=outfile)
ggplot(z, aes(x=nb, y=counts, fill=win)) + geom_bar(stat="identity", position="dodge") + ylab("Peaks Number") + xlab("Motifs Per Peak") + scale_fill_manual(values=c("grey77", "firebrick4","lightsteelblue3"), name="")+
  theme(axis.text=element_text(size=9), axis.title=element_text(face="bold", size=10), legend.position="bottom")
dev.off()

## ------------------------- ##

pmotif <- colSums(pap)
upmotif <- colSums(paup)
dwnmotif <- colSums(padw)
z <- data.frame(motif=names(upmotif), peak=pmotif/sum(pmotif)*100, up=upmotif/sum(upmotif)*100, dwn=dwnmotif/sum(upmotif)*100)
mypal <- colorRampPalette( brewer.pal( 9 , "RdBu" ) )
p1 <- ggplot(z, aes(x = motif, y=peak, fill=motif)) + geom_bar(width = 1, stat="identity") + coord_polar(theta="x") + xlab("") + ylab("") +
  scale_fill_manual(values=mypal(length(levels(z$motif))), name="") +  theme(axis.text=element_text(size=9), axis.title=element_text(face="bold", size=10), legend.text=element_text(size=6))
p2 <- ggplot(z, aes(x = motif, y=up, fill=motif)) + geom_bar(width = 1, stat="identity") + coord_polar(theta="x") + xlab("") + ylab("") +
  scale_fill_manual(values=mypal(length(levels(z$motif))), name="") +  theme(axis.text=element_text(size=9), axis.title=element_text(face="bold", size=10), legend.text=element_text(size=6))
p3 <- ggplot(z, aes(x = motif, y=dwn, fill=motif)) + geom_bar(width = 1, stat="identity") + coord_polar(theta="x") + xlab("") + ylab("") +
  scale_fill_manual(values=mypal(length(levels(z$motif))), name="") +  theme(axis.text=element_text(size=9), axis.title=element_text(face="bold", size=10), legend.text=element_text(size=6))

outfile <- sub(".bed$", "_peak_polyAsites.pdf", peakfile)
pdf(file=outfile)
plot(p1)
dev.off()

outfile <- sub(".bed$", "_up_polyAsites.pdf", peakfile)
pdf(file=outfile)
plot(p2)
dev.off()

outfile <- sub(".bed$", "dwn_polyAsites.pdf", peakfile)
pdf(file=outfile)
plot(p3)
dev.off()

