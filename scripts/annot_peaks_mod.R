rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

## annotatePeaks()
## Annotate gr by overlapping with annot
##
##
annotatePeaks <- function(x, annot, outfile, minover=1){
  ov <- as.list(findOverlaps(x, annot, minoverlap=minover))
  x.index <- which(sapply(ov, length)!=0)
  x.annot <- sapply(ov[x.index], function(idx){paste0(unique(annot$name[idx]), collapse="|")})
  x$name[x.index] <- paste0(x$name[x.index],"|",x.annot)

 # export(x[x.index], con=outfile)
  
  x[x.index]
}



## LoadAnnotData()
## Load exon file and annotate them
##
## con = input file
##
loadAnnotData <- function(con, random=FALSE){
  require(rtracklayer)
  message("Loading annotation file '", con,"' ...")
  
  ## Load annot data
  annot.gr <- import(con)
  annot.gr$nm <- sapply(strsplit(annot.gr$name,"|", fixed=TRUE),"[",1)
  annot.gr$symbol <- sapply(strsplit(annot.gr$name,"|", fixed=TRUE),"[",2)
  
  chrname <- seqlevels(annot.gr)
  if (is.element("chrMT",  chrname)){
    chrname[which(chrname=="chrMT")] <- "chrM"
    seqlevels(annot.gr) <- chrname
  }

  if (! random){
    annot.gr <- annot.gr[grep("random",seqnames(annot.gr), invert=TRUE)]
    seqlevels(annot.gr) <- seqlevels(annot.gr)[grep("random",seqlevels(annot.gr), invert=TRUE)]
  }
  annot.gr
}

## filterPeaksOnLastExon
## x = Peaks set to filter out
## annot = Genomic ranges of exon description
## invert = invert annotation, i.e intronic regions
##
filterPeaksOnAnnotation <- function(x, annot, invert=FALSE, ...){##, extend=5){

  if (invert){
    message("Invert Annotation ...")
    grl <- split(annot, annot$symbol)
    geneRanges <- range(grl)
    annot <- unlist(psetdiff(geneRanges, grl))
  }
  
  ## Select peaks on annotation
  subsetByOverlaps(x, annot, ...)
}


###################################################################
##
## Peak filtering 
##
###################################################################

require(RColorBrewer)
require(GenomicRanges)
require(rtracklayer)

## Import peaks as GRanges object
message("Import peaks file ...")
pois <- import(peakfile)
names(pois) <- pois$name


chrname <- seqlevels(pois)
if (is.element("chrMT",  chrname)){
  chrname[which(chrname=="chrMT")] <- "chrM"
  seqlevels(pois) <- chrname
}

message(length(pois)," loaded")

## Remove peaks which overlap several consecutive exons
message("Discard peaks which are not in last exons nor in intronic regions ...")

## ele
ele.gr <- loadAnnotData(ELEAnnotFile)

## le
le.gr <- loadAnnotData(LEAnnotFile)
is_le <- filterPeaksOnAnnotation(pois, le.gr, minoverlap=as.numeric(le_overlap))

## intron
exons.gr <- loadAnnotData(TRSAnnotFile)
is_intron<- filterPeaksOnAnnotation(pois, exons.gr, invert=TRUE,  minoverlap=as.numeric(intron_overlap))

#finalpeaks <- unique(c(is_le$name, is_intron$name))
#pois.filt <- pois[finalpeaks]

#message("Discarded peaks=", length(pois)-length(pois.filt))
#message("Peaks of interest=",length(pois.filt))

## Export list of peaks
#message("Export results ...")
#outfile <- sub(".bed$", "_notlastexons_notintron.bed", peakfile)
#export(pois[setdiff(pois$name, pois.filt$name)], con=outfile)

#outfile <- sub(".bed$", "_finallist.bed", peakfile)
#export(pois.filt, format="bed", con=outfile)


###################################################################
##
## Peak Annotation 
##
###################################################################

## ## PolyA site
## pois.filt.fseq=list()
## pois.filt.fseq$up=fseq$up[names(pois.filt)]
## pois.filt.fseq$dw=fseq$dw[names(pois.filt)]

## polyA <- read.csv(polyAfile, header=FALSE)
## polyA <- as.list(as.character(polyA[,1]))
## names(polyA) <- unlist(polyA)


## ## PolyA motifs
## message("Look for known polyA site ...")
## paup <- containsPolyAsignal(pois.filt.fseq$up, polyA)
## rownames(paup) <- names(pois.filt.fseq$up)
## padw <- containsPolyAsignal(pois.filt.fseq$dw, polyA)
## rownames(padw) <- names(pois.filt.fseq$dw)

## print(colSums(paup))
## print(colSums(padw))

## outfile <- sub(".bed$", "_up_polyAsites.csv", peakfile)
## write.csv(paup, file=outfile, quote=FALSE)
## outfile <- sub(".bed$", "_dw_polyAsites.csv", peakfile)
## write.csv(padw, file=outfile, quote=FALSE)                  


## Load annotation
message("Annotate Peaks ...")
message("Minimum Overlap = ",minover)

pois.le.nr <- annotatePeaks(pois, le.gr[grep("NR",le.gr$nm)], outfile=sub(".bed$", "_NR_LE_annot.bed", peakfile), minover=as.numeric(minover))
pois.le.nm <- annotatePeaks(pois, le.gr[grep("NM",le.gr$nm)],  outfile=sub(".bed$", "_NM_LE_annot.bed", peakfile), minover=as.numeric(minover))
pois.ele.nr <- annotatePeaks(pois, ele.gr[grep("NR",ele.gr$nm)], outfile=sub(".bed$", "_NR_ELE_annot.bed", peakfile), minover=as.numeric(minover))
pois.ele.nm <- annotatePeaks(pois, ele.gr[grep("NM",ele.gr$nm)], outfile=sub(".bed$", "_NM_ELE_annot.bed", peakfile), minover=as.numeric(minover))



ov <- as.list(findOverlaps(subject=pois.le.nm,query=pois,type="equal",select="all"))
x.index <- which(sapply(ov, length)!=0)
pois$name = ""
pois$status = ""
pois$name[x.index] <-   pois.le.nm$name
pois$status[x.index] <- "LE"

ov <- as.list(findOverlaps(subject=pois.le.nr,query=pois,type="equal",select="all"))
x.index <- which(sapply(ov, length)!=0)
pois$name[x.index] <- paste(pois.le.nr$name, pois$name[x.index],sep="|")
pois$status[x.index] <- paste("LE", pois$status[x.index],sep="|")

ov <- as.list(findOverlaps(subject=pois.ele.nm,query=pois,type="equal",select="all"))
x.index <- which(sapply(ov, length)!=0)
pois$name[x.index] <- paste0(pois.ele.nm$name, pois$name[x.index],sep="|")
pois$status[x.index] <- paste("ELE", pois$status[x.index],sep="|")

ov <- as.list(findOverlaps(subject=pois.ele.nr,query=pois,type="equal",select="all"))
x.index <- which(sapply(ov, length)!=0)
pois$name[x.index] <- paste0(pois.ele.nr$name, pois$name[x.index],sep="|")
pois$status[x.index] <- paste("ELE", pois$status[x.index],sep="|")

outfile=sub(".bed$", "_annot", peakfile)
write.table(as.data.frame(pois), file=outfile,quote=FALSE,row.names=FALSE, sep="\t")
