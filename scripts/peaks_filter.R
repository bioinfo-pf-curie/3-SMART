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


message("Extract flanking regions [",wsizedown,"] ...")
fseq <- getFlankingRegions(rois, wdwn=as.numeric(wsizedown), wup=NA, genome)

## A stretch
message("Look for 'A' stretch in downstream regions [",nstretch," - ",mism,"]...")
rois.stretch <- containsStretch(fseq$dwn, stretch="A", slen=as.numeric(nstretch), mm=as.numeric(mism))
message("Discarding ",length(which(rois.stretch==TRUE)), " peaks")
message("Look for 'A' stretch in downstream regions [",nstretchcons," - ",mism,"]...")
rois.stretch2 <- containsStretch(fseq$dwn, stretch="A", slen=as.numeric(nstretchcons), mm=0)
message("Discarding ",length(which(rois.stretch2==TRUE)), " peaks")

pname <- rois[which(!unlist(rois.stretch) & !unlist(rois.stretch2))]$name

## Keep LE peaks
if (keep_le_peaks == 1){
  le.gr <- loadAnnotData(LEAnnotFile)
  is_le <- filterPeaksOnAnnotation(rois, le.gr)$name
  pname <- union(pname, is_le)
}

pois <- rois[pname]
peak2discard <- rois[setdiff(names(rois), pname)]
message("Discarded peaks=",length(peak2discard))

## Export list of peaks 
message("Export results ...")
outfile <- sub(".bed$", "_genomstretch.bed", peakfile)
export(peak2discard, format="bed", con=outfile)

outfile <- sub(".bed$", "_filt.bed", peakfile)
export(pois, format="bed", con=outfile)



