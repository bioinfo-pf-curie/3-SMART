rm(list = ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0) {
  for (i in 1:la)
    eval(parse(text = args[[i]]))
}


###################################################################
##
## Peak filtering 
##
###################################################################
source(polyA_lib)

require(RColorBrewer)
require(GenomicRanges)
require(rtracklayer)

## Import peaks as GRanges object
message("Import peaks file ...")
pois <- import(peakfile)
names(pois) <- pois$name

chrname <- seqlevels(pois)
if (is.element("chrMT", chrname)) {
  chrname[which(chrname == "chrMT")] <- "chrM"
  seqlevels(pois) <- chrname
}

message(length(pois)," loaded")

## Remove peaks which overlap several consecutive exons
message("Discard peaks which are not in last exons nor in intronic regions ...")

## ele
ele.gr <- loadAnnotData(ELEAnnotFile)

## le
le.gr <- loadAnnotData(LEAnnotFile)
is_le <- filterPeaksOnAnnotation(pois, le.gr, minoverlap = as.numeric(le_overlap))

## intron
exons.gr <- loadAnnotData(TRSAnnotFile)
is_intron <- filterPeaksOnAnnotation(pois, exons.gr, invert = TRUE,  minoverlap = as.numeric(intron_overlap))

finalpeaks <- unique(c(is_le$name, is_intron$name))
pois.filt <- pois[finalpeaks]

message("Discarded peaks=", length(pois) - length(pois.filt))
message("Peaks of interest=", length(pois.filt))

## Export list of peaks
message("Export results ...")
outfile <- sub(".bed$", "_notlastexons_notintron.bed", peakfile)
export(pois[setdiff(pois$name, pois.filt$name)], con = outfile)

outfile <- sub(".bed$", "_finallist.bed", peakfile)
export(pois.filt, format = "bed", con = outfile)

## Load annotation
message("Annotate Peaks ...")
message("Minimum Overlap = ",minover)

pois.le.nr <- annotatePeaks(pois.filt, le.gr[grep("NR", le.gr$nm)], outfile = sub(".bed$", "_NR_LE_annot.bed", peakfile), minover = as.numeric(minover))
pois.le.nm <- annotatePeaks(pois.filt, le.gr[grep("NM", le.gr$nm)],  outfile = sub(".bed$", "_NM_LE_annot.bed", peakfile), minover = as.numeric(minover))
pois.ele.nr <- annotatePeaks(pois.filt, ele.gr[grep("NR", ele.gr$nm)], outfile = sub(".bed$", "_NR_ELE_annot.bed", peakfile), minover = as.numeric(minover))
pois.ele.nm <- annotatePeaks(pois.filt, ele.gr[grep("NM", ele.gr$nm)], outfile = sub(".bed$", "_NM_ELE_annot.bed", peakfile), minover = as.numeric(minover))


## Combine Annotation
ov <- as.list(findOverlaps(subject = pois.le.nm, query = pois.filt, type = "equal", select = "all"))
x.index <- which(sapply(ov, length) != 0)
pois.filt$name = ""
pois.filt$status = ""
pois.filt$name[x.index] <- pois.le.nm$name
pois.filt$status[x.index] <- "LE"

ov <- as.list(findOverlaps(subject = pois.le.nr, query = pois.filt, type = "equal", select = "all"))
x.index <- which(sapply(ov, length) != 0)
pois.filt$name[x.index] <- paste(pois.le.nr$name, pois.filt$name[x.index], sep = "|")
pois.filt$status[x.index] <- paste("LE", pois.filt$status[x.index], sep = "|")

ov <- as.list(findOverlaps(subject = pois.ele.nm, query = pois.filt, type = "equal", select = "all"))
x.index <- which(sapply(ov, length) != 0)
pois.filt$name[x.index] <- paste0(pois.ele.nm$name, pois.filt$name[x.index], sep = "|")
pois.filt$status[x.index] <- paste("ELE", pois.filt$status[x.index], sep = "|")

ov <- as.list(findOverlaps(subject = pois.ele.nr, query = pois.filt, type = "equal", select = "all"))
x.index <- which(sapply(ov, length) != 0)
pois.filt$name[x.index] <- paste0(pois.ele.nr$name, pois.filt$name[x.index], sep = "|")
pois.filt$status[x.index] <- paste("ELE", pois.filt$status[x.index], sep = "|")

outfile = sub(".bed$", "_finallist_annot.tsv", peakfile)
write.table(as.data.frame(pois.filt), file = outfile, quote = FALSE, row.names = FALSE, sep = "\t")

