## getFlankingRegions()
## Extract flanking region upstream/downstream a peak
## Return a list of upstream and downstream sequences
##
## x = GRange object describing a peak
## wdwn = downstream window size
## wup = upstream window size
## updwseq = upstream-downstream window -> Keep a new window size [wup;wdwn] at the end of peak
## genome = BSgenome object
##

getFlankingRegions <- function(x, wdwn = NA, wup = NA, genome) {

  chrname <- seqlevels(x)
  if (is.element("chrMT",  chrname)) {
    chrname[which(chrname == "chrMT")] <- "chrM"
    seqlevels(x) <- chrname
  }

  ## extract peaks regions
  peakseq <- getSeq(genome, x)
  names(peakseq) <- elementMetadata(x)$name

  ## flanking regions are based from end of peak
  x <- resize(x, width = 1, fix = "end")
  useq <- dseq <- NULL
  
  if (!is.na(wdwn)) {
    dwgr <- trim(flank(x, width = as.numeric(wdwn), start = FALSE))
    dseq <- getSeq(genome, names = dwgr)
    names(dseq) <- elementMetadata(x)$name
  }
  
  if (!is.na(wup)) {
    upgr <- trim(flank(x, width = as.numeric(wup), start = TRUE))
    useq <- getSeq(genome, names = upgr)
    names(useq) <- elementMetadata(x)$name
  }

  if (!is.na(wdwn) & !is.na(wup)) {
    dwgr <- trim(flank(x, width = as.numeric(wdwn), start = FALSE))
    updwpgr <- resize(dwgr, width = 1, fix = "end")
    updwpgrT <- trim(flank(updwpgr, width = as.numeric(wdwn + wup), start = TRUE))
    updwseq <- getSeq(genome, names = updwpgrT)
    names(updwseq) <- elementMetadata(x)$name
  }
 
  ##
  ## Extract sequences and reverse complement if necessary
  if (is.na(wdwn) | is.na(wup)) {
    return(list(peak = peakseq, up = useq, dwn = dseq))
  }else{
    return(list(peak = peakseq, up = useq, dwn = dseq, updw = updwseq))
  }
}

## containsStretch()
## Look for stretch in sequence
## Return a boolean, True if a strech is found
##
## x = DNAstring sequence
## stretch = The sequence to be stretched
## slen = length of stretch
## mm = Mismatch number in stretch
##
containsStretch <- function(x, stretch, slen, mm) {
  cnt <- vcountPattern(paste0(rep(stretch, slen), collapse = ""), x, max.mismatch = mm)
  ifelse(cnt > 0, TRUE, FALSE)
}

## containsStretch()
## Look for polyA signal in sequence
## Return a vector with the number of times each motif is found
##
## x = DNAstring sequence
## motifs = a list of motif to search
##
containsPolyAsignal <- function(x, motifs) {
  msearch <- sapply(motifs, function(m) {
     vcountPattern(m, x)
  })
  names(msearch) <- names(motifs)
  msearch
}


## annotatePeaks()
## Annotate gr by overlapping with annot
##
##
annotatePeaks <- function(x, annot, outfile, minover = 1) {
  ov <- as.list(findOverlaps(x, annot, minoverlap = minover))
  x.index <- which(sapply(ov, length) != 0)
  x.annot <- sapply(ov[x.index], function(idx){paste0(unique(annot$name[idx]), collapse = "|")})
  x$name[x.index] <- paste0(x$name[x.index],"|",x.annot)

  export(x[x.index], con = outfile)
  
  x[x.index]
}

## LoadAnnotData()
## Load exon file and annotate them
##
## con = input file
##
loadAnnotData <- function(con, random = FALSE) {
  require(rtracklayer)
  message("Loading annotation file '", con,"' ...")
  
  ## Load annot data
  annot.gr <- import(con)
  annot.gr$nm <- sapply(strsplit(annot.gr$name, "|", fixed = TRUE),"[",1)
  annot.gr$symbol <- sapply(strsplit(annot.gr$name, "|", fixed = TRUE),"[",2)
  
  chrname <- seqlevels(annot.gr)
  if (is.element("chrMT",  chrname)) {
    chrname[which(chrname == "chrMT")] <- "chrM"
    seqlevels(annot.gr) <- chrname
  }

  if (!random) {
    annot.gr <- annot.gr[grep("random", seqnames(annot.gr), invert = TRUE)]
    seqlevels(annot.gr) <- seqlevels(annot.gr)[grep("random",seqlevels(annot.gr), invert = TRUE)]
  }
  annot.gr
}

## filterPeaksOnLastExon
## x = Peaks set to filter out
## annot = Genomic ranges of exon description
## invert = invert annotation, i.e intronic regions
##
filterPeaksOnAnnotation <- function(x, annot, invert = FALSE, ...){##, extend=5){

  if (invert) {
    message("Invert Annotation ...")
    grl <- split(annot, annot$symbol)
    geneRanges <- range(grl)
    annot <- unlist(psetdiff(geneRanges, grl))
  }
  
  ## Select peaks on annotation
  subsetByOverlaps(x, annot, ...)
}


PeakIntron_LastPeak <- function(x){
  y <- x[grep("^LE|\\|LE", x$status),]
  if (nrow(y) > 0) {
    if (x[1,"strand"] == "+") {
      idLast <- y[which(y$end == max(y$end)),]
    } else{idLast <- y[which(y$start == min(y$start)),]} 
    NLE <- x[-which(rownames(idLast) == rownames(x)), ]
    NLE <- NLE[grep("ELE",NLE$status),c(1:3,6:(ncol(x) - 1))]
    if (nrow(NLE) > 0) {
      LE <- x[which(rownames(idLast) == rownames(x)), c(10:(ncol(x) - 1)) ]
      colnames(LE) <- paste0("LE.", colnames(LE))
      return(cbind(NLE, LE))
    }else return(NULL)
  } else return(NULL)
}

PeakIntron_LastPeakSum <- function(x) {
    y <- x[grep("^LE|\\|LE", x$status),]
      if (nrow(y) > 0) {
                NLE <- x[grep("^LE|\\|LE", x$status, invert = TRUE), ]
                NLE <- NLE[grep("ELE", NLE$status), c(1:3,6:(ncol(x) - 1))]
          if (nrow(NLE) > 0) {
                LE <- colSums(y[,c(10:(ncol(x) - 1)) ])
                names(LE) <- paste0("LE.", names(LE))
                return(cbind(NLE, as.list(LE)))
              }else return(NULL)
          } else return(NULL)
}

PeakLastExon_LastPeak <- function(x) {
      y <- x[grep("^LE|\\|LE", x$status),]
            if ( nrow(y) > 0) {
              if (x[1, "strand"] == "+") {
                      idLast <- y[which(y$end == max(y$end)),]
              }else{idLast <- y[which(y$start == min(y$start)),]}
              NLE <- y[-which(rownames(idLast) == rownames(y)), c(1:3,6:(ncol(x) - 1)) ]
if (nrow(NLE) > 0) {
              LE <- x[which(rownames(idLast) == rownames(y)),c(10:(ncol(x) - 1)) ]
              colnames(LE) <- paste0("LE.", colnames(LE))
              return(cbind(NLE, LE))
            }else return(NULL)
             }else return(NULL)
    }

hclustW <- function(x, ...) hclust(x, method = "ward", ...)
