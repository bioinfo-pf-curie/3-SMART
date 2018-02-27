library(plyr)
library(GenomicRanges)

rm(list = ls())
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0) {
for (i in 1:la)
eval(parse(text = args[[i]]))
}


refseq <- read.table(infile, header = TRUE, comment.char = "~", check.names = FALSE, stringsAsFactors = FALSE)

# Distinction between NM and NR genes
refseq[grep("NM_", refseq[,"name"]), "name2"] <- paste(refseq[grep("NM_", refseq[,"name"]), "name2"], "NM", sep = "_")
refseq[grep("NR_", refseq[,"name"]), "name2"] <- paste(refseq[grep("NR_", refseq[,"name"]), "name2"], "NR", sep = "_")

refseq$gene <- paste(refseq[,"chrom"], refseq[,"name2"], refseq[,"strand"], sep = "_") # New gene name depending on the strand

# Genes with only one transcripts
gene1NM <- refseq[which(refseq$gene %in% names(which(table(refseq$gene) == 1))),]

# Genes with several transcripts
geneNMs <- refseq[which(refseq$gene %in% names(which(table(refseq$gene) > 1))),]
geneNMs <- geneNMs[order(geneNMs[,"name2"], geneNMs[,"txStart"]),]

# Differenciation between genes with overlapping transcripts (same gene) and genes without overlapping transcripts (transcripts are then considered different genes)

res <- ddply(geneNMs,.(gene), function(tab) {
	if (sum(diff(as.numeric(tab[,"txStart"]))) != 0) # if every transcripts have the same start, they overlap
		{
		if (sum(diff(as.numeric(tab[,"txEnd"]))) != 0) { # if every transcripts have the same end, they overlap
			tab_gr <- GRanges(seqnames = tab[,"chrom"], ranges = IRanges(start = tab[,"txStart"], end = tab[,"txEnd"]), strand = tab[,"strand"])
			overlap = which(countOverlaps(tab_gr) == 1)
			if (length(overlap) != 0) {tab[overlap,"gene"] = paste(tab[overlap,"gene"], 1:length(overlap), sep = "_")}
		}
	}	
	return(tab)
},.progress = "text")

refseq_polyA <- rbind(gene1NM, res)
write.table(refseq_polyA, paste(out_dir, "refseq_polyA.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Creation of the 3 files : whole gene (min 5' and max 3', can be chimeric), last exon (most 3', biggest one), whole gene without last exon 

pc = proc.time()
res_files <- lapply(unique(refseq_polyA[,"gene"]), function(gene) {
	
	tab = refseq_polyA[which(refseq_polyA[,"gene"] == gene),]

	exon_start = strsplit(tab[,"exonStarts"], ",")
	exon_end = strsplit(tab[,"exonEnds"], ",")

	exon1s = sapply(exon_start, "[",1) #exon1s_min=min(exon1s)
	exon1e = sapply(exon_end, "[",1)[which(exon1s == min(exon1s))] #exon1e_max=max(exon1e)
	chosen_tr1 = tab[which(exon1e == max(exon1e))[1], "name"]
	nb_exon_tr1 = tab[which(exon1e == max(exon1e))[1], "exonCount"]
	
	exonfe = sapply(exon_end,function(x){x[length(x)] }) 	#exonfe_max=max(exonfe)
 	exonfs = sapply(exon_start,function(x){x[length(x)]})[which(exonfe == max(exonfe))]	#exonfs_min=min(exonfs)
	chosen_trf = tab[which(exonfs == min(exonfs))[1], "name"] 
	nb_exon_trf = tab[which(exonfs == min(exonfs))[1], "exonCount"]
	
	all = tab[1, c("chrom", "txStart", "txEnd", "name2", "exonCount", "strand")] ;  
  	colnames(all) = c("chr", "start", "end", "name", "nb_tr", "strand") ; 
  	last_exon = all ; wo_le = c()

  	all[,"start"] = min(exon1s)  ; all[,"end"] = max(exonfe) ; all[,"nb_tr"] = nrow(tab)	
  	last_exon[,"nb_tr"] = nrow(tab)	

 	if (tab[1,"strand"] == "-") {
  		last_exon[,"start"] = min(exon1s) ; last_exon[,"end"] = max(exon1e) ; last_exon["name"] = paste(gsub("_NM|_NR", "", last_exon[,"name"]), chosen_tr1, sep = "_") 
		if (nb_exon_tr1 != 1) { wo_le = all ; wo_le[,"start"] = max(exon1e)  ; wo_le[,"end"] = max(exonfe) ; wo_le["name"] = paste(wo_le[,"name"], sep = "_") ; wo_le[,"nb_tr"] = nrow(tab)}
   	}else{
		last_exon[,"start"] = min(exonfs) ; last_exon[,"end"] = max(exonfe); last_exon["name"] = paste(gsub("_NM|_NR", "", last_exon[,"name"]), chosen_trf, sep = "_") 
 		if (nb_exon_trf != 1) { wo_le = all ; wo_le[,"start"] = min(exon1s)  ; wo_le[,"end"] = min(exonfs) ; wo_le["name"] = paste(wo_le[,"name"], sep = "_") ; wo_le[,"nb_tr"] = nrow(tab)}
   	}
	return(list(all,last_exon,wo_le))
})
proc.time() - pc

res_files_new <- matrix(unlist(sapply(res_files, "[",1)), ncol = 6, byrow = TRUE)
#change the name like : chr_start_end_name_NM_strand
name_res_files<- paste(res_files_new[,1],res_files_new[,2],res_files_new[,3],res_files_new[,4],res_files_new[,6], sep="_")
res_files_new_a <- data.frame(res_files_new[,1], res_files_new[,2], res_files_new[,3], name_res_files, res_files_new[,5], res_files_new[,6])
write.table(res_files_new_a, paste(out_dir, "whole_gene.bed", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

res_last_exon <- matrix(unlist(sapply(res_files, "[",2)), ncol = 6, byrow = TRUE)
res_last_exon_a <- data.frame(res_last_exon[,1], res_last_exon[,2], res_last_exon[,3], name_res_files, res_last_exon[,5], res_last_exon[,6])
write.table(res_last_exon_a, paste(out_dir, "last_exon_gene.bed", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

res_files_wgwole <- matrix(unlist(sapply(res_files, "[",3)), ncol = 6, byrow = TRUE)
name_res_wgwole<- paste(res_files_wgwole[,1],res_files_wgwole[,2],res_files_wgwole[,3],res_files_wgwole[,4],res_files_wgwole[,6], sep="_")
res_files_wgwole_a <- data.frame(res_files_wgwole[,1],res_files_wgwole[,2],res_files_wgwole[,3],name_res_wgwole,res_files_wgwole[,5], res_files_wgwole[,6])
write.table(res_files_wgwole_a, paste(out_dir, "whole_gene_wo_last_exon.bed", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

