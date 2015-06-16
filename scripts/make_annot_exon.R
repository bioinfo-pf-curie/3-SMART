#############################################################
############ Make file which contains all exons #############
#############################################################
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

library(plyr)

###############################################################################################################
###############################################################################################################


#Data download 
#Data provided of UCSC Refseq Genes table RefGene
f=read.table(txfile, sep="\t")
#Addition of a name for the 13th colonne
colnames(f)[13]="gene"

#Modification of the type of colonne in character
f[,10]=as.character(f[,10])
f[,11]=as.character(f[,11])

#Add NM or NR according of the type of the gene
liste<-sapply(strsplit(as.character(f[,2]),"_"),"[",1)
#Creation of a new name paste of chr name2 and strand
f[,13]<-paste(f[,3],f[,13],f[,4],liste,sep="_")

#Order the table with genes that have more that one transcript according to the gene names then transcripts start
f<-f[order(as.character(f[,13]),as.numeric(f[,5])),]#data.frame

genenames<-names(which(table(as.character(f[,"gene"]))==1))
f1<-f[which(f[,13]%in%genenames),]
f2<-f[which(!f[,13]%in%genenames),]

res=ddply(f2,.(gene),function(tab){diff=tab[2:dim(tab)[1],5]-tab[1:(dim(tab)[1]-1),6]; return(length(which(diff>0)))},.progress="text")

f3<-f2[which(f2[,13]%in%res[which(res[,2]>0),1]),]

res2=ddply(f3,.(gene),function(tab){
 diff=tab[2:dim(tab)[1],5]-tab[1:(dim(tab)[1]-1),6]
 diff0=c(0,which(diff>0),dim(tab)[1])
 j=0
 for(i in 1:(length(diff0)-1)){
   j=j+1
   tab[(diff0[i]+1):(diff0[i+1]),"gene"]=paste(tab[(diff0[i]+1):(diff0[i+1]),"gene"],j,sep="_")
 }
 return(tab)
},.progress="text")

alltx=rbind(f1,f2[which(!f2[,13]%in%res[which(res[,2]>0),1]),],res2)

alltx<-alltx[order(as.character(alltx[,13]),as.numeric(alltx[,5])),]#data.frame

###################################################################################################################
############################# Make table with all the exons on each transcript ####################################
####################################################################################################################

numberExon<-rep(alltx[,9],alltx[,9])
txstart<-unlist(strsplit(alltx[,10],","))
txend<-unlist(strsplit(alltx[,11],","))
nameTranscript<-rep(alltx[,2],alltx[,9])
gene=rep(alltx[,13],alltx[,9])
chr=rep(as.character(alltx[,3]),alltx[,9])
strand=rep(as.character(alltx[,4]),alltx[,9])
numE=paste("E",
          unlist(apply(alltx,1,function(tab){if(tab[4]=="+"){return(1:tab[9])}
                                             else{return(tab[9]:1)}})),sep="")
transcripts=cbind(chr,txstart,txend,paste(nameTranscript,gene,numberExon,numE,sep="|"),strand)

###########################################################################################################################################
############################################### Distribution if exon are the last exon of the gene ########################################
###########################################################################################################################################

transcripts<-as.data.frame(transcripts)
#List of last exon of the gene
le<-read.table(lefile)

#Under table update for the join
LE<-cbind(le[,c(1,2,3)],"1")
#Put colnames for the join fonction
colnames(transcripts)[c(1,2,3,5)]=c("chr","start","end","strand")
colnames(LE)[c(1,2,3)]=c("chr","start","end")
#Put LE if coordinates of the exon in  "LE" correspond to the coordinates of an exon in table transcripts 
transcripts<-join(transcripts,LE,by=c("chr","start","end"),match="first")
#Add a level to put IE (internal Exon) 
levels(transcripts[,6])<-c("1","0")
#Add the level for peaks which are not LE
transcripts[is.na(transcripts[,6]),6] <- rep("0",length(transcripts[is.na(transcripts[,6]),6]))
#Save the table

write.table(transcripts[,c(1,2,3,4,6,5)],paste(out,'/transcriptsType.bed',sep=""),quote=F,col.names=F,row.names=F,sep="\t")
#transcripts<-read.table("~/PolyAseq/tableTranscrits/transcriptsType.bed",sep="\t")
#colnames(transcripts)<-c("chr","start","end","name","strand","type")

