#!/bin/bash

## Init
dir=$(dirname $0)

NORMAL="\\033[0;39m" 
RED="\\033[1;31m"

die() {
    echo -e "$RED""$*""$NORMAL" 1>&2
    exit 1
}

## Get args
usage() {
    echo "usage: $0 -c CONFIG -i INPUT_FILE -o OUTPUT_FILE"
}

if [ $# -eq 0 ]; then
    usage
    exit
fi

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) CONF=$2; shift;;
	(-i) INPUT=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) die "$0: error - unrecognized option $1" 1>&2;;
	(*)  break;;
    esac
    shift
done

if [[ -z ${CONF} || -z ${INPUT} || -z ${OUTPUT} ]]; then
    usage
    exit
fi

## Read configuration files
source ${CONF}

## Ouptut file
if [ ! -d ${OUTPUT} ]; then
    mkdir -p ${OUTPUT}
fi

## ANNOTATION
if [ $BUILD_ANNOT == 1 ]
then
    echo "Build annotation files ..."             
    ## Create gene file
    echo "${R_PATH}/R --vanilla CMD BATCH \"--args txFile='${ANNOT_DIR}/${ORG}/${UCSC_EXPORT}' out='${ANNOT_DIR}/${ORG}'\" ${SCRIPTS}/make_annot_gene.R ${LOGS}/make_annot_gene.Rout"
    ${R_PATH}/R --vanilla CMD BATCH "--args txfile='${ANNOT_DIR}/${ORG}/${UCSC_EXPORT}' out='${ANNOT_DIR}/${ORG}'" ${SCRIPTS}/make_annot_gene.R ${LOGS}/make_annot_gene.Rout
    ## Create exon file
    echo "${R_PATH}/R --vanilla CMD BATCH \"--args txfile='${ANNOT_DIR}/${ORG}/${UCSC_EXPORT}' lefile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' out='${ANNOT_DIR}/${ORG}'\" ${SCRIPTS}/make_annot_exon.R ${LOGS}/make_annot_exon.Rout"
    ${R_PATH}/R --vanilla CMD BATCH "--args txfile='${ANNOT_DIR}/${ORG}/${UCSC_EXPORT}' lefile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' out='${ANNOT_DIR}/${ORG}'" ${SCRIPTS}/make_annot_exon.R ${LOGS}/make_annot_exon.Rout
fi 

if [[ ! -e ${ANNOT_DIR}/${ORG}/last_exon_gene.bed || ! -e ${ANNOT_DIR}/${ORG}/whole_gene.bed || ! -e ${ANNOT_DIR}/${ORG}/whole_gene_wo_last_exon.bed || ! -e ${ANNOT_DIR}/${ORG}/transcriptsType.bed ]]
then
    die "Error in annotation process"
fi

## Quality Control
echo "Filter low quality reads (MAPQ < ${MIN_MAPQ}) ..."
BAM_PX=${OUTPUT}/`basename ${INPUT} | sed -e "s/.bam$/_q${MIN_MAPQ}/"`
BAMQ=${BAM_PX}.bam

#${SAMTOOLS_PATH}/samtools view -b -q ${MIN_MAPQ} ${INPUT} | ${SAMTOOLS_PATH}/samtools sort -  ${BAM_PX} # @ option for parallelism ?
#${SAMTOOLS_PATH}/samtools index ${BAMQ}

# Add some statistics about the quality filtering (half of the library removed, full of AAA..)

## Detect peaks
echo "Peak Detection ..."
#BEDQ=`echo ${BAMQ} | sed -e "s/.bam$/.bed/"`
BED_PEAKS=`echo ${BAMQ} | sed -e "s/.bam$/_peaks.bed/"`
#${BEDTOOLS_PATH}/bedtools bamtobed -i  ${BAMQ} > ${BEDQ}
#${BEDTOOLS_PATH}/bedtools merge -s -n -d ${MAX_DIST_MERGE} -i ${BEDQ} | ${AWK_PATH}/awk -v thr=${MIN_NB_READS_PER_PEAK} 'BEGIN{OFS="\t";c=1}($4>=thr){print $1,$2,$3,"peak_"c,$4,$5;c=c+1}' > ${BED_PEAKS}
#${BEDTOOLS_PATH}/bedtools bamtobed -i  ${BAMQ} | ${BEDTOOLS_PATH}/bedtools merge -s -n -d ${MAX_DIST_MERGE} -i stdin | ${AWK_PATH}/awk -v thr=${MIN_NB_READS_PER_PEAK} 'BEGIN{OFS="\t";c=1}($4>=thr){print $1,$2,$3,"peak_"c,$4,$5;c=c+1}' > ${BED_PEAKS}


## Add statistics on peak size /reads distribution 
## TODO

## Filter peaks - stretchA + notLE/intro
echo "Filter peaks ..."
#echo "${R_PATH}/R --vanilla CMD BATCH \"--args peakfile='${BED_PEAKS}' LEAnnotFile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' polyAfile='${ANNOT_DIR}/polyA_signal.csv' wsizeup='${WINSIZE_UP}' wsizedown='${WINSIZE_DOWN}' nstretch='${NB_STRETCH_POLYA}' mism='${MISM}' nstretchcons='${NB_STRETCH_CONSECUTIVE}' org='${ORG}'\" ${SCRIPTS}/peaks_filter.R ${LOGS}/peaks_filter.Rout"
#${R_PATH}/R --vanilla CMD BATCH "--args peakfile='${BED_PEAKS}' LEAnnotFile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' TRSAnnotFile='${ANNOT_DIR}/${ORG}/transcriptsType.bed' polyAfile='${ANNOT_DIR}/polyA_signal.csv' wsizeup='${WINSIZE_UP}' wsizedown='${WINSIZE_DOWN}' nstretch='${NB_STRETCH_POLYA}' mism='${MISM}' nstretchcons='${NB_STRETCH_CONSECUTIVE}' org='${ORG}'" ${SCRIPTS}/peaks_filter.R ${LOGS}/peaks_filter.Rout

FILT_PEAKS=`echo ${BED_PEAKS} | sed -e 's/.bed$/_filt.bed/'`
if [[ ! -e ${FILT_PEAKS} ]]; then
    die "Error : Filtered peaks file not found"
fi

## Annotate peaks - Annotation based on Gene level
echo "Annotate peaks ..."
echo "${R_PATH}/R --vanilla CMD BATCH \"--args peakfile='${FILT_PEAKS}' annotDir='${ANNOT_DIR}/${ORG}' minover='${MIN_OVERLAP}'\" ${SCRIPTS}/annot_peaks.R ${LOGS}/annot_peaks.Rout"
<<<<<<< .mine
${R_PATH}/R --vanilla CMD BATCH "--args peakfile='${FILT_PEAKS}' annotDir='${ANNOT_DIR}/${ORG}' minover='${MIN_OVERLAP}'" ${SCRIPTS}/annot_peaks.R ${LOGS}/annot_peaks.Rout
=======
${R_PATH}/R --vanilla CMD BATCH "--args peakfile='${FILT_PEAKS}'  minover='${MIN_OVERLAP}' annotDir='${ANNOT_DIR}/${ORG}'" ${SCRIPTS}/annot_peaks.R ${LOGS}/annot_peaks.Rout
>>>>>>> .r26185

ANNO_PEAKS=`echo ${FILT_PEAKS} | sed -e 's/.bed$/_NM_ELE_annot.bed/'`
if [[ ! -e ${ANNO_PEAKS} ]]; then
    die "Error : Annotated peaks file not found"
fi

echo "done !!"

# 	echo "Création des sites des potentiels sites de polyAdénylation pour la comparaison U1Amo/4T1"
# 	#File with U1 rep1, U1 rep2, 4T1 rep1 and 4T1 rep2 samples
# 	cat $fileout/'A162T7/A162T7peaksfiltes.bed' $fileout/'A162T8/A162T8peaksfiltes.bed' $fileout/'A162T9/A162T9peaksfiltes.bed' $fileout/'A162T12/A162T12peaksfiltes.bed' > $fileout/'4T1_U1.bed'
# 	$bedtools sort -i $fileout/'4T1_U1.bed' > $fileout/'4T1_U1.sorted.bed'
# 	$bedtools merge -s -n -i $fileout/'4T1_U1.sorted.bed' > $fileout/'4T1_U1_peaks.bed'
# 	#File with U1 rep1, U1 rep2, 4T1 rep1 and 4T1 rep2 samples
# 	cat $fileout/'A162T7/A162T7peaksfiltes.bed' $fileout/'A162T8/A162T8peaksfiltes.bed' $fileout/'A162T1/A162T1peaksfiltes.bed' $fileout/'A162T2/A162T2peaksfiltes.bed' > $fileout/'4T1_67NR.bed'
# 	$bedtools sort -i $fileout/'4T1_67NR.bed' > $fileout/'4T1_67NR.sorted.bed'
# 	$bedtools merge -s -n -i $fileout/'4T1_67NR.sorted.bed' > $fileout/'4T1_67NR_peaks.bed'
	
# 	$R --vanilla -q --args $fileout/ 4T1_U1_ < $me/addname_1.R	
# 	$R --vanilla -q --args $fileout/ 4T1_67NR_ < $me/addname_1.R

# 	echo "Annotation des sites des potentiels sites de polyAdénylation"
# 	##Find peaks corresponding to genes
# 	#For NM_genes
# 	$bedtools intersect -a $fileout/'genesTables/all_except_last_exon_NM.bed' -b $fileout/'4T1_U1_peaksname.bed' -loj -s > $fileout/'e_le_4T1_U1_NM.bed'
# 	$bedtools intersect -a $fileout/'genesTables/all_except_last_exon_NM.bed' -b $fileout/'4T1_67NR_peaksname.bed' -loj -s > $fileout/'e_le_4T1_67NR_NM.bed'
# 	$bedtools intersect -a $fileout/'genesTables/last_exon_NM.bed' -b $fileout/'4T1_U1_peaksname.bed' -loj -s > $fileout/'le_4T1_U1_NM.bed'
# 	$bedtools intersect -a $fileout/'genesTables/last_exon_NM.bed' -b $fileout/'4T1_67NR_peaksname.bed' -loj -s > $fileout/'le_4T1_67NR_NM.bed'
# 	#For NM_genes
# 	$bedtools intersect -a $fileout/'genesTables/all_except_last_exon_NR.bed' -b $fileout/'4T1_U1_peaksname.bed' -loj -s > $fileout/'e_le_4T1_U1_NR.bed'
# 	$bedtools intersect -a $fileout/'genesTables/all_except_last_exon_NR.bed' -b $fileout/'4T1_67NR_peaksname.bed' -loj -s > $fileout/'e_le_4T1_67NR_NR.bed'
# 	$bedtools intersect -a $fileout/'genesTables/last_exon_NR.bed' -b $fileout/'4T1_U1_peaksname.bed' -loj -s > $fileout/'le_4T1_U1_NR.bed'
# 	$bedtools intersect -a $fileout/'genesTables/last_exon_NR.bed' -b $fileout/'4T1_67NR_peaksname.bed' -loj -s > $fileout/'le_4T1_67NR_NR.bed'
	
# 	$R --vanilla -q --args $fileout/ e_le 4T1_U1 < $me/locPeaks_NM_9.r
# 	$R --vanilla -q --args $fileout/ le 4T1_67NR < $me/locPeaks_NM_9.r
# 	$R --vanilla -q --args $fileout/ e_le 4T1_67NR < $me/locPeaks_NM_9.r
# 	$R --vanilla -q --args $fileout/ le 4T1_U1 < $me/locPeaks_NM_9.r

# 	$R --vanilla -q --args $fileout/ e_le 4T1_U1 < $me/locPeaks_NR_10.r
# 	$R --vanilla -q --args $fileout/ le 4T1_67NR < $me/locPeaks_NR_10.r
# 	$R --vanilla -q --args $fileout/ e_le 4T1_67NR < $me/locPeaks_NR_10.r
# 	$R --vanilla -q --args $fileout/ le 4T1_U1 < $me/locPeaks_NR_10.r
	
	
	
# 	echo "Analyse DESeq"
# 	$R --vanilla -q --args $fileout/ e_le 4T1_U1 < $me/DESeq_Lastexon_11.r
# 	$R --vanilla -q --args $fileout/ e_le 4T1_U1 < $me/DESeq_LastPic_12.r
