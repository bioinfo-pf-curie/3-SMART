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
usage()
{
    echo "usage: $0 -p CONF_FILE -i INPUT -o OUTPUT_FILE"
}

if [ $# -eq 0 ]; then
    usage
    exit
fi

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) CONF_FILE=$2; shift;;
	(-i) INPUT=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) die "$0: error - unrecognized option $1" 1>&2;;
	(*)  break;;
    esac
    shift
done


if [[ -z ${CONF_FILE} || -z ${INPUT} || -z ${OUTPUT} ]]; then
    usage
    exit
fi

## Read configuration files
source ${CONF_FILE}

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

## Ouptut directory from per sample analysis
if [ ! -d ${OUTPUT} ]; then
    echo "${OUTPUT} directory not found"
    exit
else
    out=`basename ${INPUT} | sed -e 's/.bam//'`
    FILT_PEAKS=${OUTPUT}/${out}"_q20_peaks_filt.bed"
fi


## Annotate peaks - Annotation based on Gene level

echo "Annotate peaks ..."

cmd="${R_PATH}/R --vanilla CMD BATCH \"--args peakfile='${FILT_PEAKS}' LEAnnotFile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' ELEAnnotFile='${ANNOT_DIR}/${ORG}/whole_gene_wo_last_exon.bed' TRSAnnotFile='${ANNOT_DIR}/${ORG}/transcriptsType.bed'  le_overlap='${MIN_LE_OV}' intron_overlap='${MIN_INTRON_OV}'  minover='${MIN_ANNOT_OV}'\" ${SCRIPTS}/annot_peaks.R ${LOGS}/annot_peaks.Rout"
echo $cmd
eval $cmd

FINAL_PEAKS=`echo ${FILT_PEAKS} | sed -e 's/.bed$/_finallist.bed/'`
echo $FINAL_PEAKS

if [[ ! -e ${FINAL_PEAKS} ]]; then
    die "Error : final peaks file not found"
fi


## Describe peaks

echo "Peaks description ..."

cmd="${R_PATH}/R --vanilla CMD BATCH \"--args org='${ORG}' peakfile='${FINAL_PEAKS}' polyAfile='${POLYA_MOTIF}'  wsizeup='${WINSIZE_UP}' wsizedown='${WINSIZE_DOWN}'\" ${SCRIPTS}/peaks_descriptor.R ${LOGS}/peaks_descriptor.Rout"
echo $cmd
eval $cmd
