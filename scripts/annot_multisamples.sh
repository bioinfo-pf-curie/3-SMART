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
    echo "usage: $0 -c CONF_FILE -l INPUT_LIST  -s SAMPLES_TO_COMBINE -o OUTPUT_DIR"
}

if [ $# -eq 0 ]; then
    usage
    exit
fi

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) CONF_FILE=$2; shift;;
	(-l) INPUT_LIST=$2; shift;;
	(-s) COMB=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) die "$0: error - unrecognized option $1" 1>&2;;
	(*)  break;;
    esac
    shift
done


if [[ -z ${CONF_FILE} || -z ${INPUT_LIST} || -z ${OUTPUT} ]]; then
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
fi

## Read sample plan and put it in a hash table
declare -A SAMPLES
while read line  
do  
    id=`echo $line | ${AWK_PATH}/awk '{print $1}'`
    f=`echo $line | ${AWK_PATH}/awk '{print $2}'`
    SAMPLES[${id}]="${f}"
done < ${INPUT_LIST}


## For all combinaison
## Output for all samples comparison
odir=`echo $COMB | tr "," "_"`
OUTPUT_ALL=${OUTPUT}/${odir}
mkdir -p ${OUTPUT_ALL} || die "Cannot create output folder"
if [ -e ${OUTPUT_ALL}/merged_peaks.bed ]; then
    rm ${OUTPUT_ALL}/merged_peaks.bed
fi
    
for i in `echo $COMB | tr "," " "`
do  
    fn=${SAMPLES["$i"]}
    echo -e "--$fn--"
    dname=`basename $fn | sed -e 's/.bam//'`
    input=${OUTPUT}/$dname/*_filt.bed
    cat $input >> ${OUTPUT_ALL}/all_peaks.bed
done

${BEDTOOLS_PATH}/sortBed -i ${OUTPUT_ALL}/all_peaks.bed | ${BEDTOOLS_PATH}/bedtools merge -s -n -i - | ${AWK_PATH}/awk  'BEGIN{OFS="\t";c=1}{print $1,$2,$3,"mpeak_"c,$4,$5;c=c+1}' >> ${OUTPUT_ALL}/merged_peaks.bed
rm ${OUTPUT_ALL}/all_peaks.bed

## Annotate peaks - Annotation based on Gene level
echo "Annotate peaks ..."
cmd="${R_PATH}/R --vanilla CMD BATCH \"--args peakfile='${OUTPUT_ALL}/merged_peaks.bed' LEAnnotFile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' ELEAnnotFile='${ANNOT_DIR}/${ORG}/whole_gene_wo_last_exon.bed' TRSAnnotFile='${ANNOT_DIR}/${ORG}/transcriptsType.bed'  le_overlap='${MIN_LE_OV}' intron_overlap='${MIN_INTRON_OV}' minover='${MIN_ANNOT_OV}'\" ${SCRIPTS}/annot_peaks.R ${LOGS}/annot_peaks.Rout"
echo $cmd
eval $cmd

ANNO_PEAKS=`echo merged_peaks.bed | sed -e 's/.bed$/_NM_ELE_annot.bed/'`
if [[ ! -e ${OUTPUT_ALL}/${ANNO_PEAKS} ]]; then
    die "Error : Annotated peaks file not found"
fi

## Annotate peaks - Annotation based on Gene level
echo "Peaks description ..."
cmd="${R_PATH}/R --vanilla CMD BATCH \"--args org='${ORG}' peakfile='${OUTPUT_ALL}/merged_peaks_finallist.bed' polyAfile='${POLYA_MOTIF}'  wsizeup='${WINSIZE_UP}' wsizedown='${WINSIZE_DOWN}'\" ${SCRIPTS}/peaks_descriptor.R ${LOGS}/peaks_descriptor.Rout"
echo $cmd
eval $cmd





