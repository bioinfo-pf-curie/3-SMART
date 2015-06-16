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

echo "Get Peaks ..."

## Quality Control
echo "Filter low quality reads (MAPQ < ${MIN_MAPQ}) ..."
BAM_PX=${OUTPUT}/`basename ${INPUT} | sed -e "s/.bam$/_q${MIN_MAPQ}/"`
BAMQ=${BAM_PX}.bam

${SAMTOOLS_PATH}/samtools view -b -q ${MIN_MAPQ} ${INPUT} | ${SAMTOOLS_PATH}/samtools sort -  ${BAM_PX} # @ option for parallelism ?
${SAMTOOLS_PATH}/samtools index ${BAMQ}

# Add some statistics about the quality filtering (half of the library removed, full of AAA..)

## Detect peaks
echo "Peak Detection ..."
BED_PEAKS=`echo ${BAMQ} | sed -e "s/.bam$/_peaks.bed/"`
${BEDTOOLS_PATH}/bedtools bamtobed -i  ${BAMQ} | ${BEDTOOLS_PATH}/bedtools merge -s -n -d ${MAX_DIST_MERGE} -i stdin | ${AWK_PATH}/awk -v thr=${MIN_NB_READS_PER_PEAK} 'BEGIN{OFS="\t";c=1}($4>=thr){print $1,$2,$3,"peak_"c,$4,$5;c=c+1}' > ${BED_PEAKS}

## Add statistics on peak size /reads distribution 
## TODO

## Filter peaks - stretchA 
echo "Filter peaks ..."
cmd="${R_PATH}/R --vanilla CMD BATCH \"--args peakfile='${BED_PEAKS}' keep_le_peaks='${KEEP_LE_PEAKS}'  LEAnnotFile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' wsizeup='${WINSIZE_UP}' wsizedown='${WINSIZE_DOWN}' nstretch='${NB_STRETCH_POLYA}' mism='${MISM}' nstretchcons='${NB_STRETCH_CONSECUTIVE}' org='${ORG}'\" ${SCRIPTS}/peaks_filter.R ${LOGS}/peaks_filter.Rout"
echo $cmd
eval $cmd

FILT_PEAKS=`echo ${BED_PEAKS} | sed -e 's/.bed$/_filt.bed/'`
if [[ ! -e ${FILT_PEAKS} ]]; then
    die "Error : Filtered peaks file not found"
fi

