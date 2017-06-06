#!/bin/bash

set -o pipefail  # trace ERR through pipes
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value

SOFT="3-SMART"
VERSION="2.0"

################
##  Function  ##

function evalecho {
    echo 
#    echo $1
    echo
    eval $1
}
################

##  Initialization  ##
dir=$(dirname $0)

NORMAL="\\033[0;39m" 
RED="\\033[1;31m"

die() {
    echo -e "$RED""$*""$NORMAL" 1>&2
    exit 1
}

##  Get args  ##
usage()
{
    echo "usage: $0 -c CONF_FILE -l INPUT_LIST -s STEP -o OUTPUT_DIR"
    echo "Use option -h for more information"
    exit
}

help()
{
    echo " ***$SOFT  $VERSION ***"
    echo
    echo "OPTIONS"
    echo "    -c CONFIG : configuration file for 3-SMART processing"
    echo "    -l LIST : input list file; 3 columns with the sample id, path and condition of samples to compare"
    echo "    -s STEP : run all or only a subset of the $SOFT workflow"
    echo "	  all : run all workflow"
    echo "	  merge_peaks : Concatenation of all peak files (3-SMART data)"
    echo " 	  annotate_peaks : Annotate peaks - Annotation based on Gene level"
    echo " 	  peaks_description : Create pictures with motifs and peaks"
    echo "    -o OUTPUT : output folder"
    echo "    [-h] : help"
    echo "    [-v] : version"
    exit

}

version()
{
    echo "$SOFT version $VERSION"
    exit
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
	(-s) STEP=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) help;;
	(-v) version;;
	(--) shift; break;;
	(-*) die "$0: error - unrecognized option $1" 1>&2;;
	(*)  break;;
    esac
    shift
done


if [[ -z ${CONF_FILE} || -z ${INPUT_LIST} || -z ${STEP} || -z ${OUTPUT} ]]; then
    usage
    exit
fi

##  Read configuration files  ##
source ${CONF_FILE}

##  Ouptut directory from per sample analysis  ##
if [ ! -d ${OUTPUT} ]; then
    echo "${OUTPUT} directory not found"
    exit
fi

##  For all combinaison  ##
##  Output for all samples comparison  ##
odir=`echo ${COMBINE_SAMPLE} | tr "," "_"`
OUTPUT_ALL=${OUTPUT}/${odir}
mkdir -p ${OUTPUT_ALL} || die "Cannot create output folder"
if [ -e ${OUTPUT_ALL}/merged_peaks.bed ]; then
    rm ${OUTPUT_ALL}/merged_peaks.bed
fi

##  Read sample plan and put it in a hash table  ##
declare -A SAMPLES
declare -A SAMPLES_NAME
while read line
do 
	id=`echo ${line} | ${AWK_PATH}/awk '{print $1}'`
	f=`echo ${line} | ${AWK_PATH}/awk '{print $2}'`
	n=`echo ${line} | ${AWK_PATH}/awk '{print $3}'`
	SAMPLES[${id}]="${f}"
	SAMPLES_NAME[${id}]="${n}"
done < ${INPUT_LIST}

##  STEP option  ##
for i in `echo ${STEP} | tr "," " " `
do
    NAME_STEP=${i}

##  logs folder  ##
    LOGS=${OUTPUT_ALL}/logs
    if [ ! -d ${LOGS} ]; then
	mkdir -p ${LOGS}
    fi


##  Print sample names to combine  ##
    echo "Combine samples: ${COMBINE_SAMPLE}"

########################################################
###  Concatenation of all peak files (3-SMART data)  ###
########################################################

## Concatenation of filtered peaks from different samples. After that, we decided to merge the overlapping peaks.

    if [[ ${NAME_STEP} == "merge_peaks" || ${NAME_STEP} == "all" ]]; then
	echo "Concatenate filtered peak files ..."
	if [ -e ${OUTPUT_ALL}/all_peaks.bed ]; then
	    rm ${OUTPUT_ALL}/all_peaks.bed
	fi
	for i in `echo ${COMBINE_SAMPLE} | tr "," " "`
	do
	    fn=${SAMPLES["${i}"]}
	    input="${fn}"
	    cat ${input} >> ${OUTPUT_ALL}/all_peaks.bed
	done

	name=`echo ${COMBINE_SAMPLE} | tr "," " "`
	if [ ! -f ${OUTPUT_ALL}/all_peaks.bed ]; then
	    die " error: the all_peaks.bed file doesn't exist" 1>&2
        fi
	echo "Merge Peaks ..."
        cmd="${BEDTOOLS_PATH}/sortBed -i ${OUTPUT_ALL}/all_peaks.bed | ${BEDTOOLS_PATH}/bedtools merge -s -n -i - | ${AWK_PATH}/awk 'BEGIN{OFS=\"\t\";c=1}{print \$1,\$2,\$3,\"mpeak_\"c,\$4,\$5;c=c+1}' > ${OUTPUT_ALL}/merged_peaks.bed"
        evalecho "$cmd"
    fi

#########################################################
###  Annotate peaks - Annotation based on Gene level  ###
#########################################################

## Annotation of the merged peaks

    if [[ ${NAME_STEP} == "annotate_peaks" || ${NAME_STEP} == "all" ]]; then
        if [ ! -f ${OUTPUT_ALL}/merged_peaks.bed ]; then
	    die " error: the merged_peaks.bed file doesn't exist" 1>&2
        fi
        echo "Annotate peaks ..."
        cmd="${R_PATH}/R CMD BATCH \"--args peakfile='${OUTPUT_ALL}/merged_peaks.bed' polyA_lib='${SCRIPTS}/polyA_lib.R' LEAnnotFile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' ELEAnnotFile='${ANNOT_DIR}/${ORG}/whole_gene_wo_last_exon.bed' TRSAnnotFile='${ANNOT_DIR}/${ORG}/transcriptsType.bed'  le_overlap='${MIN_LE_OV}' intron_overlap='${MIN_INTRON_OV}' minover='${MIN_ANNOT_OV}'\" ${SCRIPTS}/annot_peaks.R ${LOGS}/annot_peaks.Rout"
        evalecho "$cmd"
    fi

###########################
###  Peaks description  ###
###########################

## Creation of pictures and table about the merged peaks

    if [[ ${NAME_STEP} == "peaks_description" || ${NAME_STEP} == "all" ]]; then
        if [ ! -f ${OUTPUT_ALL}/merged_peaks_finallist.bed ]; then
	    die " error: the merged_peaks_finallist.bed file doesn't exist" 1>&2
        fi
	echo "Peaks description ..."
	cmd="${R_PATH}/R CMD BATCH \"--args org='${ORG}' peakfile='${OUTPUT_ALL}/merged_peaks_finallist.bed' polyA_lib='${SCRIPTS}/polyA_lib.R' polyAfile='${POLYA_MOTIF}' wsizeup='${WINSIZE_UP}' wsizedown='${WINSIZE_DOWN}'\" ${SCRIPTS}/peaks_descriptor.R ${LOGS}/peaks_descriptor.Rout"
	evalecho "$cmd"
    fi

done

