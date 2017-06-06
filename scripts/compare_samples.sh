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

die()
{
    echo -e "$RED""$*""$NORMAL" 1>&2
    exit 1
}

##  Get args  ##
usage()
{
    echo "usage: $SOFT $0 -c CONF_FILE -l INPUT_LIST -s STEP -o OUTPUT_DIR"
    echo "Use option -h for more information"
    exit
}

help()
{
    echo " *** $SOFT  $VERSION ***"
    echo
    echo "OPTIONS"
    echo "    -c CONFIG : configuration file for 3-SMART processing"
    echo "    -l LIST : input list file; 3 columns with the sample id, path of bam file and condition of samples to compare"
    echo "    -s STEP : run all or only a subset of the $SOFT workflow"
    echo "	  all : run all workflow"
    echo "	  quantification : create a table of counts of samples to compare and filter them"
    echo " 	  differential_analysis : do a differential analysis with samples to compare"
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
	(-s) STEP=$2; shift;;
	(-l) INPUT_LIST=$2; shift;;
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

##  Ouptut file  ##
if [ ! -d ${OUTPUT} ]; then
    mkdir -p ${OUTPUT}
fi

##  Read sample plan and put it in a hash table  ##
declare -A SAMPLES
declare -A SAMPLES_NAME
while read line  
do  
    id=`echo $line | ${AWK_PATH}/awk '{print $1}'`
    f=`echo $line | ${AWK_PATH}/awk '{print $2}'`
    n=`echo $line | ${AWK_PATH}/awk '{print $3}'`
    SAMPLES[${id}]="${f}"
    SAMPLES_NAME[${id}]="${n}"
done < ${INPUT_LIST}

##  Output for all samples comparison  ##
odir=`echo ${COMBINE_SAMPLE} | tr "," "_"`
OUTPUT_ALL=${OUTPUT}/${odir}

##  STEP option  ##
for i in `echo ${STEP} | tr "," " " `
do
    NAME_STEP=${i}

##  logs folder  ##
    LOGS=${OUTPUT_ALL}/logs
    if [ ! -d ${LOGS} ]; then
	mkdir -p ${LOGS}
    fi

########################################
###  Quantification of merged peaks  ###
########################################

## Creation of a table of counts about merged peaks (to have the number of reads for each merged peaks for each sample)
## After that, we remove peaks which are on different genes. Here we differencie NM and NR genes. If one peak is present on NM gene and NR gene, we remove this peak.

    if [[ ${NAME_STEP} == "quantification" || ${NAME_STEP} == "all" ]]; then
        if [ ! -f ${OUTPUT_ALL}/merged_peaks_finallist_annot.tsv ]; then
	    die " error: the merged_peaks_finallist_annot.tsv file doesn't exist" 1>&2
        fi
	echo "Peak quantification ..."
	file=""
	for i in `echo ${COMBINE_SAMPLE} | tr "," " "`
    	do  
            fn=${SAMPLES["$i"]}
    	    dname=`basename $fn | sed -e 's/.bam//'`
	    input=${OUTPUT}/${dname}/${dname}_MAPQ_sort.bam
	    file="$file ${fn}"
        done

        cmd="${AWK_PATH}/awk '(NR>1){OFS=\"\t\";print \$1,\$2,\$3,\$4,\$4,\$5,\$6,\$7,\$8}' ${OUTPUT_ALL}/merged_peaks_finallist_annot.tsv  | grep "mpeak" | ${BEDTOOLS_PATH}/bedtools multicov -s -bams $file -bed - > ${OUTPUT_ALL}/merged_peaks_finallist_annot_qt.bed"
        evalecho "$cmd"


	if [ ! -f ${OUTPUT_ALL}/merged_peaks_finallist_annot_qt.bed ]; then
	    die " error: the merged_peaks_finallist_annot_qt.bed file doesn't exist" 1>&2
        fi
	echo "Remove peaks which are on different genes ..."
	name=`echo ${COMBINE_SAMPLE} | tr "," "\t"`
	cmd="${AWK_PATH}/awk -v name=\"\$name\" 'BEGIN{OFS=\"\t\";print \"chr\",\"start\",\"end\",\"length\",\"length2\",\"strand\",\"peak\",\"merge\",\"status\",name,\"gene\";}{name=\"\";p=\"T\";split(\$7,a,\"|\");for(x in a){if(a[x] ~ /^chr/){split(a[x],b,\"_\");if(name==\"\"){name=b[2]\"_\"b[3];}else{if(p==\"T\" && name != b[2]\"_\"b[3]){p=\"F\";}}}};  if(p==\"T\"){OFS=\"\t\";print \$0,\$1\"_\"name;} }' ${OUTPUT_ALL}/merged_peaks_finallist_annot_qt.bed > ${OUTPUT_ALL}/merged_peaks_finallist_IanDif.bed"
        evalecho "$cmd"
    fi

###############################
###  Differential analysis  ###
###############################
    if [[ ${NAME_STEP} == "differential_analysis" || ${NAME_STEP} == "all" ]]; then
        if [ ! -f ${OUTPUT_ALL}/merged_peaks_finallist_IanDif.bed ]; then
	    die " error: the merged_peaks_finallist_IanDif.bed file doesn't exist" 1>&2
        fi
	echo "Differential analysis ..."
	cmd="${R_PATH}/R CMD BATCH \"--args peakfile='${OUTPUT_ALL}/merged_peaks_finallist_IanDif.bed' polyA_lib='${SCRIPTS}/polyA_lib.R' input_list='${INPUT_LIST}' group='${COMPARE_SAMPLE}' min_count_cond='${MIN_COUNT_PER_COND}' ratioValue='${PERCENT_VALUE_RATIO}'\" ${SCRIPTS}/compare.R ${LOGS}/compare.Rout"
	evalecho "$cmd"
    fi

done


