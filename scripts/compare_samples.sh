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
    echo "usage: $0 -c CONF_FILE -l INPUT_LIST -s SAMPLE_TO_COMBINE -g SAMPLE_GROUP -o OUTPUT_DIR"
}

if [ $# -eq 0 ]; then
    usage
    exit
fi

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) CONF_FILE=$2; shift;;
	(-s) COMB=$2; shift;;
	(-l) INPUT_LIST=$2; shift;;
	(-g) GROUP=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) die "$0: error - unrecognized option $1" 1>&2;;
	(*)  break;;
    esac
    shift
done


if [[ -z ${CONF_FILE} || -z ${INPUT_LIST} || -z ${GROUP} || -z ${COMB} || -z ${OUTPUT} ]]; then
    usage
    exit
fi

## Read configuration files
source ${CONF_FILE}

## Ouptut file
if [ ! -d ${OUTPUT} ]; then
    mkdir -p ${OUTPUT}
fi

## Read sample plan and put it in a hash table

## Read sample plan and put it in a hash table
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

## Output for all samples comparison
omdir=`echo $COMB | tr "," "_"`
OUTPUT_MERGE=${OUTPUT}/${omdir}

## Create condition list
#COND_LIST=$(cut -f2 ${INPUT_LIST} | awk '{split($0,a," "); FS=",";for(x in a){OFS=";";print "\42"a[x]"\42"} }')
#COND_LIST=$(cut -f2 ${INPUT_LIST} )
#COND_LIST=$(echo $COND_LIST | sed "s/ /,/g")

### Annotate peaks
#echo "Annotate peaks ..."
#${R_PATH}/R --vanilla CMD BATCH "--args peakfile='${OUTPUT_MERGE}/merged_peaks_finallist.bed' LEAnnotFile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' ELEAnnotFile='${ANNOT_DIR}/${ORG}/whole_gene_wo_last_exon.bed' TRSAnnotFile='${ANNOT_DIR}/${ORG}/transcriptsType.bed'  l#e_overlap='${MIN_LE_OV}' intron_overlap='${MIN_INTRON_OV}' minover='${MIN_ANNOT_OV}'" ${SCRIPTS}/annot_peaks_mod.R ${LOGS}/annot_peaks.Rout

### quantification of merged peaks
echo "Peak quantification ..."
file=""
##name=""
for i in `echo $COMB | tr "," " "`
do  
    fn=${SAMPLES["$i"]}
    dname=`basename $fn | sed -e 's/.bam//'`
    input=${OUTPUT}/${dname}/${dname}_q20.bam
    file="$file ${OUTPUT}/${dname}/${dname}_q20.bam"
    ##name="$name\t"${SAMPLES_NAME["$i"]}
done

name=`echo $COMB | tr "," "\t"`

${AWK_PATH}/awk '(NR>1){OFS="\t";print $1,$2,$3,$4,$4,$5,$6,$7,$8}' ${OUTPUT_MERGE}/merged_peaks_finallist_annot.tsv  | grep "mpeak" | ${BEDTOOLS_PATH}/bedtools multicov -s -bams $file -bed - > ${OUTPUT_MERGE}/merged_peaks_finallist_annot_qt.bed

## Remove peaks which are on different genes
${AWK_PATH}/awk -v name="$name" 'BEGIN{OFS="\t";print "chr","start","end","length","length2","strand","peak","merge","status",name,"gene";}{name="";p="T";split($7,a,"|");for(x in a){if(a[x] ~ /^chr/){split(a[x],b,"_");if(name==""){name=b[2];}else{if(p=="T" && name != b[2]){p="F";}}}};  if(p=="T"){OFS="\t";print $0,$1"_"name;} }' ${OUTPUT_MERGE}/merged_peaks_finallist_annot_qt.bed > ${OUTPUT_MERGE}/merged_peaks_finallist_IanDif.bed

### differential analysis
echo "Differential analysis ..."
cmd="/bioinfo/local/build/R/R-3.1.0/bin/R --vanilla CMD BATCH \"--args peakfile='${OUTPUT_MERGE}/merged_peaks_finallist_IanDif.bed' input_list='${INPUT_LIST}' group='${GROUP}' min_count_cond='${MIN_COUNT_PER_COND}'\" ${SCRIPTS}/compare.R ${LOGS}/compare.Rout"
echo $cmd
eval $cmd

##/bioinfo/local/build/R/R-3.1.0/bin/R --vanilla CMD BATCH "--args peakfile='${OUTPUT_MERGE}/merged_peaks_finallist_IanDif.bed' input_list= '${INPUT_LIST}'condfile='${GROUP}'" ${SCRIPTS}/compare.R ${LOGS}/compare.Rout
