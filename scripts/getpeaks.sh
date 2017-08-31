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
usage() {
    echo "usage: $0 -c CONF_FILE -i INPUT_FILE -s STEP -n SAMPLE -o OUTPUT_DIR"
    echo "Use option -h for more information"
    exit
}

help()
{
    echo " ***$SOFT  $VERSION ***"
    echo
    echo "OPTIONS"
    echo "    -c CONFIG : configuration file for 3-SMART processing"
    echo "    -i INPUT : input file; fastq.gz format"
    echo "    -s STEP : run all or only a subset of the $SOFT workflow"
    echo "	  all : run all workflow"
    echo "	  trimming : Remove reads which have a problem with sequencer and remove adapters"
    echo " 	  remove_stretchA : Remove A stretch in 3' side and select a read size"
    echo " 	  fastqc : Quality control for fastq file"
    echo " 	  mapping : Map reads on the reference genome with Bowtie2"
    echo " 	  duplicate_reads : Remove or only mark the duplicated reads"
    echo " 	  peak_calling : Detect peaks"
    echo " 	  filter_peaks : Select peaks without stretch of genomic A in reads"
    echo "    -n SAMPLE : Sample id (number or letter; no symbols)"
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
	(-c) CONF=$2; shift;;
	(-i) INPUT=$2; shift;;
	(-s) STEP=$2; shift;;
	(-n) SAMPLE=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) help;;
	(-v) version;;
	(--) shift; break;;
	(-*) die "$0: error - unrecognized option $1" 1>&2;;
	(*)  break;;
    esac
    shift
done

if [[ -z ${CONF} || -z ${INPUT} || -z ${STEP} || -z ${SAMPLE} || -z ${OUTPUT} ]]; then
    usage
    exit
fi

##  Read configuration files  ##
source ${CONF}

##  Output file  ##
if [ ! -d ${OUTPUT} ]; then
    mkdir -p ${OUTPUT}
fi

if [ ! -d ${OUTPUT}/${SAMPLE} ]; then
    mkdir -p ${OUTPUT}/${SAMPLE}
fi

##  STEP option  ##
for i in `echo ${STEP} | tr "," " " `
do
    NAME_STEP=${i}

##  logs folder  ##
    LOGS=${OUTPUT}/${SAMPLE}/logs
    if [ ! -d ${LOGS} ]; then
	mkdir -p ${LOGS}
    fi


################################
###  Clean fastq + Trimming  ###
################################

## clean_fastq step removes reads which have had a problem with the sequencer
## trimming step removes adapters in 5' and 3' side of reads

    OUTPUT_CLEAN=${OUTPUT}/${SAMPLE}/`basename ${INPUT} | sed -e "s/.fastq.gz$/_clean.fastq/"`
    OUTPUT_TRIMMING=${OUTPUT}/${SAMPLE}/`basename ${OUTPUT_CLEAN} | sed -e 's/_clean.fastq$/_trimming_filter.fastq/'`
    if [[ ${NAME_STEP} == "trimming" || ${NAME_STEP} == "all" ]]; then
	echo "Clean fastq ..."
	cmd="${PYTHON_PATH}/python ${SCRIPTS}/selectFile.py -i ${INPUT} -o ${OUTPUT_CLEAN}"
	evalecho "$cmd"

        if [ ! -f ${OUTPUT_CLEAN} ]; then
	    die " error: the *_clean.fastq file doesn't exist" 1>&2
	fi
	echo "Trimmed data ..."
	cmd="${CUTADAPT_PATH}/cutadapt -g ${RD_1_SP} -a ${INDEX_SP} -n 2 -o ${OUTPUT_TRIMMING} ${OUTPUT_CLEAN} 2>&1 > ${LOGS}/cutadapt.log"
	evalecho "$cmd"
    fi

##############################################################################
###  Remove stretches of A in 3' + Select minimum of nucleotides per read  ###
##############################################################################

## remove_stretchA step permits to remove the A stretch (with a minimum of 3 A in 3' side) and select reads with a minimum of nucleotides (${MIN_LENGTH_READS})

    OUTPUT_STRETCH=${OUTPUT}/${SAMPLE}/`basename ${OUTPUT_TRIMMING} | sed -e 's/_trimming_filter.fastq$/_remove_stretchA.fastq/'`
    if [[ ${NAME_STEP} == "remove_stretchA" || ${NAME_STEP} == "all" ]]; then
## File verification
	if [ ! -f ${OUTPUT_TRIMMING} ]; then
	    die " error: the *_trimming_filter.fastq file doesn't exist" 1>&2
	fi
	echo "Remove A stretch in 3' side..."
	cmd="${PYTHON_PATH}/python ${SCRIPTS}/removeStretchPolyA.py -input ${OUTPUT_TRIMMING} -output ${OUTPUT_STRETCH} -length ${MIN_LENGTH_READS}"
	evalecho "$cmd"
    fi

################
###  FASTQC  ###
################
    if [[ ${NAME_STEP} == "fastqc" || ${NAME_STEP} == "all" ]]; then
	if [ ! -f ${OUTPUT_STRETCH} ]; then
	    die " error: the *_remove_stretchA.fastq file doesn't exist" 1>&2
	fi
	echo "Fastqc ..."
        cmd="${FASTQC_PATH}/fastqc ${OUTPUT_STRETCH}"
        evalecho "$cmd"
    fi

#################
###  Mapping  ###
#################

## In this step, we use Bowtie2 as mapper. We decided to apply the MAPQ and after that, sort and index this file.

    OUTPUT_SORT=${OUTPUT}/${SAMPLE}/`basename ${OUTPUT_STRETCH} | sed -e 's/_remove_stretchA.fastq$/_MAPQ_sort.bam/'`

    if [[ ${NAME_STEP} == "mapping" || ${NAME_STEP} == "all" ]]; then
	if [ ! -f ${OUTPUT_STRETCH} ]; then
	    die " error: the *_remove_stretchA.fastq file doesn't exist" 1>&2
	fi
	## mapping_Tophat folder
	MAPPING=${OUTPUT}/${SAMPLE}/mapping_Bowtie2
	if [ ! -d ${MAPPING} ]; then
	    mkdir -p ${MAPPING}
	fi
	echo "Mapping with Bowtie2 ..."
	cmd="${BOWTIE2_PATH}/bowtie2 -x ${BOWTIE2_INDEX} -q ${OUTPUT_STRETCH} -N ${MAX_MISMATCH_SEED} -L ${SEED_LENGTH} -k ${RANDOM_HIT} -S ${MAPPING}/accepted_hits.sam"
	evalecho "$cmd"

	if [ ! -f ${MAPPING}/accepted_hits.sam ]; then
	    die " error: the accepted_hits.sam file in mapping_Bowtie2 folder doesn't exist" 1>&2
	fi
	echo "Mapping Quality >= ${MIN_MAPQ} ..."
	echo "Sorting ..."
	cmd="${SAMTOOLS_PATH}/samtools view -b -q ${MIN_MAPQ} ${MAPPING}/accepted_hits.sam | ${SAMTOOLS_PATH}/samtools sort -O bam -T prefix.bam - -o ${OUTPUT_SORT}"
	evalecho "$cmd"
	cmd="${SAMTOOLS_PATH}/samtools index ${OUTPUT_SORT}"
	evalecho "$cmd"
    fi

################################
###  Remove Duplicate reads  ###
################################

## duplicate_reads step permits to remove OR mark duplicated reads. (Duplicated read is a reads with the same sequence and the same length)

    OUTPUT_RM_DUPLICATED_READS=${OUTPUT}/${SAMPLE}/`basename ${OUTPUT_SORT} | sed -e 's/_MAPQ_sort.bam$/_rm_duplicated_reads.bam/'`

    if [[ ${NAME_STEP} == "duplicate_reads" || ${NAME_STEP} == "all" ]]; then
	if [ ! -f ${OUTPUT_SORT} ]; then
	    die " error: the *_MAPQ_sort.bam file doesn't exist" 1>&2
	fi
	if [ "${REMOVE_DUPLICATES}" = 1 ]; then
	    echo "Remove Duplicated reads ..."
	    cmd="${JAVA_PATH}/java -Xmx30G -jar ${PICARD_TOOLS}/MarkDuplicates.jar I=${OUTPUT_SORT} O=${OUTPUT_RM_DUPLICATED_READS} METRICS_FILE=${OUTPUT}/${SAMPLE}/metrix.txt REMOVE_DUPLICATES=TRUE"
	    evalecho "$cmd"
	else
	    echo "Don't remove Duplicated reads"
	    cmd="${JAVA_PATH}/java -jar ${PICARD_TOOLS}/MarkDuplicates.jar I=${OUTPUT_SORT} O=${OUTPUT_RM_DUPLICATED_READS} METRICS_FILE=${OUTPUT}/${SAMPLE}/metrix.txt REMOVE_DUPLICATES=FALSE"
	    evalecho "$cmd"
	fi
    fi

######################
###  Detect peaks  ###
######################

## peak_calling step permits to detect peaks

    OUTPUT_PEAKS=${OUTPUT}/${SAMPLE}/`basename ${OUTPUT_SORT} | sed -e 's/_MAPQ_sort.bam$/_peaks.bed/'`

    if [[ ${NAME_STEP} == "peak_calling" || ${NAME_STEP} == "all" ]]; then
	if [ ! -f ${OUTPUT_SORT} ]; then
	    die " error: the *_MAPQ_sort.bam file doesn't exist" 1>&2
	fi
	echo "Get Peaks ..."
        echo "Peak Detection ..."
        cmd="${BEDTOOLS_PATH}/bedtools bamtobed -i  ${OUTPUT_SORT} | ${BEDTOOLS_PATH}/bedtools merge -s -d ${MAX_DIST_MERGE} -c 4 -o count -i stdin | ${AWK_PATH}/awk -v thr=\${MIN_NB_READS_PER_PEAK} 'BEGIN{OFS=\"\t\";c=1}(\$5>=thr){print \$1,\$2,\$3,\"peak_\"c,\$5,\$4;c=c+1}' > ${OUTPUT_PEAKS}"
	evalecho "$cmd"
    fi

################################
###  Build Annotation files  ###
################################
    if [[ ${BUILD_ANNOT} == 1 ]]; then
        echo "Build annotation files ..."             
        ## Create gene file
	cmd="${R_PATH}/R CMD BATCH \"--args infile='${ANNOT_DIR}/${ORG}/${UCSC_EXPORT}' out_dir='${ANNOT_DIR}/${ORG}'\" ${SCRIPTS}/make_annot_gene.R ${LOGS}/make_annot_gene.Rout"
        evalecho "$cmd"
        ## Create exon file
        cmd="${R_PATH}/R --vanilla CMD BATCH \"--args txfile='${ANNOT_DIR}/${ORG}/${UCSC_EXPORT}' lefile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' out='${ANNOT_DIR}/${ORG}'\" ${SCRIPTS}/make_annot_exon.R ${LOGS}/make_annot_exon.Rout"
       evalecho "$cmd"
    fi 

#################################
###  Filter peaks - stretchA  ###
#################################

## In this step, we remove peaks which present a stretch of genomic A (*_peaks_genomstretch.bed) and keep filtered peaks (*_peaks_filt.bed)

    OUTPUT_FILT_PEAKS=${OUTPUT}/${SAMPLE}/`basename ${OUTPUT_PEAKS} | sed -e 's/.bed$/_filt.bed/'`

    if [[ ${NAME_STEP} == "filter_peaks" || ${NAME_STEP} == "all" ]]; then
        echo "Filter peaks ..."
        cmd="${R_PATH}/R CMD BATCH \"--args peakfile='${OUTPUT_PEAKS}' keep_le_peaks='${KEEP_LE_PEAKS}' polyA_lib='${SCRIPTS}/polyA_lib.R' LEAnnotFile='${ANNOT_DIR}/${ORG}/last_exon_gene.bed' wsizeup='${WINSIZE_UP}' wsizedown='${WINSIZE_DOWN}' nstretch='${NB_STRETCH_POLYA}' mism='${MISM}' nstretchcons='${NB_STRETCH_CONSECUTIVE}' lexogen='${LEXOGEN}' org='${ORG}'\" ${SCRIPTS}/peaks_filter.R ${LOGS}/peaks_filter.Rout"
        evalecho "$cmd"

        if [[ ! -e ${OUTPUT_FILT_PEAKS} ]]; then
            die "error : the *_peaks_filt.bed file doesn't exist"
        fi
    fi


done


