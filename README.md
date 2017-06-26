3-SMART Quick Start Guide
======================

See NEWS for information about changes in this and previous versions

See LOGBOOK for details about the 3-SMART developments

What is the 3-SMART pipeline ?
-------------------------------

3-SMART (3'-Seq Mapping Annotation and Regulation Tool) was setup to process 3' sequencing data from sequencing reads.
It includes the trimming, mapping, peak detection, filtering, annotation and differential analysis.

Contact
-------

For any questions about the pipeline, please contact <mandy.cadix@curie.fr>

How to install it ?
-------------------

The following dependencies are required :

* [R] (https://www.r-project.org/) (version 3.2.3) with the *RColorBrewer*, *ggplot2*, *rtracklayer*, *DESeq2*, *magrittr*, *dplyr*, *gplots*, *plyr* and *GenomicRanges* packages

* [Samtools] (http://samtools.sourceforge.net) (version 1.1)

* [BEDTools] (http://bedtools.readthedocs.io/en/latest/content/installation.html) (version 2.17.0)

* [Python] (https://www.python.org/downloads/release/python-279) (version 2.7.9) with the *os*, *argparse*, *re* and *gzip* packages

* [Cutadapt] (http://cutadapt.readthedocs.io/en/stable/installation.html) (version 1.8.2)

* [Bowtie2] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version 2.2.5)

* [Java/jdk] (http://www.oracle.com/technetwork/java/javase/downloads/index-jsp-138363.html) (version 1.7.0\_45)

* [Picard\_tools] (https://broadinstitute.github.io/picard/)

* [Fastqc] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.10.1)


To install the 3-SMART pipeline, simply extract the archive and set up the configuration file with the paths to dependencies.

    unzip 3-SMART-2.0.zip
    cd 3-SMART-2.0


Annotation Files
----------------

The pipeline is using a couple of annotation files with gene annotation and last exons information. These files are based on UCSC Refseq gene.
In order to generate all required annotation files, please set the ANNOT\_DIR, ORG and UCSC\_EXPORT in the configuration file.


    BUILD_ANNOT=1
    ORG=hg19
    UCSC_EXPORT=refseq_export_hg19.csv

3-SMART will start by creating the annotation files in the forder ANNOT\_DIR/ORG/
The annotation only have to be generated once. Human annotation are provided with the pipeline.

Trimming
--------

To treat reads with adapter, it is necessary to know adapters 5' (P5) and 3' (P7):


    ADAPTERP5="GTTCAGAGTTCTACAGTCCGACGATC"
    ADAPTERP7="TGGAATTCTCGGGTGCCAAGG"


How to use it ?
---------------

**PART 1 : reads processing** 3-SMART can be used for a single sample. In order to use the pipeline, please set up the configuration according to your analysis, and run the following command to do the preprocessing data:


    /script/getpeaks.sh -c CONFIG -i INPUT_FILE -s STEP -n SAMPLE -o OUTPUT_DIR [-h] [-v]

    -c: The configuration file

    -i: The fastq format file of the CLIP-seq sample. CAUTION: only the files in fastq.gz format are accepted.

    -s: The steps of this pipeline that is to be launched. There are 6 differents steps: 
	trimming: Remove reads which have had problems with the sequencer and remove adapters
	remove_stretchA: Remove A stretch in 3' side and select a read size
	fastqc: Quality control for fastq file
	mapping: Map reads with bowtie2 and use the mapping quality to select reads. On this step we create a sort and index file.
	duplicate_reads: Remove (or only mark) the duplicated reads
	peak_calling: Detect peaks
	filter_peaks: Select peaks without stretch of genomic A in reads
	all: Launch all previous step

    -n: The sample identifiers (number or letter; no symbols)

    -o: The output directory

    [-h]: Help

    [-v]: Version



**PART 2 : annotation multisamples** Annotate peaks present in multiple samples. Run the following command to do the annotation of multiple samples:


    /script/annot_multisamples.sh -c CONFIG -l INPUT_LIST -s STEP -o OUTPUT_DIR [-h] [-v]

    -c: The configuration file

    -l: The list of BED files (SAMPLE_NUMBER_peaks_filt.bed), obtained from the PART 1 (See input_list_annotMultisamples.txt)

    -s: The steps of this pipeline that is to be launched. There are 3 differents steps: 
	merge_peaks: Concatenate, sort and merge all peaks present in different samples
	annotate_peaks: Annotate the merged peaks
	peaks_description: Create different pictures and table
	all: Launch all previous step

    -o: The output directory

    [-h]: Help message

    [-v]: Version 

This is a input list of BED files obtained thanks to the PART 1, for the PART 2, *input_list_annotMultisamples.txt*


    Sample2        path/Sample2/Sample2_test_peaks_filt.bed        Test
    Sample4        path/Sample4/Sample4_test_peaks_filt.bed        Test
    Sample1        path/Sample1/Sample1_ctrl_peaks_filt.bed        Ctrl
    Sample3        path/Sample3/Sample3_ctrl_peaks_filt.bed        Ctrl



**PART 3 : differential analysis** If you want to compare treated VS untreated please enter treated samples before untreated samples. Run the following command to do the differential analysis:


    /script/compare_samples.sh -c CONFIG -l INPUT_LIST -s STEP -o OUTPUT_DIR

    -c: The configuration file

    -l: The list of BAM files, obtained from the PART 1 (See input_list_compare.txt)

    -s: The steps of this pipeline that is to be launched. There are 2 differents steps:
	quantification: Quantify merged peaks and remove peaks which are on different genes
	differential_analysis: Make the differential analysis
	all: Launch all previous step

    -o: The output directory

    [-h]: Help message

    [-v]: Version

This is a input list of BAM files obtained thanks to the PART 1, for the PART 3, *input_list_compare.txt*


    Sample2        path/Sample2/Sample2_test_MAPQ_sort.bam        Test
    Sample4        path/Sample4/Sample4_test_MAPQ_sort.bam        Test
    Sample1        path/Sample1/Sample1_ctrl_MAPQ_sort.bam        Ctrl
    Sample3        path/Sample3/Sample3_ctrl_MAPQ_sort.bam        Ctrl


