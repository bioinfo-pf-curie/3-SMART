.. Nicolas Servant
.. PolyA-seq pipeline
.. v1.0.0
.. 15-04-01


SMART Quick Start Guide
*************************

See NEWS for information about changes in this and previous versions

See LOGBOOK for details about the HiC-Pro developmens

What is the SMART pipeline ?
============================

SMART (3'-Seq Mapping Annotation and Regulation Tool) was set up to process 3' sequencing data from aligned sequencing reads.
It includes the peak detection, filtering, and annotation.
If you use it, please cite :
*Boldina et al. 2015*


How to install it ?
===================

The following dependancies are required :
* 
* R with the *RColorBrewer*, *ggplot2*, *rtracklayer*, *DESeq2*, *magrittr*, *dplyr*, *qplots* and *GenomicRanges* packages
* Samtools (>0.1.18)
* BEDTools (>2.17.0)

To install the SMART pipeline, simply extract the archive and set up the configuration file with the paths to dependencies.

.. code-block:: guess

  tar -zxvf smart_1.0.0.tar.gz
  cd smart_1.0.0


Annotation Files
================

The pipeline is using a couple of annotation files with gene annotation and last exons information. These files are based on UCSC Refseq gene.
In order to generate all required annotation files, please set the ANNOT_DIR, ORG and UCSC_EXPORT in the configuration file

.. code-block:: guess

  BUILD_ANNOT=1
  ORG=mm9
  UCSC_EXPORT=refseq_export_mm9.csv

SMART will start by creating the annotation files in the forder ANNOT_DIR/ORG/
The annotation only have to be generated once. Mouse annotation are provided with the pipeline.


How to use it ?
===============

SMART can be used both for a single sample or for a list of samples.
In the case of sample list, peaks detected in all samples are merged before annotation. 
In oder to use the pipeline, please set up the configuration according to your analysis, and run the following command :

.. code-block:: guess

   /bin/smart -c CONFIG [-i INPUT_BAM] [-l INPUT_LIST] -o OUTPUT_DIR

