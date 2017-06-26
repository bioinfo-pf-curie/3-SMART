#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import gzip
#import bz2

parser = argparse.ArgumentParser(description = 'SCRIPT TO SELECT THE SEQUENCES WITHOUT PROBLEM OF SEQUENCING')

parser.add_argument('-i','--input', help = 'input file (fastq casava)', required = True)
parser.add_argument("-o", "--output", help = "output file - fastq with N sequences", required = True) 
args = parser.parse_args()


input_file = args.input
output_file = args.output

input_file=gzip.open(input_file,'rb')
#input_file=bz2.BZ2File(input_file,'rb')

#input_file = open(input_file, 'rb')
output_file = open(output_file, 'w')


flag = 0
y_counter = 0

for line in input_file:	

	if (" " in line) and (":N:" in line):
		flag = 1

	if (" " in line) and (":Y:" in line):
		flag = 0
		y_counter += 1

	if flag == 1:
		output_file.write(line)

print "Discard ", y_counter, "reads"

input_file.close()
output_file.close()


## import pysam
## for read in samfile.fetch(XXX){
## read
## read.tag
## read.pos
## read.strand
## read.is_mapped ?
##}


