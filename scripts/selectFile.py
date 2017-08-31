#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import gzip

parser = argparse.ArgumentParser(description = 'SCRIPT TO SELECT THE SEQUENCES WITHOUT PROBLEM OF SEQUENCING')

parser.add_argument('-i','--input', help = 'input file (fastq casava)', required = True)
parser.add_argument("-o", "--output", help = "output file - fastq with N sequences", required = True) 
args = parser.parse_args()


input_file = args.input
output_file = args.output


input_file=gzip.open(input_file,'rb')
output_file = open(output_file, 'w')


flag = 0
y_counter = 0
n_counter = 0

for line in input_file:	

	if (" " in line) and (":N:" in line):
		flag = 1
		n_counter += 1

	if (" " in line) and (":Y:" in line):
		flag = 0
		y_counter += 1

	if flag == 1:
		output_file.write(line)

print "total reads", "\t", "reads discarded", "\n", n_counter + y_counter, "\t", y_counter 

input_file.close()
output_file.close()


