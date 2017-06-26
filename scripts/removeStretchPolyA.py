#!/usr/bin/python
# -*- coding: utf-8 -*-

###########################################################
#Mandy Cadix & Jocelyn Brayet
#
#usage: removeStretchPolyA.py [-h] -input <INPUT_FILE> -output <OUTPUT_FILE> -length <MIN_LENGTH_READS>
#
#Remove stretches of A in 3' side. 
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -input <INPUT_FILE>, --inputFile <INPUT_FILE>
#                        Input file (fastq format).
#  -output <OUTPUT_FILE>, --outputFile <OUTPUT_FILE>
#                        Output file (path/to/file/nameFile.fastq).
#  -length MIN_LENGTH_READS, --minLengthReads MIN_LENGTH_READS
#			 Select minimum number of nucleotides to keep reads
#
###########################################################

###########################################################
## Import
import argparse
import os
import re

## Tested with python 2.7.9
removeStretchVersion = '1.2 - 03/04/2017'

if __name__ == '__main__':

	########### arguments ####################
	parser = argparse.ArgumentParser(description='remove stretch of A', epilog='Version '+removeStretchVersion)
	
	parser.add_argument('-input', '--inputFile', type=str, action="store", default="", help='Input file (fastq format).', required=True)
	parser.add_argument('-output', '--outputFile', type=str, action="store", default="", help='Output file (path/to/file/nameFile.fastq).', required=True)
	parser.add_argument('-length', '--minLengthReads', type=str, action="store", default="", help='minimum length of reads', required=True)

	dargs = vars(parser.parse_args())

	inputFile = open(dargs["inputFile"],"r")
	outputFile = open(dargs["outputFile"],"w")
	minLengthReads = int(dargs["minLengthReads"])


	nbr_line = 0
	after = 0
	polyA_counter = 0
	writeRead = False
	minReads = 0
	totalReads = 0

	for line in inputFile:

		nbr_line=nbr_line+1

		if (nbr_line==1):
			line1 = line
			totalReads +=1

		if (nbr_line==3):
			line3 = line

		if (nbr_line==2 and re.match(r'(.*)A{3,}$', line)):
			after_string = re.findall(r'A{3,}$', line)
			after = len(after_string[0])+1
			line = re.sub(r'A{3,}$','',line)
			line2 = line
			polyA_counter += 1

		else:
			if (nbr_line==2):
				line2 = line

		if (nbr_line==2 and len(line2)-1 >= minLengthReads):
			writeRead = True

		if (nbr_line==4):
			line = line[0:(len(line)-after)]
			line4 = line.replace("\n","")+"\n"
			nbr_line = 0
			after = 0

			if (writeRead == True):
				outputFile.write(line1 + line2 + line3 + line4)
				minReads += 1
			writeRead = False

	print "On", totalReads, "total reads,", polyA_counter, "reads have a polyA tail removed."
	print "After removed polyA tail, there are", (totalReads-minReads), "too short reads which have been removed." 

	inputFile.close()
	outputFile.close()


