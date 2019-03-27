#!/usr/bin/python2.7
#
# subset_fastq.py
#
# Pull out the reads you want
# 
# Z Norgaard
#
# Last Modified: 10nov2015
#

# Arguments
from sys import argv
script, infile, outfile = argv

# Modules
from decimal import *
from cmath import *

# Open Out-File
out = open(outfile, 'w')

# Parse In-File
with open(infile) as fh:
	while True:
		# Store Info, End of File Check
		hdr = fh.readline()
		if not hdr: break
		read = fh.readline()
		spcr = fh.readline()
		qscor = fh.readline()

		# Calculate Average Read Score
		splitscor = [ord(c) for c in qscor.rstrip('\n')]
		splitscor = [c - 33 for c in splitscor]
		read_len = len(splitscor)
		avgqscor = sum(splitscor)/Decimal(read_len)
		
		if avgqscor >= 30:
			out.write('%s%s%s%s' % (hdr, read, spcr, qscor))

out.close()
