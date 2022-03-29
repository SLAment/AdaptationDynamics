#!/usr/bin/env python
# encoding: utf-8

# ================== getvariantspool.py : Filter the vcf file using coverage values =================
# ==================================================
# Sandra Lorena Ament Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021/09/28
# +++++++++++++++++++++++++++++++++++++++++++++++++


import argparse # For the fancy options
import sys  # To exit the script, and to pipe out
import gzip # For unzipping files
import os  # For the input name
from datetime import datetime

# ------------------------------------------------------
version = 1.0
versiondisplay = "{0:.2f}".format(version)

# Input from console
try:
	vcffile = sys.argv[1]
except:
	print("Usage: python " + sys.argv[0] + " myfile.vcf.gz")
	print("Version " + versiondisplay)
	sys.exit(1)

vcfopen = gzip.open(vcffile, 'rt') # 't' is text mode, to interpret the tabs and new lines
# vcfopen = open(vcffile, 'r')

# ------------------------------------------------------
# Define functions
# ------------------------------------------------------
def getcov( covstring ):
	if covstring == "./.:0,0": # Some crappy sites have this format
		return(0)
	else:
		alleles = [feature for feature in covstring.split(":")]	
		if alleles[1] == ".": # No read is valid
			return(0)
		else:
			AD1, AD2 = [float(feature) for feature in alleles[1].split(",")]	 
			cov = AD1 + AD2
			if cov == 0:
				return(0)
			else:
				maf = min(AD1, AD2)/cov
				return(maf)

# ------------------------------------------------------
# Finally filter
# ------------------------------------------------------

for line in vcfopen:
	if "##" in line: # header
		sys.stdout.write(line) 
		# pass	
	elif "#CHROM" in line: # column names
		# pass
		# sys.stdout.write('##FILTER=<ID=getvariantspool,Description="Sites surviving filtering based on coverage">\n') # TODO: make this more informative

		now = datetime.now()
		sys.stdout.write('##getvariantspool.py=v. ' + versiondisplay + ' on ' + now.strftime("%d/%m/%Y %H:%M:%S") + '\n')

		sys.stdout.write(line) 
	else:
		keep = 0 # Can this site survive the filtering?

		# Extract information from the line
		columns = line.rstrip("\n").split('\t')
		site = columns[9:]

		for i in range(0, len(site)): # Sum all the MAF of all samples. If they are all 0, then this site is monomorphic
			cov = getcov(site[i])
			keep += cov
		
		# Print the surviving sites
		if keep > 0: 
			sys.stdout.write(line)
