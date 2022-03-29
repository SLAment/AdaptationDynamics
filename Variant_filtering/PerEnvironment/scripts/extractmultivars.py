#!/usr/bin/env python
# encoding: utf-8

# ================== extractmultivars.py =================
# Separate multiallelic sites (after splitting with GATK as 
# `LeftAlignAndTrimVariants --split-multi-allelics`) from normal 
# biallelic sites in different files
# ==================================================
# Sandra Lorena Ament Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021/10/29
# +++++++++++++++++++++++++++++++++++++++++++++++++

import argparse # For the fancy options
import sys  # To exit the script, and to pipe out
import gzip # For unzipping files
import os  # For the input name
from datetime import datetime
from collections import defaultdict #This tells python that the dictionary contains a list so you can freely append things to the value of each key

# ------------------------------------------------------
version = 1.2
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Filter the vcf file using coverage values *", epilog="") # Create the object using class argparse

# Add options
parser.add_argument('vcf', help="Standard vcf file")
parser.add_argument("--outputdir", "-o", help="Path to output files (default is folder of the vcf file)")
# parser.add_argument('--uncompressed', '-z', help="VCF file is not compressed", default=False, action='store_true')

try:
	# ArgumentParser parses arguments through the parse_args() method. You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	if '.gz' in args.vcf:
		vcfopen = gzip.open(args.vcf, 'rt') # 't' is text mode, to interpret the tabs and new lines
	else:
		vcfopen = open(args.vcf, 'r')

except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()

# ------------------------------------------------------
# Define output
# ------------------------------------------------------

if '.gz' in args.vcf:
	filename = args.vcf.replace(".vcf.gz", "")
else:
	filename = args.vcf.replace(".vcf", "") # Taking out the prefix of the file

if args.outputdir:
	filename = f"{args.outputdir}/{os.path.basename(filename)}"

outbis = open(f"{filename}_bi.vcf", "w")
outmultis = open(f"{filename}_multi.vcf", "w")


# ------------------------------------------------------
# Try to get sites
# ------------------------------------------------------
sites = defaultdict(list) #This tells python that the dictionary contains a list so you can freely append things to the value of each key
for line in vcfopen:
	if "##" in line: # header
		outbis.write(line)
		outmultis.write(line)
	elif "#CHROM" in line: # column names
		now = datetime.now()

		outbis.write('##extractmultivars.py=v. ' + versiondisplay + ' to get biallelic sites on ' + now.strftime("%d/%m/%Y %H:%M:%S") + '\n')
		outmultis.write('##extractmultivars.py=v. ' + versiondisplay + ' to get split multiallelic sites on ' + now.strftime("%d/%m/%Y %H:%M:%S") + '\n')

		outbis.write(line)
		outmultis.write(line)
	else:
		# I need to save the sites in memory to identify the multiallelic ones
		# Extract information from the line
		columns = line.rstrip("\n").split('\t')
		columns = line.rstrip("\n").split('\t')
		coordinate = columns[0] + "_" + columns[1]
		
		sites[coordinate].append(line)

countmulti = 0
countibiall = 0
for site in sites.keys():
	if len(sites[site]) == 1:
		countibiall += 1
		outbis.write(sites[site][0])
	else:
		countmulti += 1
		for mul in sites[site]:
			outmultis.write(mul)

print(f"There are {countibiall} biallelic sites.")
print(f"There are {countmulti} multiallelic sites.")
