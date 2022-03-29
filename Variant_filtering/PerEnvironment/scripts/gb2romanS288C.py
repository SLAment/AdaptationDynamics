#!/usr/bin/env python
# encoding: utf-8

# ================== gb2romanS288C.py =================
# Change chromosome names in a vcf file from the Genbank accession numbers to
# the Roman numbers to match the database in snpEff
# ==================================================
# Sandra Lorena Ament Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021/10/27
# +++++++++++++++++++++++++++++++++++++++++++++++++

import argparse # For the fancy options
import sys  # To exit the script, and to pipe out
import re
from datetime import datetime

# ------------------------------------------------------
version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* gb2romanS288C v." + versiondisplay + " *", epilog="Change chromosome names in a vcf file from the Genbank accession numbers to\nthe Roman numbers to match the database in snpEff.") # Create the object using class argparse

# Add basic options
parser.add_argument('vcf', help="Standard vcf file (not compressed")
parser.add_argument('--rom2gb', '-r', help="Change roman numbers to Genbank (default is the inverse)", default=False, action='store_true')
parser.add_argument('--chr', '-c', help="The chromosome names have a format like 'chrX' rather than 'X' ", default=False, action='store_true')

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)


try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	vcfopen = open(args.vcf, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()
# ============================

# ------------------------------------------------------
# Prepare functions
# ------------------------------------------------------

if args.rom2gb:
	if args.chr:
		chrdic = {"chrI": "BK006935.2", "chrII": "BK006936.2", "chrIII": "BK006937.2", "chrIV": "BK006938.2", "chrV": "BK006939.2", "chrVI": "BK006940.2", "chrVII": "BK006941.2", "chrVIII": "BK006934.2", "chrIX": "BK006942.2", "chrX": "BK006943.2", "chrXI": "BK006944.2", "chrXII": "BK006945.2", "chrXIII": "BK006946.2", "chrXIV": "BK006947.3", "chrXV": "BK006948.2", "chrXVI": "BK006949.2", "chrmt": "AJ011856.1"}
	else:
		chrdic = {"I": "BK006935.2", "II": "BK006936.2", "III": "BK006937.2", "IV": "BK006938.2", "V": "BK006939.2", "VI": "BK006940.2", "VII": "BK006941.2", "VIII": "BK006934.2", "IX": "BK006942.2", "X": "BK006943.2", "XI": "BK006944.2", "XII": "BK006945.2", "XIII": "BK006946.2", "XIV": "BK006947.3", "XV": "BK006948.2", "XVI": "BK006949.2", "Mito": "AJ011856.1"}
else:
	chrdic = {"BK006935.2": "I", "BK006936.2": "II", "BK006937.2": "III", "BK006938.2": "IV", "BK006939.2": "V", "BK006940.2": "VI", "BK006941.2": "VII", "BK006934.2": "VIII", "BK006942.2": "IX", "BK006943.2": "X", "BK006944.2": "XI", "BK006945.2": "XII", "BK006946.2": "XIII", "BK006947.3": "XIV", "BK006948.2": "XV", "BK006949.2": "XVI", "AJ011856.1": "Mito"}


chrheadregex = re.compile("(##contig=<ID=)([\w\.]+)(.*)")

# ------------------------------------------------------
# Finally filter
# ------------------------------------------------------

for line in vcfopen:
	if "##contig" in line:
		chrhead = chrheadregex.search(line)
		if chrhead:
			newline = f"##contig=<ID={chrdic[chrhead.group(2)]}{chrhead.group(3)}\n" 
			sys.stdout.write(newline) 
		else:
			sys.stdout.write(line) 
	elif "##" in line: # header
		sys.stdout.write(line) 
		pass	
	elif "#CHROM" in line: # column names
		now = datetime.now()
		sys.stdout.write('##gb2romanS288C.py=v.' + versiondisplay + ' on ' + now.strftime("%d/%m/%Y %H:%M:%S") + '\n')
		sys.stdout.write(line) 
	else:
		# Extract information from the line
		columns = line.rstrip("\n").split('\t')
		CHROM = columns[0]

		newline = f"{chrdic[CHROM]}\t" + '\t'.join(columns[1:]) + "\n"
		sys.stdout.write(newline)


