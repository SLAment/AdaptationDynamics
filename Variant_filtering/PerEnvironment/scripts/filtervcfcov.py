#!/usr/bin/env python
# encoding: utf-8

# ================== filtervcfcov: Filter the vcf file using coverage values =================
# ==================================================
# Sandra Lorena Ament Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021/09/27
# +++++++++++++++++++++++++++++++++++++++++++++++++


import argparse # For the fancy options
import sys  # To exit the script, and to pipe out
import os  # For the input name
from datetime import datetime

# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Filter the vcf file using coverage values *", epilog="") # Create the object using class argparse

# Add options
parser.add_argument('vcf', help="Standard vcf file compressed with bgzip")
parser.add_argument('--coveragesummary', '-s', help="A coverage summary produced by vcf4adaptation_vcfR_plotter.R with samples as rows and culumns called q0.25 and q0.95 for filtering.")
parser.add_argument("--min-coverage", "-m", help="Sites with less than this value gets excluded. If used, then coveragesummary q0.25 is ignored.", dest="min", type=int) # nargs='+' All, and at least one, argument
parser.add_argument("--max-coverage", "-M", help="Sites with more than this value gets excluded. If used, then coveragesummary q0.95 is ignored.", dest="max", type=int) # nargs='+' All, and at least one, argument
parser.add_argument("--basemin", "-b", help="If set, sites will be filtered out if their coverage is smaller than max(q0.025, BASEMIN)", type=int) # nargs='+' All, and at least one, argument

# More options
parser.add_argument("--ignorecols", "-c", help="String of sample column numbers in base 0 to be ignored, separated by commas and starting at 0. E.g. 5,6,7.", type=str) # nargs='+' All, and at least one, argument
parser.add_argument("--ignoremedian", "-i", help="Ignore this sample during filtering if its median is below this value. It only works with --coveragesummary.", type=float) # nargs='+' All, and at least one, argument

# parser.add_argument('--uncompressed', '-z', help="VCF file is not compressed", default=False, action='store_true')

try:
	# ArgumentParser parses arguments through the parse_args() method. You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	vcfopen = open(args.vcf, 'r')

except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()

# ------------------------------------------------------
# Define functions
# ------------------------------------------------------

def getcov( covstring ):
	if covstring == "./.:0,0": # Some crappy sites have this format
		return(0)
	else:
		# alleles = [feature for feature in covstring.split(":")]	
		# DP = alleles[2] # Approximate read depth (reads with MQ=255 or with bad mates are filtered)

		# if DP == ".": # No read is valid
		# 	return(0)
		# else:
		# 	return(int(DP))
		alleles = [feature for feature in covstring.split(":")]	
		if alleles[1] == ".": # No read is valid
			return(0)
		else:
			AD1, AD2 = [float(feature) for feature in alleles[1].split(",")]	 
			DP = AD1 + AD2
			return(DP)

# ------------------------------------------------------
# Massage the data
# ------------------------------------------------------

if args.coveragesummary: # There is an input file with individual coverage min and max values to use for filtering
	coveragesummary = open(args.coveragesummary, 'r')
	covsum = [line.rstrip("\n").split("\t") for line in coveragesummary]

	# Make a dictionary with the coverage minimum and maximum of each sample
	samplesdic = {}
	for line in covsum[1:]: # the first line is the header, so skip
		samplesdic[line[0]] = (float(line[3]), float(line[4]), float(line[1])) # Min, Max, and median

	# samples = list(samplesdic.keys())

elif args.min == None and args.max == None: # There is nothing to filter with
	print("No filtering minimum (--min-coverage) or maximum coverage (--max-coverage) set for filtering, and no CoverageSummary.txt file given (--coveragesummary).")
	sys.exit(1)
else: # There is only one parameter to filter
	if args.min != None: # The min was set, but not the max
		maxcov = 10000000000000 # A crazy large number that won't be surpassed ever
	elif args.max != None: # The max was set, but not the min
		mincov = 0


# Are there samples to be ignored?
badcolums = []
if args.ignorecols != None: 
	badcolums = [int(col) for col in args.ignorecols.split(",")] # Make list of columns

# ------------------------------------------------------
# Finally filter
# ------------------------------------------------------

for line in vcfopen:
	if "##" in line: # header
		sys.stdout.write(line) 
		# pass	
	elif "#CHROM" in line: # column names
		# Get samples 
		columns = line.rstrip("\n").split('\t')
		samples = [i for i in columns[9:]]
		# pass
		sys.stdout.write('##FILTER=<ID=filtervcfcov,Description="Sites surviving filtering based on coverage">\n') # TODO: make this more informative

		now = datetime.now()
		sys.stdout.write('##filtervcfcov.py=v. ' + versiondisplay + ' on ' + now.strftime("%d/%m/%Y %H:%M:%S") + '\n')

		sys.stdout.write(line) 
	else:
		keep = True # Can this site survive the filtering?

		# Extract information from the line
		columns = line.rstrip("\n").split('\t')
		site = columns[9:]

		for i in range(0, len(site)):
			if i not in badcolums: # If it's not in the list, then check the coverage
				cov = getcov(site[i])

				# Is the coverage of this site within the ranges for this sample?
				if args.coveragesummary:
					mincov, maxcov, median = samplesdic[samples[i]]

				# Redefined the min and max coverage for this sample if user asks for it
				if args.min != None:
					mincov = args.min

				if args.max != None:
					maxcov = args.max

				if args.basemin != None: # Use the biggest number, the q25/set minimum or the BASEMIN value	
					if args.basemin > mincov:	
						mincov = args.basemin

				# Filter
				if (cov < mincov) or (cov > maxcov):
					if args.coveragesummary and args.ignoremedian != None: # Are you desperate and want to ignore samples with low median coverage?
						if median >= args.ignoremedian: # This sample is still good enough to filter
							keep = False # Remove this site
					else:
						keep = False # Remove this site
				# print(samples[i], cov, mincov, maxcov, keep, median)
		# print(line)
		# sys.exit(1)

		# Print the surviving sites
		if keep: sys.stdout.write(line) 
