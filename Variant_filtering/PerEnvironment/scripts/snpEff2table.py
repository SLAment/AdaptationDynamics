#!/usr/bin/env python
# encoding: utf-8

# ================== snpEff2table.py =================
# Extract the first and major effect of the annotation of SnpEff in a vcf file
# ==================================================
# Sandra Lorena Ament Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021/11/03
# +++++++++++++++++++++++++++++++++++++++++++++++++

import sys  # To exit the script, and to pipe out

# ------------------------------------------------------
version = 1.0
versiondisplay = "{0:.2f}".format(version)

# Input from console
try:
	vcffile = sys.argv[1]
except:
	print("Usage: python " + sys.argv[0] + " myfile.vcf")
	print("Version " + versiondisplay)
	sys.exit(1)

# vcfopen = gzip.open(vcffile, 'rt') # 't' is text mode, to interpret the tabs and new lines
vcfopen = open(vcffile, 'r')

# ------------------------------------------------------
# Read through the file
# ------------------------------------------------------

for line in vcfopen:
	if "##" in line: # header
		pass
	elif "#CHROM" in line: # column names
		head = "Contig\tposition\tALT_allele\ttype\tputative_effect\tgene_name\tGene_ID\tFeature_type\tFeature_ID\tTranscript_biotype\tRank_total\tHGVS.c\tHGVS.p\tcDNA_pos_len\tCDS_pos_len\tProtein_pos_len\tDistance_to_feature\tNotes\n"
		sys.stdout.write(head) 
	else:
		# Extract information from the line
		columns = line.rstrip("\n").split('\t')
		CONTIG = columns[0]
		POS = columns[1]
		INFO = columns[7].split(';')

		for inf in INFO:
			if 'ANN=' in inf:
				rawannotation = inf

		annotation = rawannotation.replace('ANN=', '').split(',')

		# For now let's just report the first annotation
		focalannotation = annotation[0].split('|')
		
		focalannotation_string = ''
		for element in focalannotation:
			if element == '':
				focalannotation_string += '\tNA'
			else:
				focalannotation_string += f"\t{element}"

		# focalannotation_string = '\t'.join(focalannotation)

		newline = f'{CONTIG}\t{POS}\t{focalannotation_string}\n'
		sys.stdout.write(newline) 



