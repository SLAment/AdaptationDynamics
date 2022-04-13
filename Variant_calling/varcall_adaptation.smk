# -*- snakemake -*-

### varcall_adaptation.smk: Calling variants for the adaptation experiment
#############################################################################

# https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected

# Nima Rafati had an early version of this pipeline, so I have some notes here and there from that

#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021-10-08
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

from glob import glob
import os


# module load bioinfo-tools snakemake/5.30.1 GATK/4.1.4.1 samtools/1.9

# -------------------------------------------------
# Data
PATH2gvcfs = "path/2/gvcfs" # CHANGE
REFGenome = "path/2/reference/genome.fa" # CHANGE
REFGenomeBED = "path/2/reference/genome.fa.fai.bed" # CHANGE

# Global variables
ENVIRONMENTS = ["NaCl", "Ethanol", "LiAc0.01", "LiAc0.02", "SK1", "Y55", "N_Founder_R1", "LE_Founder_R1"]

# Samples to be removed
BADLIST = [
"NaCl_G1000_R4",
"NaCl_G1000_R5",
"NaCl_G30_R5",
"NaCl_G60_R5",
"NaCl_G100_R5",
"NaCl_G200_R5",
"NaCl_G300_R5",
"NaCl_G400_R5",
"NaCl_G500_R5",
"NaCl_G700_R5",
"LiAc0.02_G60_R2",
"NaCl_G1000_R1",
"LiAc0.01_G30_R3",
"LiAc0.01_G60_R3",
"LiAc0.01_G100_R3",
"LiAc0.01_G200_R3",
"LiAc0.01_G300_R3",
"LiAc0.01_G400_R3",
"LiAc0.01_G500_R3",
"LiAc0.01_G700_R3",
"LiAc0.01_G1000_R3",
"LiAc0.01_G1000_R1",
"LiAc0.01_G1000_R2",
"LiAc0.01_G1000_R4",
"LiAc0.01_G1000_R5",
"LiAc0.02_G1000_R1",
"LiAc0.02_G1000_R2",
"LiAc0.02_G1000_R3",
"LiAc0.02_G1000_R4",
"LiAc0.02_G1000_R5"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: getPASSsites
# ----------


## Make a list of the environments (key) and their gvcfs (values)
gvcfls = [(glob(PATH2gvcfs + "/*/{environment}*g.vcf.gz".format(environment=environment))) for environment in ENVIRONMENTS]

# That created a list of lists for each environment, so flatten the list of gvcf files
# https://datascienceparichay.com/article/python-flatten-a-list-of-lists-to-a-single-list/
flat_gvcfls = [item for sublist in gvcfls for item in sublist]

# Some samples are problematic, so remove them
good_gvcfs = []
for gvcf in flat_gvcfls:
	good = True
	for bad in BADLIST:
		if bad in gvcf:
			good = False
	if good: good_gvcfs.append(gvcf)

rule all:
	input:
		"All_Variants_parents_PASS.vcf",
		expand("Raw_Call_parents_split_multi_{typevar}_FILTER.vcf.gz", typevar = ["SNP", "INDEL"]),


rule GenomicsDBImport: # Equivalent to Nima's Sbatch_GenomicsDB.script
	""" Import single-sample GVCFs into GenomicsDB before joint genotyping """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.1/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php
	# https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport
	input:
		intervals = REFGenomeBED,
		gvcfs = good_gvcfs
	output:
		directory("genomicsdb")
	params:
		time = "3:00:00",
		threads = 2,
		JavaMem = int(2 * 6), # A Rackham node contains 128 GB of RAM and 20 compute cores (each core gets at most 6.8 GB).
	run:
		# Create a string in the format --variant path/to/gvcf/sample1 --variant path/to/gvcf/sample2 etc...
		variantlist = ""
		for sample in input.gvcfs:
			variantlist += "--variant " + sample + " "
		
		shell("mkdir -p tmp")
		# Notice GATK will create the output directory
		gatkcommand = f'gatk GenomicsDBImport --java-options "-Xmx{params.JavaMem}G" --tmp-dir=tmp --genomicsdb-workspace-path {output[0]} -L {input.intervals} {variantlist}'
		shell(gatkcommand) # execute
		

rule GenotypeGVCFs: # Equivalent to Nima's Sbatch_Genotyping.script
	""" Perform joint genotyping """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php
	# https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport
	input:
		my_database = "genomicsdb",
		ref = REFGenome
	output:
		rawvcf = "Raw_Call_parents.vcf"
	params:
		time = "2-00:00:00",
		threads = 8,
		JavaMem = int(8 * 6.8),
	shell:
		"""
		gatk --java-options "-Xmx{params.JavaMem}G" GenotypeGVCFs \\
		-R {input.ref} \\
		-V gendb://{input.my_database} \\
		-O {output.rawvcf} \\
		-G StandardAnnotation --new-qual 
		"""
# --create-output-variant-index	If true, create a VCF index when writing a coordinate-sorted VCF file.
# --use-new-qual-calculator / -new-qual 	Use the new AF model instead of the so-called exact model. Default: true
# By default, GATK HaplotypeCaller and GenotypeGVCFs do not emit variants with QUAL < 10, controlled with -stand-call-conf

# Nima added 
# -G StandardAnnotation
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622
# https://gatk.broadinstitute.org/hc/en-us/articles/360050814612-HaplotypeCaller
# --annotation-group, -G:String  One or more groups of annotations to apply to variant calls. This argument may be specified 0 or more times. Default value: null. Possible Values: {AlleleSpecificAnnotation, AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation, StandardHCAnnotation, StandardMutectAnnotation}

# Note: This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy for non-diploid organisms.

rule SplitMultiallelic: # Part of Sbatch_Extracting_SNP_INDEL_and_stats.script
	""" Following Nima, split multiallelic sites into biallelic """
	input:
		ref = REFGenome,
		vcf = "Raw_Call_parents.vcf"
	output:
		vcf = "Raw_Call_parents_split_multi.vcf"
	params:
		time = "4:00:00",
		threads = 2,
	shell:
		"gatk LeftAlignAndTrimVariants --split-multi-allelics --dont-trim-alleles -R {input.ref} -V {input.vcf} -O {output.vcf}"

rule SplitVariants: # Part of Sbatch_Extracting_SNP_INDEL_and_stats.script
	""" Separate the SNP/INDEL """
	input:
		ref = REFGenome,
		vcf = "Raw_Call_parents_split_multi.vcf"
	output:
		vcf = "Raw_Call_parents_split_multi_{vartype}.vcf"
	params:
		time = "1:00:00",
		threads = 1,
	shell:
		"gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include {wildcards.vartype} --restrict-alleles-to BIALLELIC -O {output.vcf}"

		# https://gatk.broadinstitute.org/hc/en-us/articles/360036362532-SelectVariants
		# Could it be Nima used --restrict-alleles-to BIALLELIC because the multiallelic 

rule MarkFilter_SNP: # Part of Sbatch_Filtering_Manually.script
	""" Mark sites that pass the filter using Nima's predetermined values """
	# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
	input:
		vcf = "Raw_Call_parents_split_multi_SNP.vcf"
	output:
		vcf = "Raw_Call_parents_split_multi_SNP_FILTER.vcf"
	params:
		time = "30:00",
		threads = 1, # for memory?
	shell:
		"""
		gatk VariantFiltration \\
		--variant {input.vcf} \\
		-O {output.vcf} \\
		-filter "QD < 2.0" --filter-name "QD2" \\
		-filter "FS > 10.0" --filter-name "FS10" \\
		-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3" \\
		-filter "MQRankSum  < -6.0 " --filter-name "MQRankSum-6" \\
		-filter "SOR > 3.0" --filter-name "SOR-3" \\
		-filter "MQ < 40.0" --filter-name "MQ-40" 
		"""
		# Nima had "MQ > 40.0" but I believe that is wrong (see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)

rule MarkFilter_INDEL: # Part of Sbatch_Filtering_Manually.script
	""" Mark sites that pass the filter using Nima's predetermined values """
	input:
		vcf = "Raw_Call_parents_split_multi_INDEL.vcf"
	output:
		vcf = "Raw_Call_parents_split_multi_INDEL_FILTER.vcf"
	params:
		time = "30:00",
		threads = 1, # for memory?
	shell:
		"""
		gatk VariantFiltration \\
		--variant {input.vcf} \\
		-O {output.vcf} \\
		-filter "QD < 10.0" --filter-name "QD10" \\
		-filter "FS > 10.0" --filter-name "FS10" \\
		-filter "ReadPosRankSum < -4.0" --filter-name "ReadPosRankSum-4" \\
		-filter "MQRankSum  < -8.0 " --filter-name "MQRankSum-8" \\
		-filter "SOR > 4.0" --filter-name "SOR-4" \\
		-filter "MQ < 40.0" --filter-name "MQ-40"
		"""

rule getPASSsites: # Part of Sbatch_Filtering_Manually.script
	""" Remove sites rejected in VariantFiltration and put them together """
	input:
		snps = "Raw_Call_parents_split_multi_SNP_FILTER.vcf",
		indels = "Raw_Call_parents_split_multi_INDEL_FILTER.vcf",
	output:
		vcf = "All_Variants_parents_PASS.vcf"
	shell:
		"""
		# Extract PASS variants and save them all in one single file
		nline=$(grep -m1 -n CHROM {input.snps} | sed 's/:.*//')
		head -$nline {input.snps} > {output.vcf}
		awk '($7 == "PASS")' {input.snps} {input.indels} | sort -k1,1 -k2n,2 >> {output.vcf}
		"""

rule compressAndIndex:
	input:
		bigvcf = "All_Variants_parents_PASS.vcf", # Just so it does this after making the big file
		vcf = "Raw_Call_parents_split_multi_{typevar}_FILTER.vcf"
	output:
		vcf = "Raw_Call_parents_split_multi_{typevar}_FILTER.vcf.gz"
	params:
		time = "4:00:00",
		threads = 1
	shell:
		"bgzip {input.vcf} && tabix -p vcf {output.vcf}"

