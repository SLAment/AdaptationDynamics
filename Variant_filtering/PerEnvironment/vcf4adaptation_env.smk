# -*- snakemake -*-

### vcf4adaptation_env: Preparing the vcf file for the analysis of allele frequencies, splitting by environment
#############################################################################

# The variant calling was done previously by Nima Rafati, including quality
# filtering.
# He aso modified the vcf file to split multiallelic sites in such a way that
# they are represented in more than one line. For example:

# $ grep -P 'BK006934.2\t1495\t' Raw_Call-SNP.vcf | cut -f1-10
# BK006934.2	1495	.	G	A	32531.97	.	AC=166;AF=0.241;AN=690	GT:AD:DP:GQ:PL	0/1:92,12:109:99:136,0,3034
# BK006934.2	1495	.	G	C	32531.97	.	AC=12;AF=0.017;AN=690	GT:AD:DP:GQ:PL	0/0:92,5:109:99:0,150,3197

# In the original file, the site looks like so:

# $ grep -P 'BK006934.2\t1495\t' Raw_Call.vcf | cut -f1-10
# BK006934.2	1495	.	G	A,C	32531.97	.	AC=166,12;AF=0.241,0.017;AN=690;BaseQRankSum=0.00;DP=38775;ExcessHet=139.8646;FS=37.983;InbreedingCoeff=-0.3996;MLEAC=183,10;MLEAF=0.265,0.014;MQ=33.44;MQRankSum=0.00;QD=1.89;ReadPosRankSum=0.00;SOR=4.220	GT:AD:DP:GQ:PL	0/1:92,12,5:109:99:136,0,3034,286,2915,3333

#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021-09-24-29
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# $ module load bioinfo-tools snakemake/5.30.1 bcftools/1.12 vcftools/0.1.16 GATK/4.1.4.1 
import os

# -------------------------------------------------
# Data
BASEVCF = "path/to/All_Variants_parents_PASS.vcf"
REFGenome = "data/genome.fa"
bad_mutations = "data/Manually_curated_bad_mutations.txt"

# Scripts to massage the data
vcf4adaptation_vcfR_plotter = "scripts/vcf4adaptation_vcfR_plotter_env.R"
filtervcfcov = "scripts/filtervcfcov.py"
getvariantspool = "scripts/getvariantspool.py"
extractmultivars = "scripts/extractmultivars.py"
snpEff2table = "scripts/snpEff2table.py"

# Plotting
Y55_along_contigs = "scripts/Y55_along_contigs.R"
gb2romanS288C = "scripts/gb2romanS288C.py"
DeNovoMutations = "scripts/DeNovoMutations.R"

# Global variables
ENVIRONMENTS = ["NaCl", "Ethanol", "LiAc0.01", "LiAc0.02"]

# -------------------------------------------------
# Define the output name from the input vcf file
# outputname = os.path.basename(BASEVCF).rstrip('.vcf')
outputname = os.path.splitext(os.path.basename(BASEVCF))[0]

# ----------
# Rules not submitted to a job
localrules: definesamples, listSamplesByEnv, splitvcf_by_allelic, changecontignames, changecontigback, getSnpEffTable, genmapindex, genmap, filtergenmap, bedtools
# ----------

rule all:
	input:
		# The results of de novo mutations
		"results/deNovo_all_0.1_CoolGenes.pdf",
		"results/deNovoFixers_0.35_genes_4stringdb.txt",
		"results/deNovoFixers_0.35_SNPeff.pdf",
		
		# Mention these two, so they both go through the same rule
		expand("figures/{outputname}_{env}_var_bi_miss0_{typevar}_Coverage.png", env = ENVIRONMENTS, outputname = outputname, typevar = ["SNPs", "INDELs"]),
		expand("figures/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_map1_Coverage.png", env = ENVIRONMENTS, outputname = outputname, typevar = ["SNPs", "INDELs"]),

		# Just to make sure they both get compressed
		expand("vcfs/{outputname}_{env}_var_{allelic}.vcf.gz", env = ENVIRONMENTS, outputname = outputname, allelic = ["bi", "multi"]), 
		

rule definesamples:
	""" Get a list of samples to extract from the original vcf file """
	input:
		BASEVCF
	output:
		"data/samples.txt"
	shell:
		"bcftools view -h {input} | grep '#CHROM' | cut -f10- | tr '\\t' '\\n' > {output}"
		# "bcftools view -h {input} | grep '#CHROM' | cut -f10- | tr '\\t' '\\n' | grep -vE '^30C_T|^H' | grep -vE '^NaCl_G[0-9]+_R5|NaCl_G1000_R4|LiAc0.01_G1000_R[1,2,4,5]|LiAc0.02_G60_R2' > {output}"
		# Exclude the hibryd samples, contaminated treatments and time points, and low coverage ones


rule listSamplesByEnv:
	""" Make lists of samples by environment """
	input:
		list = "data/samples.txt"
	output:
		list = "data/samples_{env}.txt"
	run:
		if wildcards.env == "NaCl":
			founder = "N_Founder"
		else:
			founder = "LE_Founder"
		cmd = "grep -E '^" + f"{wildcards.env}|^{founder}|^SK1|^Y55" + "' {input.list} > {output.list}"
		shell(cmd)

rule extractsamplesfromvcf:
	""" Extract samples from the vcf file """
	input:
		vcf = BASEVCF,
		lista = "data/samples_{env}.txt"
	output:
		"vcfs/{outputname}_{env}.vcf.gz"
	wildcard_constraints:
		outputname = outputname,
		env = '|'.join([x for x in ENVIRONMENTS])
	params:
		# threads = 1,
		time = "1:00:00",
	shell:
		"vcftools --vcf {input.vcf} --recode --recode-INFO-all --stdout --keep {input.lista} | bgzip > {output}"

rule filter4polymorphic:
	""" Remove sites that are no longer polymorphic (approx) """
	input:
		vcf = "vcfs/{outputname}_{env}.vcf.gz"
	output:
		vcf = "vcfs/{outputname}_{env}_var.vcf.gz"
	params:
		# threads = 1,
		time = "1:00:00",
	shell:
		"python {getvariantspool} {input.vcf} | bgzip > {output.vcf}"
	# I wrote my own script to remove invarible sites based on the alelle frequencies, not the genotypes of GATK

	# 	"vcftools --gzvcf {input} --recode --recode-INFO-all --stdout --maf 0.0005 | bgzip > {output}"
	# 	# --maf <float>: Include only sites with a Minor Allele Frequency greater than or equal to the "--maf"
	# 	# 0.0005 -> 1 read for 2000x
	# 	# Using --mac 1 would have been better but it's the same result
	# 	# --mac <integer>: Include only sites with Minor Allele Count greater than or equal to the "--mac" value

rule splitvcf_by_allelic:
	""" Split the biallelic sites from multiallelic sites in different files """
	input:
		"vcfs/{outputname}_{env}_var.vcf.gz"
	output:
		bi = "vcfs/{outputname}_{env}_var_bi.vcf",
		multi = "vcfs/{outputname}_{env}_var_multi.vcf", # We didn't use them in the end
		log = "logs/{outputname}_{env}_extractmultivars.txt"
	shell:
		"python {extractmultivars} {input} > {output.log}"

rule compressvcfs:
	input:
		"vcfs/{outputname}_{env}_var_{allelic}.vcf"
	output:
		"vcfs/{outputname}_{env}_var_{allelic}.vcf.gz"
	wildcard_constraints:
		outputname = outputname,
		env = '|'.join([x for x in ENVIRONMENTS]),
		allelic = 'bi|multi'
	params:
		time = "20:00",
		# threads = 1,
	shell:
		"bgzip {input}"

rule filter4missingdata:
	""" Remove sites with missing data """
	input:
		"vcfs/{outputname}_{env}_var_bi.vcf.gz"
	output:
		"vcfs/{outputname}_{env}_var_bi_miss0.vcf.gz"
	params:
		# threads = 1,
		time = "1:00:00",
	shell:
		"vcftools --gzvcf {input} --recode --recode-INFO-all --stdout --max-missing 1 | bgzip > {output}"
		# --max-missing <float>: Exclude sites on the basis of the proportion
		#   of missing data (defined to be between 0 and 1, where 0 allows
		#   sites that are completely missing and 1 indicates no missing data
		#   allowed).

rule makeindex_vcf:
	""" Index vcf file """
	input:
		"vcfs/{outputname}_{env}_var_bi_miss0.vcf.gz"
	output:
		"vcfs/{outputname}_{env}_var_bi_miss0.vcf.gz.tbi"
	params:
		# threads = 1,
		time = "30:00",
	shell:
		"tabix -p vcf {input}"

rule getSNPs:
	""" Make separate vcf file for just the SNPs """
	input:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0.vcf.gz",
		index = "vcfs/{outputname}_{env}_var_bi_miss0.vcf.gz.tbi"
	output:
		"vcfs/{outputname}_{env}_var_bi_miss0_SNPs.vcf"
	params:
		# threads = 1,
		time = "1:00:00",		
	shell:
		"gatk SelectVariants -R {REFGenome} -V {input.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {output}"
		# "gatk SelectVariants -R {REFGenome} -V {input.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {output}"
		# This however will still keep sites where the alternative is "*", which means a deletion of that base.
		# Nima had --restrict-alleles-to BIALLELIC but because the sites are allready split, I don't think it works as intended
		# https://gatk.broadinstitute.org/hc/en-us/articles/360037260511-SelectVariants#--restrict-alleles-to

rule getINDELs:
	""" Make separate vcf file for just the INDELs """
	input:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0.vcf.gz",
		index = "vcfs/{outputname}_{env}_var_bi_miss0.vcf.gz.tbi"
	output:
		"vcfs/{outputname}_{env}_var_bi_miss0_INDELs.vcf"
	params:
		# threads = 1,
		time = "1:00:00",		
	shell:
		"gatk SelectVariants -R {REFGenome} -V {input.vcf} --select-type-to-include INDEL --restrict-alleles-to BIALLELIC -O {output}"

# ------- Find sites that are not good -------

rule PlotCoverage:
	""" Plot coverage distribution for all samples and report """
	input:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}.vcf"
	output:
		plot = "figures/{outputname}_{env}_var_bi_miss0_{typevar}_Coverage.png",
		table = "data/{outputname}_{env}_var_bi_miss0_{typevar}_CoverageSummary.tab"
	params:
		time = "10:00:00",
		threads = 10, # It needs some memory
		lquantile = 0.25,
		uquantile = 0.95,
	# conda: 
	# 	"envs/statsplot.yaml"
	script:
		vcf4adaptation_vcfR_plotter

# ------- Filter -------

rule filtervcfcov:
	""" Filter vcf file based on the coverage distribution per sample """
	input:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}.vcf",
		table = "data/{outputname}_{env}_var_bi_miss0_{typevar}_CoverageSummary.tab"
	output:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p.vcf",
	params:
		time = "3:00:00",
		# threads = 1,
	shell:
		"python {filtervcfcov} {input.vcf} --coveragesummary {input.table} --min-coverage 25 > {output.vcf}"
		# "python {filtervcfcov} {input.vcf} --coveragesummary {input.table} --min-coverage 25 --ignoremedian 50 > {output.vcf}"
		# Minimum coverage of 25x
		# As the maximum is not specified, the 95 quantile coverage of each sample in the table will be used instead'
		# Samples with a median coverage less than 50x will be ignored during filtering, making the filtering less harsh. (not used in the end)


# ------- Run GenMap -------

rule genmapindex:
	""" Make index for Genmap """ # super fast
	input:
		REFGenome
	output:
		"genmap/index/index.ids.concat"
	conda: 
		"envs/genmap.yaml"
	shell:
		"rm -r genmap/index; " # Snakemake makes the folder index, but then genmap tries to make it again and it fails
		"genmap index -F {input} -I genmap/index"

rule genmap:
	""" Run GenMap """ # super fast
	input:
		"genmap/index/index.ids.concat",
	output:
		"genmap/S288C.bedgraph"
	conda: 
		"envs/genmap.yaml"
	# params:
	# 	TMPDIR = "temp"
	shell:
		"mkdir -p temp; "
		"export TMPDIR=$PWD/temp; "
		"genmap map -K 100 -E 1 -I genmap/index -O genmap/S288C -t -bg"

rule filtergenmap:
	""" Get only good intevals of the genome """
	input:
		"genmap/S288C.bedgraph"
	output:
		"genmap/S288C_map1.bed"
	shell:
		"awk '$4==1' {input} | tr ' ' '\\t' > {output}"

rule bedtools:
	""" Get SNPs with good mappability """
	input:
		bed = "genmap/S288C_map1.bed",
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p.vcf"
	output:
		"vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_map1.vcf"
	shell:
		"bcftools view -h {input.vcf} > {output}; " #Get the header
		"bedtools intersect -a {input.vcf} -b {input.bed} >> {output}"

# ------- Plot the frequency of one of the parentales -------

rule plotY55freq:
	""" Plot the allele frequency of the Y55 allele along chromosomes for founder and G700 """
	input:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_map1.vcf",
	output:
		tablesnp = "data/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_map1_Y55freq_var.tab", # Used to find the common set in the vcf4adaptation.smk pipeline
		tablesnpall = "data/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_map1_Y55freq_all_var.tab", # Including sites that are monomorphic so freqY55 = 1 and de novo mutations used in vcf4adaptation_env.smk

		# Extra for exploring data, but didn't make it to the paper
		plot = "figures/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_map1_Y55freqG700.png",
		tablewin = "data/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_map1_Y55freq_win.tab",

	# conda: 
	# 	"envs/plot.yaml"
	params:
		time = "30:00",
		threads = 1,
		gen = "G700"
	script:
		Y55_along_contigs

# ------- Annotate with snpEff -------

rule prepare_snpEff:
	""" In order to use snpEff in conda I need to install the database """
	output:
		"logs/snpeff_installation.log"
	conda: 
		"envs/snpEff.yaml"
	shell:
		"snpEff download R64-1-1.99 && touch {output}"	

rule changecontignames:
	""" The contig names used in the reference genome are the Genbank accession numbers, but the snpEff database has the Roman letters instead for the chromosome names """
	input:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p.vcf",
	output:
		vcf = temp("vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_roman.vcf"),
	params:
		# threads = 1,
		time = "10:00",
	shell:
		"python {gb2romanS288C} {input.vcf} > {output.vcf}"

rule snpEff:
	""" Annotate with snpEff """
	input:
		log = "logs/snpeff_installation.log",
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_roman.vcf",
	output:
		vcf = temp("vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_roman_snpEff.vcf"),
		csv = "snpEff/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_roman_snpEff.csv",
	params:
		time = "3:00:00",
		threads = 2,
		ram = int(2 * 6.8) # A Rackham node contains 128 GB of RAM and 20 compute cores (each core gets at most 6.8 GB).
	conda: 
		"envs/snpEff.yaml"
	shell:
		"snpEff R64-1-1.99 -csvStats {output.csv} -nodownload -canon -v {input.vcf} > {output.vcf}"
		# "java -jar /sw/bioinfo/snpEff/4.3t/rackham/snpEff.jar S288C -csvStats All_Variants_PASS.csv -nodownload -canon -v All_Variants_PASS.vcf -c ~/Software/snpEff/snpEff.config > All_Variants_PASS_snpEff.vcf  "

rule changecontigback:
	""" The contig names used in the reference genome are the Genbank accession numbers, but the snpEff database has the Roman letters instead for the chromosome names """
	input:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_roman_snpEff.vcf",
	output:
		vcf = "vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_snpEff.vcf",
	params:
		# threads = 1,
		time = "10:00",
	shell:
		"python {gb2romanS288C} {input.vcf} --rom2gb > {output.vcf}"

rule getSnpEffTable:
	""" Make a table of the variant major effect from snpEff """
	input:
		"vcfs/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_snpEff.vcf",
	output:
		"snpEff/{outputname}_{env}_var_bi_miss0_{typevar}_25x95p_snpEff_majoreffect.tab",
	# params:
	# 	# threads = 1,
	# 	time = "10:00",
	shell:
		"python {snpEff2table} {input} > {output}"

# ------- Plot the de novo mutations -------

rule plotDeNovoMutations: 
	""" Plot the de novo mutations dynamics """
	input:
		EtOH_snp = f"data/{outputname}_Ethanol_var_bi_miss0_SNPs_25x95p_map1_Y55freq_all_var.tab",
		NaCl_snp = f"data/{outputname}_NaCl_var_bi_miss0_SNPs_25x95p_map1_Y55freq_all_var.tab",
		LiAc01_snp = f"data/{outputname}_LiAc0.01_var_bi_miss0_SNPs_25x95p_map1_Y55freq_all_var.tab",
		LiAc02_snp = f"data/{outputname}_LiAc0.02_var_bi_miss0_SNPs_25x95p_map1_Y55freq_all_var.tab",

		EtOH_indel = f"data/{outputname}_Ethanol_var_bi_miss0_INDELs_25x95p_map1_Y55freq_all_var.tab",
		NaCl_indel = f"data/{outputname}_NaCl_var_bi_miss0_INDELs_25x95p_map1_Y55freq_all_var.tab",
		LiAc01_indel = f"data/{outputname}_LiAc0.01_var_bi_miss0_INDELs_25x95p_map1_Y55freq_all_var.tab",
		LiAc02_indel = f"data/{outputname}_LiAc0.02_var_bi_miss0_INDELs_25x95p_map1_Y55freq_all_var.tab",	

		EtOH_snp_snpeff = f"snpEff/{outputname}_Ethanol_var_bi_miss0_SNPs_25x95p_snpEff_majoreffect.tab",	
		NaCl_snp_snpeff = f"snpEff/{outputname}_NaCl_var_bi_miss0_SNPs_25x95p_snpEff_majoreffect.tab",	
		LiAc01_snp_snpeff = f"snpEff/{outputname}_LiAc0.01_var_bi_miss0_SNPs_25x95p_snpEff_majoreffect.tab",	
		LiAc02_snp_snpeff = f"snpEff/{outputname}_LiAc0.02_var_bi_miss0_SNPs_25x95p_snpEff_majoreffect.tab",

		EtOH_indel_snpeff = f"snpEff/{outputname}_Ethanol_var_bi_miss0_INDELs_25x95p_snpEff_majoreffect.tab",	
		NaCl_indel_snpeff = f"snpEff/{outputname}_NaCl_var_bi_miss0_INDELs_25x95p_snpEff_majoreffect.tab",	
		LiAc01_indel_snpeff = f"snpEff/{outputname}_LiAc0.01_var_bi_miss0_INDELs_25x95p_snpEff_majoreffect.tab",	
		LiAc02_indel_snpeff = f"snpEff/{outputname}_LiAc0.02_var_bi_miss0_INDELs_25x95p_snpEff_majoreffect.tab",

		bad_mutations = bad_mutations
	output:
		infogenes = "results/deNovoGenesINFO_fixers_0.35.tab", # Raw mutations that reach more than 0.35 frequency at some point and the genes affected
		deNovofixers_trajectories = "results/deNovoFixers_0.35_trajectories_raw.tab", # Trajectories of raw mutations that reach more than 0.35 frequency at some point
		deNovofixers_trajectories_curated = "results/deNovoFixers_0.35_trajectories_curated.tab", # Trajectories of curated mutations that reach more than 0.35 frequency at some point (and hence considered adaptive)
		plotgenes = "results/deNovo_all_0.1_CoolGenes.pdf",
		genes4stringdb = "results/deNovoFixers_0.35_genes_4stringdb.txt",
		denovoeffectsp = "results/deNovoFixers_0.35_SNPeff.pdf",
	# conda: 
	# 	"envs/plot.yaml"
	params:
		time = "30:00",
		threads = 1,
	script:
		DeNovoMutations 

