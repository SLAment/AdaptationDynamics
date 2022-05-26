### vcf4adaptation: Preparing the vcf file for the analysis of allele frequencies
#############################################################################

# The initial steps of variant calling were done previously by Nima Rafati. The 
# rest was done with the pipeline `varcall_adaptation.smk`.
# He originally modified the vcf file to split multiallelic sites in such a way that
# they are represented in more than one line. For example:

# $ grep -P 'BK006934.2\t1495\t' Raw_Call-SNP.vcf | cut -f1-10
# BK006934.2	1495	.	G	A	32531.97	.	AC=166;AF=0.241;AN=690	GT:AD:DP:GQ:PL	0/1:92,12:109:99:136,0,3034
# BK006934.2	1495	.	G	C	32531.97	.	AC=12;AF=0.017;AN=690	GT:AD:DP:GQ:PL	0/0:92,5:109:99:0,150,3197

# In the original file, the site looks like so:

# $ grep -P 'BK006934.2\t1495\t' Raw_Call.vcf | cut -f1-10
# BK006934.2	1495	.	G	A,C	32531.97	.	AC=166,12;AF=0.241,0.017;AN=690;BaseQRankSum=0.00;DP=38775;ExcessHet=139.8646;FS=37.983;InbreedingCoeff=-0.3996;MLEAC=183,10;MLEAF=0.265,0.014;MQ=33.44;MQRankSum=0.00;QD=1.89;ReadPosRankSum=0.00;SOR=4.220	GT:AD:DP:GQ:PL	0/1:92,12,5:109:99:136,0,3034,286,2915,3333

# I kept this design but in the end I didn't use these multiallelic sites.


#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022-03-30 - 2022-05-26
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 2

import os

# -------------------------------------------------
# External files
BASEVCF = "../../Variant_calling/All_Variants_parents_PASS.vcf" # Result of that pipeline
REFGenome = "../PerEnvironment/data/genome.fa"
chroms = "data/S288C_contigs.txt"
denovomutations = "../PerEnvironment/results/deNovoFixers_0.35_trajectories_curated.tab"

# Scripts to massage the data
vcf4adaptation_vcfR_plotter = "scripts/vcf4adaptation_vcfR_plotter_all.R"
filtervcfcov = "../PerEnvironment/scripts/filtervcfcov.py"
getvariantspool = "../PerEnvironment/scripts/getvariantspool.py"
extractmultivars = "../PerEnvironment/scripts/extractmultivars.py"

# Plotting
Y55_along_contigs = "../PerEnvironment/scripts/Y55_along_contigs.R"
Y55_findcommonset = "scripts/Y55_findcommonset.R"
haploselect = "scripts/haplocluster.R"
haplotrajectories = "scripts/haplotrajectories.R"
Y55trajectorySNPs = "scripts/Y55trajectorySNPs_stats.R"
ParallelFixation = "scripts/ParallelFixation.R"
ParallelFixationPlot = "scripts/ParallelFixationPlot.R"

# Global variables
ENVIRONMENTS = ["NaCl", "Ethanol", "LiAc0.01", "LiAc0.02"]

# -------------------------------------------------
# Define the output name from the input vcf file
# outputname = os.path.basename(BASEVCF).rstrip('.vcf')
outputname = os.path.splitext(os.path.basename(BASEVCF))[0]

# ----------
# Rules not submitted to a job
localrules: splitvcf_by_allelic, genmapindex, genmap, filtergenmap, bedtools, definesamples, listSamplesByEnv
# ----------

rule all:
	input:
		# Mention these two, so they both go through the same rule
		expand("figures/{outputname}_bi_miss0_{typevar}_Coverage.pdf", outputname = outputname, typevar = ["SNPs", "INDELs"]),
		expand("figures/{outputname}_bi_miss0_{typevar}_25x95p_map1_Coverage.pdf", outputname = outputname, typevar = ["SNPs", "INDELs"]),
		
		# -- Figures
		# Results figures: SGV raw data trajectories
		expand("results/Y55trajectory_{typevar}_novo_NaCl.pdf", typevar = "SNPs"),
		expand("results/Y55trajectory_{typevar}_novo_Ethanol.pdf", typevar = "SNPs"),
		expand("results/Y55trajectory_{typevar}_novo_LiAc0.01.pdf", typevar = "SNPs"),
		expand("results/Y55trajectory_{typevar}_novo_LiAc0.02.pdf", typevar = "SNPs"),
		# Results figures: loss of heterozygosity
		expand("results/SigmoidalPropFix_{typevar}.pdf", typevar = "SNPs"),
		expand("results/CoverageVsPropFixation_delta_{typevar}.pdf", typevar = "SNPs"),
		# Results figures: SGV in clusters
		f"results/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_shared_clusters_trajectories.pdf",
		f"results/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_shared_clusters_bigbois.pdf", # numbers
		# Results figures: parallel fixation of parental alleles
		"results/Parallel_heatmap_G100vsG700.png",
		"results/Parallel_count_G100vsG700.png",
		"results/Parallel_heatmap_G30.png",


rule splitvcf_by_allelic:
	""" Split the biallelic sites from multiallelic sites in different files """
	input:
		BASEVCF
		# "vcfs/{outputname}_goodsamples.vcf"
	output:
		bi = "vcfs/{outputname}_bi.vcf",
		multi = "vcfs/{outputname}_multi.vcf", # We didn't use them in the end
		log = "logs/{outputname}_extractmultivars.txt"
	shell:
		"python {extractmultivars} {input} --outputdir 'vcfs' --outname {wildcards.outputname} > {output.log}"

rule compressvcfs:
	input:
		"vcfs/{outputname}_bi.vcf"
	output:
		"vcfs/{outputname}_bi.vcf.gz"
	wildcard_constraints:
		outputname = outputname,
	params:
		time = "20:00",
		# threads = 1,
	shell:
		"bgzip {input}"

rule filter4missingdata:
	""" Remove sites with missing data """
	input:
		"vcfs/{outputname}_bi.vcf.gz"
	output:
		"vcfs/{outputname}_bi_miss0.vcf.gz"
	params:
		threads = 1,
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
		"vcfs/{outputname}_bi_miss0.vcf.gz"
	output:
		"vcfs/{outputname}_bi_miss0.vcf.gz.tbi"
	params:
		threads = 1,
		time = "30:00",
	shell:
		"tabix -p vcf {input}"

rule getSNPs:
	""" Make separate vcf file for just the SNPs """
	input:
		vcf = "vcfs/{outputname}_bi_miss0.vcf.gz",
		index = "vcfs/{outputname}_bi_miss0.vcf.gz.tbi"
	output:
		"vcfs/{outputname}_bi_miss0_SNPs.vcf"
	params:
		threads = 1,
		time = "1:00:00",		
	shell:
		"gatk SelectVariants -R {REFGenome} -V {input.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {output}"
		# This however will still keep sites where the alternative is "*", which means a deletion of that base.

rule getINDELs:
	""" Make separate vcf file for just the INDELs """
	input:
		vcf = "vcfs/{outputname}_bi_miss0.vcf.gz",
		index = "vcfs/{outputname}_bi_miss0.vcf.gz.tbi"
	output:
		"vcfs/{outputname}_bi_miss0_INDELs.vcf"
	params:
		threads = 1,
		time = "1:00:00",		
	shell:
		"gatk SelectVariants -R {REFGenome} -V {input.vcf} --select-type-to-include INDEL --restrict-alleles-to BIALLELIC -O {output}"

# ------- Find sites that are not good -------

rule PlotCoverage:
	""" Plot coverage distribution for all samples and report """
	input:
		vcf = "vcfs/{outputname}_bi_miss0_{typevar}.vcf"
	output:
		plot = "figures/{outputname}_bi_miss0_{typevar}_Coverage.pdf", # All four environments at the same time
		table = "data/{outputname}_bi_miss0_{typevar}_CoverageSummary.tab"
	params:
		time = "10:00:00",
		threads = 10, # It needs some memory
		lquantile = 0.25,
		uquantile = 0.95,
	script:
		vcf4adaptation_vcfR_plotter

# ------- Filter -------

rule filtervcfcov:
	""" Filter vcf file based on the coverage distribution per sample """
	input:
		vcf = "vcfs/{outputname}_bi_miss0_{typevar}.vcf",
		table = "data/{outputname}_bi_miss0_{typevar}_CoverageSummary.tab"
	output:
		vcf = "vcfs/{outputname}_bi_miss0_{typevar}_25x95p.vcf",
	params:
		time = "3:00:00",
		threads = 1,
	shell:
		"python {filtervcfcov} {input.vcf} --coveragesummary {input.table} --min-coverage 25 > {output.vcf}"
		# Minimum coverage of 25x
		# As the maximum is not specified, the 95 quantile coverage of each sample in the table will be used instead'

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
		vcf = "vcfs/{outputname}_bi_miss0_{typevar}_25x95p.vcf",
	output:
		"vcfs/{outputname}_bi_miss0_{typevar}_25x95p_map1.vcf"
	shell:
		"bcftools view -h {input.vcf} > {output}; " #Get the header
		"bedtools intersect -a {input.vcf} -b {input.bed} >> {output}"

# ------- Divide by environments -------

rule definesamples:
	""" Get a list of samples to extract from the original vcf file """
	input:
		BASEVCF
	output:
		"data/samples.txt"
	shell:
		"bcftools view -h {input} | grep '#CHROM' | cut -f10- | tr '\\t' '\\n' > {output}"

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
		vcf = "vcfs/{outputname}_bi_miss0_{typevar}_25x95p_map1.vcf",
		lista = "data/samples_{env}.txt"
	output:
		"vcfs/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}.vcf.gz"
	wildcard_constraints:
		outputname = outputname,
		env = '|'.join([x for x in ENVIRONMENTS])
	params:
		# threads = 1,
		time = "1:00:00",
	shell:
		"vcftools --vcf {input.vcf} --recode --recode-INFO-all --stdout --keep {input.lista} | bgzip > {output}; "
		"tabix -p vcf {output}"

# ------- Plot the frequency of one of the parentales -------

rule plotY55freq:
	""" Plot the allele frequency of the Y55 allele along chromosomes for founder and G700 """
	input:
		vcf = "vcfs/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}.vcf.gz",
	output:
		tablesnp = "data/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}_Y55freq_var.tab", # Used to find the common set in the vcf4adaptation.smk pipeline
		tablesnpall = "data/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}_Y55freq_all_var.tab", # Including sites that are monomorphic so freqY55 = 1 and de novo mutations used in vcf4adaptation_env.smk

		# Extra for exploring data, but didn't make it to the paper
		plot = "figures/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}_Y55freqG700.png",
		tablewin = "data/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}_Y55freq_win.tab",
	# conda: 
	# 	"envs/plot.yaml"
	params:
		time = "30:00",
		threads = 1,
		gen = "G700"
	script:
		Y55_along_contigs

rule Y55_findcommonset:
	"""  Subset the SNPs of all environments to keep only those present in all samples """
	input:
		expand("data/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}_Y55freq_var.tab", env = ENVIRONMENTS, outputname = outputname, typevar = "SNPs"),
		expand("vcfs/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}.vcf.gz", env = ENVIRONMENTS, outputname = outputname, typevar = "SNPs"),
	output:
		expand("data/{outputname}_bi_miss0_{typevar}_25x95p_map1_allenv_Y55freq_var_shared.tab", outputname = outputname, typevar = "SNPs"),
		expand("vcfs/{outputname}_bi_miss0_{typevar}_25x95p_map1_{env}_shared.vcf.gz", env = ENVIRONMENTS, outputname = outputname, typevar = "SNPs"),
	params:
		time = "30:00",
		threads = 1,
	script:
		Y55_findcommonset

rule haploselect:
	""" Find alleles under selection and define haplotypes """ 
	input:
		chroms = chroms,
		nacl = f"vcfs/{outputname}_bi_miss0_SNPs_25x95p_map1_NaCl_shared.vcf.gz",
		ethanol = f"vcfs/{outputname}_bi_miss0_SNPs_25x95p_map1_Ethanol_shared.vcf.gz",
		L1 = f"vcfs/{outputname}_bi_miss0_SNPs_25x95p_map1_LiAc0.01_shared.vcf.gz",
		L2 = f"vcfs/{outputname}_bi_miss0_SNPs_25x95p_map1_LiAc0.02_shared.vcf.gz",
	output:
		table = f"data/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_shared_clusters.txt"
	params:
		time = "30:00",
		threads = 1,
	script:
		haploselect

rule plotSGVandmutations:
	""" Put together the SGV haplotypes and the de novo mutations """
	input:
		clusters = f"data/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_shared_clusters.txt",
		rawsnps = f"data/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_Y55freq_var_shared.tab",
		mutations = denovomutations,
	output:
		sgvplot = f"results/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_shared_clusters_trajectories.pdf", # Figure in the paper
		bigbois = f"results/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_shared_clusters_bigbois.pdf",
	params:
		time = "15:00",
		threads = 1,
	script:
		haplotrajectories		

rule Y55trajectorySNPs:
	""" Plot the allele frequency of the SK1 allele along chromosomes for all samples """
	input:
		table = expand("data/{outputname}_bi_miss0_{typevar}_25x95p_map1_allenv_Y55freq_var_shared.tab", outputname = outputname, typevar = "SNPs"),
		mutations = denovomutations,
		coverage = expand("data/{outputname}_bi_miss0_{typevar}_25x95p_map1_CoverageSummary.tab", outputname = outputname, typevar = "SNPs"),
	output:
		nacl = "results/Y55trajectory_{typevar}_novo_NaCl.pdf",
		ethanol = "results/Y55trajectory_{typevar}_novo_Ethanol.pdf",
		li1 = "results/Y55trajectory_{typevar}_novo_LiAc0.01.pdf",
		li2 = "results/Y55trajectory_{typevar}_novo_LiAc0.02.pdf",
		sigmoid = "results/SigmoidalPropFix_{typevar}.pdf", # Figure in the paper
		covcorr = "results/CoverageVsPropFixation_delta_{typevar}.pdf", # Supplementary figure in the paper
	params:
		time = "15:00",
		threads = 1,
	script:
		Y55trajectorySNPs

rule parallelfixationcalc:
	""" Find regions in the genome that fixed for the same allele """
	input:
		table = f"data/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_Y55freq_var_shared.tab",
		mutations = denovomutations,
	output:
		wins = f"data/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_Y55freq_var_shared_wins10kb_min4.tab",
	script:
		ParallelFixation

rule parallelfixationplot:
	""" Find regions in the genome that fixed for the same allele """
	input:
		wins = f"data/{outputname}_bi_miss0_SNPs_25x95p_map1_allenv_Y55freq_var_shared_wins10kb_min4.tab",
	output:
		G100vsG700 = "results/Parallel_heatmap_G100vsG700.png",
		fixwins = "results/Parallel_count_G100vsG700.png",
		G30 = "results/Parallel_heatmap_G30.png",		
	script:
		ParallelFixationPlot
