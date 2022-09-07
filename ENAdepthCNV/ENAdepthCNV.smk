### ENAdepthCNV: Inferring copy number changes in the ENA gene from PoolSeq
#############################################################################

# In this pipeline we are interested in inferring changes in the copy number
# of a particular gene: ENA. This gene exist as a single copy in many *S.
# cerevisae* strains. In others, there are 3 other homologous copies of ENA
# genes at that locus, named ENA1, ENA2, and ENA5. From these, ENA1 is
# important and characterized as "P-type ATPase sodium pump; involved in Na+
# and Li+ efflux to allow salt tolerance". The trio (ENA1, ENA2, and ENA5)
# seems to have been introgressed from *S. paradoxus* and is present in the
# European lineage of *S. cerevisae*, including the reference strain S288c.

# Since we used SK1 and Y55 in our experiment, and these strains only have one
# ENA gene, we need to re-map all reads to one of their genomes to try to
# infer if there were changes in the copy number of these gene during
# experimental evolution.

#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022-08-30 to 2022-09-06
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

from glob import glob
import statistics # for mean

# -------------------------------------------------
# External files
REFGenome = "path/to/SK1.genome.fa"

# BED file with the regions of interest
bedGENE = "data/ENA_SK1.bed"

# List of sample's IDs
filewithIDs = "data/SampleIDs_AdaptationDynamics.txt"

# Scripts
ENAdepthCNV = "scripts/ENAdepthCNV.R"

# -------------------------------------------------
# Global variables
ENVIRONMENTS = ["NaCl", "Ethanol", "LiAc0.01", "LiAc0.02"]

# Samples to be removed because of bad coverage or contamination :(
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
"LiAc0.02_G1000_R5",
"LiAc0.01_G500_R5"] # Sequencing failed
# -------------------------------------------------

# Make a dictionary with the file containing the samples' IDs
# and remove all the samples that are contaminated
tabs = [line.rstrip("\n").split("\t")[0] for line in open(filewithIDs, 'r') if line.rstrip("\n").split("\t")[0] not in BADLIST] 			# Read tab file into a list

# ----------
# Rules not submitted to a job
localrules: rename_files, linkref, indexbwa, genmapindex, genmap, filtergenmap, bedtools, depthratio, plotreldepth
# ----------

rule all:
	input:
		"results/ENA_depth.png"

rule linkref:
	""" Make a link to the reference so it doesn't put new stuff in the ref folder """
	input:
		REFGenome
	output:
		"ref/genome.fa"
	shell:
		"ln -s {input} {output}"

rule indexbwa:
	""" Index genome with BWA """
	input:
		genome = "ref/genome.fa"
	output:
		index = "ref/genome.fa.bwt"
	shell:
		"""
		bwa index {input.genome}
		"""

rule bwa_mem:
	""" Map Illumina reads with BWA """
	input:
		genome = "ref/genome.fa",
		index = "ref/genome.fa.bwt",
		read1 = "raw-data/{sample}_1.fastq.gz",
		read2 = "raw-data/{sample}_2.fastq.gz",
	output:
		bwaoutput = temp("mapping/{sample}/{sample}.bam.sorted"),
	log:
		"logs/bwa_mem/{sample}.log"
	params:
		time = "3:30:00",
		threads = 10,
		rg = lambda wildcards: "@RG\\tID:" + wildcards.sample + "\\tSM:" + wildcards.sample + "\\tPL:ILLUMINA", # https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
	version: "1"
	shell:
		"""
		(bwa mem {input.genome} {input.read1} {input.read2} -t {params.threads} -R '{params.rg}' -M | samtools view -Su - | samtools sort -l 5 -O bam -T {wildcards.sample} -@ {params.threads} > {output.bwaoutput}) 2> {log}
		# -l 5 following Doug
		"""
rule markduplicates:
	""" Mark duplicates in BAM """
	# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
	input:
		bwaoutput = "mapping/{sample}/{sample}.bam.sorted"
	output:
		mdoutput = "mapping/{sample}/{sample}.sorted.debup.bam",
		mdmetrics = "mapping/{sample}/{sample}.sorted.metrics.txt"
	params:
		time = "1:30:00",
		threads = 5, # some jobs needs a lot of memory
	version: "1"
	shell:
		"""
		# Using normal Picard
		java -jar $PICARD_ROOT/picard.jar MarkDuplicates -I {input.bwaoutput} -O {output.mdoutput} -M {output.mdmetrics} --ASSUME_SORT_ORDER coordinate --CREATE_INDEX true --REMOVE_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
		"""	
		# # VALIDATION_STRINGENCY=ValidationStringency
		# #                               Validation stringency for all SAM files read by this program.  Setting stringency to
		# #                               SILENT can improve performance when processing a BAM file in which variable-length data
		# #                               (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This
		# #                               option can be set to 'null' to clear the default value. Possible values: STRICT,
		# #                               LENIENT, SILENT
		# # CREATE_INDEX=Boolean          Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value:
		# #                               false. This option can be set to 'null' to clear the default value. Possible values:
		# #                               true, false
		# # TMP_DIR (File)  Default value: null. This option may be specified 0 or more times.
		# # --REMOVE_DUPLICATES:Boolean   If true do not write duplicates to the output file instead of writing them with
        # # 		                      appropriate flags set.  Default value: false. Possible values: {true, false}
        # # --OPTICAL_DUPLICATE_PIXEL_DISTANCE:Integer
        # #                               The maximum offset between two duplicate clusters in order to consider them optical
        # #                               duplicates. The default is appropriate for unpatterned versions of the Illumina platform.
        # #                               For the patterned flowcell models, 2500 is moreappropriate. For other platforms and
        # #                               models, users should experiment to find what works best.  Default value: 100.
        # #	--ASSUME_SORT_ORDER,-ASO:SortOrder
        # #                               If not null, assume that the input file has this order even if the header says otherwise.
        # #                               Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate,
        # #                               unknown}  Cannot be used in conjunction with argument(s) ASSUME_SORTED (AS)

# ------- Run GenMap -------

rule genmapindex:
	""" Make index for Genmap """ # super fast
	input:
		"ref/genome.fa"
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
		"genmap/refgenome.bedgraph"
	conda: 
		"envs/genmap.yaml"
	shell:
		"mkdir -p temp; "
		"export TMPDIR=$PWD/temp; "
		"genmap map -K 100 -E 1 -I genmap/index -O genmap/refgenome -t -bg"

rule filtergenmap:
	""" Get only good intevals of the genome """
	input:
		"genmap/refgenome.bedgraph"
	output:
		"genmap/refgenome_map1.bed"
	shell:
		"awk '$4==1' {input} | tr ' ' '\\t' > {output}"

rule bedtools:
	""" Get SNPs with good mappability """
	input:
		bedgenmap = "genmap/refgenome_map1.bed",
		bedgene = bedGENE,
	output:
		"data/bedgene_map1.bed"
	shell:
		"bedtools intersect -a {input.bedgene} -b {input.bedgenmap} >> {output}"

# ------- Run Samtools -------

rule samtools:
	input:
		bed = "data/bedgene_map1.bed",
		bam = "mapping/{sample}/{sample}.sorted.debup.bam",
	output:
		"coverage/{sample}_ENA.bed"
	params:
		time = "30:00",
		threads = 1,
	shell:
		"samtools depth -a {input.bam} -b {input.bed} -Q 20 -o {output}"
		# -a	Output all positions (including those with zero depth)
		# -b FILE	Compute depth at list of positions or regions in specified BED FILE.
		# -H	Write a comment line showing column names at the beginning of the output. 
		# -o FILE	Write output to FILE.
		# -Q, --min-MQ INT	Only count reads with mapping quality greater than or equal to INT

# ------- Process the depth of coverage around the ENA gene and its flanks -------

rule depthratio:
	input:
		expand("coverage/{sample}_ENA.bed", sample = tabs),
	output:
		"results/ENA_depth.txt"
	run: # The coordinates of the ENA gene is chrIV:532984-536259
		outopen = open(output[0], 'w')
		outopen.write(f"Sample\tcovflanks\tNflank\tcovENA\tNena\tratio\n")
		enastart = 532984
		enaend = 536259
		flanklen = 100000

		for file in input: 
			sample = file.replace("coverage/", "").replace("_ENA.bed", "") # get the name from the file name
			print(sample) # some reporting to see progress

			flanks = []
			ena = []
			for line in open(file, 'r'):
				tabs = line.rstrip("\n").split("\t")
				contig, pos, cov = tabs
				pos = int(pos)

				if (pos < enastart) or (pos > enaend):
					flanks.append(int(cov))
				else:
					ena.append(int(cov))
			# Record output in file	
			outopen.write(f"{sample}\t{statistics.mean(flanks)}\t{len(flanks)}\t{statistics.mean(ena)}\t{len(ena)}\t{statistics.mean(ena)/statistics.mean(flanks)}\n")

rule plotreldepth:
	""" Plot a histogram of the relative depth of ENA compared to its flanks """
	input:
		table = "results/ENA_depth.txt"
	output:
		plot = "results/ENA_depth.png"
	script:
		ENAdepthCNV

