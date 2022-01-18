# README: Call variants for the adaptation dynamics paper

For this project, pool sequencing data of diferente populations of *S. cerevisiae* evolving on four different environments were produced. The gvcf files were produced independently from this pipeline. A few samples were discarded based on contamination with *S. paradoxus*, low coverage or potential contamination from within the experiment. The excluded samples are:

	LiAc0.01_G1000_R1
	LiAc0.01_G1000_R2
	LiAc0.01_G1000_R4
	LiAc0.01_G1000_R5
	LiAc0.02_G60_R2
	NaCl_G1000_R4
	NaCl_G1000_R5
	NaCl_G100_R5
	NaCl_G200_R5
	NaCl_G300_R5
	NaCl_G30_R5
	NaCl_G400_R5
	NaCl_G500_R5
	NaCl_G60_R5
	NaCl_G700_R5

This pipeline was run on the supercomputer [UPPMAX](https://uppmax.uu.se/), which has a CentOS Linux operating system with a slurm scheduler. However, they should work fine also in other unix environments.

## Data needed for the pipeline

- Folder with gvcf files created with GATK
- The reference genome
- A BED file of the genome to be used as INTERVALS by GATK


## Building the environment

The pipeline as it is depends on the following modules in UPPMAX:

	$ module load bioinfo-tools snakemake/5.30.1 GATK/4.1.4.1 samtools/1.9


## Pipeline

This is a pipeline build on the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow manager. The paths to the necessary files are on the top of the pipeline and should be changed by the user.

First, to get an idea of how the pipeline looks like we can make a rulegraph:

    $ snakemake --snakefile varcall_adaptation.smk --rulegraph | dot -Tpng > rulegraph.png

![rulegraph](rulegraph.png "rulegraph")

To check that the files for the pipeline are in order:

	$ snakemake --snakefile varcall_adaptation.smk -pn

So the pipeline assumes the BAM files are already in said folders.


	$ screen -R varcalling
	$ module load bioinfo-tools snakemake/5.30.1 GATK/4.1.4.1 samtools/1.9
	$ snakemake --snakefile varcall_adaptation.smk -p --cluster "sbatch -A snicXXXX-X-XXX -p core -n {params.threads} -t {params.time} --mail-user xxxxxx@xxxxx.xx --mail-type=ALL" -j 30 --keep-going &> snakemake.log &


## About UPPMAX

When running this pipeline, UPPMAX had the following specifications:

	$ hostnamectl
	   Static hostname: rackham2.uppmax.uu.se
	         Icon name: computer-server
	           Chassis: server
	        Machine ID: f911affe94fa4ccb8e6deebe489bdd9b
	           Boot ID: 20854dc584e94630ae6d6860d5e94dd7
	  Operating System: CentOS Linux 7 (Core)
	       CPE OS Name: cpe:/o:centos:centos:7
	            Kernel: Linux 3.10.0-1160.41.1.el7.x86_64
	      Architecture: x86-64

	$ java -version
	java version "1.8.0_151"
	Java(TM) SE Runtime Environment (build 1.8.0_151-b12)
	Java HotSpot(TM) 64-Bit Server VM (build 25.151-b12, mixed mode)
