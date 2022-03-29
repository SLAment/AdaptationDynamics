# Adaptation Dynamics

Code associated to the project ["The dynamics of adaptation to stress from standing genetic variation"](https://www.biorxiv.org/content/10.1101/2022.03.26.485920v1) from the Stelkens Lab.

There are three main [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines the deal with the variant data, analysis and plotting. 

1. `varcall_adaptation.smk` -- For variant calling (directory `Variant_calling`)
2. `vcf4adaptation_env.smk` -- For variant filtering and analysis of de novo mutations (directory `Variant_filtering/PerEnvironment/`)
3. `vcf4adaptation.smk` -- For variant filtering and analysis of the standing genetic variation (directory `Variant_filtering/AllTogether/`) -- TODO

The repository also contains code of Ciaran Gilchrist analysing the phenotypic data.

----

I ran the pipelines in Uppsala University's supercomputer [UPPMAX](https://uppmax.uu.se/), which has a CentOS Linux operating system with a slurm scheduler. However, they should work fine also in other unix environments.
