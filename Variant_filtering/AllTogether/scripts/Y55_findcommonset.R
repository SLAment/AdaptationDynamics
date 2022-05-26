#!/usr/bin/env Rscript

# Y55_findcommonset: Subset the variants of all environments to keep only those present in all samples
#############################################################################
# Part of the Snakemake pipeline vcf4adaptation.smk 
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022-01-21
# Version 1
#############################################################################
# Load the necessary libraries
# ============================
library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
library(vcfR)

# ============================
# Set paths
# ============================
## Tables
tablemaf_NaCl <- snakemake@input[[1]]
tablemaf_EtOH <- snakemake@input[[2]]
tablemaf_LiAc01 <- snakemake@input[[3]]
tablemaf_LiAc02 <- snakemake@input[[4]]

## vcf files
vcffile_EtOH <- snakemake@input[[5]]
vcffile_NaCl <- snakemake@input[[6]]
vcffile_LiAc01 <- snakemake@input[[7]]
vcffile_LiAc02 <- snakemake@input[[8]]

## Outputs
trajectos <- snakemake@output[[1]]
vcfset_EtOH <- snakemake@output[[2]]
vcfset_NaCl <- snakemake@output[[3]]
vcfset_LiAc01 <- snakemake@output[[4]]
vcfset_LiAc02 <- snakemake@output[[5]]

# ============================
# Reading the data
# ============================
allfreq_NaCl <- read.table(tablemaf_NaCl, header=TRUE)
allfreq_EtOH <- read.table(tablemaf_EtOH, header=TRUE)
allfreq_LiAc01 <- read.table(tablemaf_LiAc01, header=TRUE)
allfreq_LiAc02 <- read.table(tablemaf_LiAc02, header=TRUE)

vcf_EtOH <- read.vcfR(vcffile_EtOH)
vcf_NaCl <- read.vcfR(vcffile_NaCl)
vcf_LiAc01 <- read.vcfR(vcffile_LiAc01)
vcf_LiAc02 <- read.vcfR(vcffile_LiAc02)

# ============================
# Modify the dataframes to assign the Founders
# ============================
# In order to have the 0 time point for every replicate I have to repeat the 
# data frame of the founder for all of them. (Salt has a different founder)

fixfounder <- function(allfreq){
  # Get the founder
  allfreq_env <- allfreq %>% filter(!sample %in% c("SK1", "Y55"))
  allfreq_env_noF <- allfreq_env %>% filter(!Generation %in% c("Founder"))
  allfreq_env_F <- allfreq_env %>% filter(Generation %in% c("Founder"))
  
  # For this environment
  allfreq_env_F$Environment <- unique(allfreq_env_noF$Environment)
  
  new_allfreq_env <- allfreq_env_noF
  for (r in allfreq_env_noF$Replicate %>% unique()) {
    allfreq_env_F$Replicate <- r
    allfreq_env_F$Replicate <- r
    new_allfreq_env <- rbind(new_allfreq_env, allfreq_env_F)
  }
  
  return(new_allfreq_env)
}

# put them together
allfreq <- rbind(fixfounder(allfreq_EtOH),
                 fixfounder(allfreq_NaCl),
                 fixfounder(allfreq_LiAc01),
                 fixfounder(allfreq_LiAc02)) %>% 
  filter(!(gen == 1000 & Environment %in% c("LiAc0.01", "LiAc0.02"))) %>% filter(!sample == "NaCl_G1000_R1") # Bad samples

# Change names of chromosomes
allfreq$Contig <- as.factor(allfreq$Contig)
allfreq$Contig <- recode_factor(allfreq$Contig, BK006935.2 = "ChrI", BK006936.2 = "ChrII", BK006937.2 = "ChrIII", BK006938.2 = "ChrIV", BK006939.2 = "ChrV", BK006940.2 = "ChrVI", BK006941.2 = "ChrVII", BK006934.2 = "ChrVIII", BK006942.2 = "ChrIX", BK006943.2 = "ChrX", BK006944.2 = "ChrXI", BK006945.2 = "ChrXII", BK006946.2 = "ChrXIII", BK006947.3 = "ChrXIV", BK006948.2 = "ChrXV", BK006949.2 = "ChrXVI")

# ============================
# Modify the dataframes to assign the Founders
# ============================
# They way I made these files, I only kept sites that were polymorphic in the 
# founders (maf>0.3), so I need to make a common set because there are two founders
LEfounder <- filter(allfreq, sample == "LE_Founder_R1", Environment == "LiAc0.01", Replicate == "R1") 
Nfounder <- filter(allfreq, sample == "N_Founder_R1", Replicate == "R1") 

# Find the common set
commonset <- intersect(LEfounder$Coordinate, Nfounder$Coordinate)
allfreq_set <- filter(allfreq, Coordinate %in% commonset) %>% select(-POS)
cat("No. of SNPs per sample:", length(commonset), "\n")
# ============================
## Save dataframe for plotting later
# ============================
write.table(allfreq_set, file = trajectos, sep = "\t", row.names=FALSE, quote=FALSE)

# ============================
# Subset
# ============================
# goodsites is a list with concatenated CHROM and POS positions (what I call Coordinate)
subsetvcf <- function(vcf, goodsites){
  chrsvcf <- vcf@fix[, "CHROM"] 
  posvcf <- vcf@fix[, "POS"] 
  coordvcf <- paste0(chrsvcf, "_", posvcf)
  
  # Subset vcf file
  subsetvcf <- vcf
  subsetvcf@fix <- vcf@fix[which(coordvcf %in% goodsites), ]
  subsetvcf@gt <- vcf@gt[which(coordvcf %in% goodsites), ]
  return(subsetvcf)
}


# Save to files
write.vcf(subsetvcf(vcf_EtOH, commonset), file = vcfset_EtOH)
write.vcf(subsetvcf(vcf_NaCl, commonset), file = vcfset_NaCl)
write.vcf(subsetvcf(vcf_LiAc01, commonset), file = vcfset_LiAc01)
write.vcf(subsetvcf(vcf_LiAc02, commonset), file = vcfset_LiAc02)

