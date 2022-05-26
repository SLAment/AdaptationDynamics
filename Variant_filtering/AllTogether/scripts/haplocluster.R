#!/usr/bin/env Rscript

# haplocluster: Cluster similar trajectories of SNPs
#############################################################################
# Modified version of the script `haplo_full_fig.R` by Alexandre Rego
# =======================================
# Modified by S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022-04-20
# Version 1
#############################################################################
# Load the necessary libraries
# ============================

library(dplyr)
library(ggplot2)
library(vcfR)
library(gtools) # for mixedsort()
library(tidyr) # For separate()

# ============================
# Set paths
# ============================
vcffile_EtOH <- snakemake@input$ethanol
vcffile_NaCl <- snakemake@input$nacl 
vcffile_LiAc01 <- snakemake@input$L1
vcffile_LiAc02 <- snakemake@input$L2

# Output
clustersfile <- snakemake@output$table

# ============================
# Reading the data
# ============================
vcfE <- read.vcfR(vcffile_EtOH, verbose = FALSE)
vcfN <- read.vcfR(vcffile_NaCl, verbose = FALSE) 
vcfL1 <- read.vcfR(vcffile_LiAc01, verbose = FALSE) 
vcfL2 <- read.vcfR(vcffile_LiAc02, verbose = FALSE)

# ============================
# Define CONSTANTs
# ============================
# Minimum number of SNPs to consider a cluster a big boi
MIN_N_VAR <- 10

# ============================
# Define helper functions
# ============================
## Function to turn AD matrix into a DP matrix
getDPfromADMatrix <- function(ADdf){
  dp <- data.frame(matrix(sapply(strsplit(ADdf,','),function(x) sum(as.numeric(x))),
                          nrow=nrow(ADdf),ncol=ncol(ADdf)))
  colnames(dp) <- colnames(ADdf)
  rownames(dp) <- rownames(ADdf)
  return(dp)
}

# Normalize allele freqs, from the HaploValidate method (step 2)
transform.af <- function(af) {
  af.sqrt <- asin(sqrt(af))
  af.transf <- t(af.sqrt)
  af.scale <- scale(af.transf, center = TRUE, scale = TRUE)
  return(af.scale)
}

# ============================
cat("Process data\n")
# ============================

intN <- extract.gt(vcfN, element="AD")
intE <- extract.gt(vcfE, element="AD")
intL1 <- extract.gt(vcfL1, element="AD")
intL2 <- extract.gt(vcfL2, element="AD")

# Get allele frequencies of the reference allele
dpN <- getDPfromADMatrix(intN)
refafN <- extract.gt(vcfN, element="AD", as.numeric=T)/dpN

dpL1 <- getDPfromADMatrix(intL1)
refafL1 <- extract.gt(vcfL1, element="AD", as.numeric=T)/dpL1

dpE <- getDPfromADMatrix(intE)
refafE <- extract.gt(vcfE, element="AD", as.numeric=T)/dpE

dpL2 <- getDPfromADMatrix(intL2)
refafL2 <- extract.gt(vcfL2, element="AD", as.numeric=T)/dpL2

## Big function to cluster all variants given at any ime point and classify them
## into clusters based on their correlations, in such a way to maximize the number 
## of haplotypes with equal or more than `MIN_N_VAR` variants
clusterall <- function(founder, replicate, refaf, dp, condition){
  # founder = 'N'; replicate = 2; refaf = refafN; dp = dpN; condition = "NaCl_R2"
  cat(founder, replicate, "\n")
  # Get names of the relevant samples
  group <- colnames(refaf)[grep(paste0(founder,'_Founder|R',replicate), colnames(refaf))] ## Replace grep input with appropriate founder and replicate number

  # ## Names of all time points and replicates of interest sorted
  prettyorder <- gtools::mixedsort(group) 
  af_set <- refaf[,prettyorder] ## data.table with sites and samples

  ##### Making haplotype blocks
  af_trans <- transform.af(af_set) ## transform data
  tradj_cors <- cor(af_trans) ## get correlations
  cordist <-  1 - as.dist(abs(tradj_cors)) ## Convert to distance format
  corclust <- hclust(cordist, method='average') # perform clustering
  
  big_boi_len <- c() # Big bois being the big haplotypes
  corr_cut <- seq(0.9,0.3,-0.1) # Explore a range of minimum correlation coefficients (as in the HaploValidate paper Step 1)
  for(i in corr_cut){ #from min corr of 0.9 to 0.3
    clust_ids <- cutree(corclust, h = 1-i)
    clust_table <- table(clust_ids)
    big_bois <- which(clust_table >= MIN_N_VAR) # MIN_N_VAR being the minimum number of SNPs present in the haplotype to be a big boi
    big_boi_len <- c(big_boi_len, length(big_bois))
  }
  
  clust_ids <- cutree(corclust, h = 1 - corr_cut[which.max(big_boi_len)]) ### 1 - max distance based on which correlation produced the most number of haplotypes. If same # then choose less strict correlation.
  # clust_ids <- cutree(corclust, h = 1 - corr_cut[which.min(big_boi_len)]) ### 1 - min distance based on which correlation produced the least number of haplotypes. If same # then choose less strict correlation.
  
  # Now we chose the correlation coefficient, get the big bois for the chosen one
  clust_table <- table(clust_ids)
  
  ## ---- Different from Alex' ----
  ### In order to make it work with my other scripts and data frames, I need the 
  # cluster ID and the number of variants in each cluster, so later I can tell 
  # which ones are big bois
  af_fin_al <- cbind(af_set, pos = row.names(af_set), clust = clust_ids) # This are all SNPs, not just big bois
  clustvariants <- data.frame(Condition = condition, Coordinate = rownames(af_fin_al), clust = af_fin_al$clust) %>% 
    mutate(N_var = clust_table[as.numeric(clust)])
  
  return(clustvariants)
}

cat("Get selected sites and cluster them...\n")
## Create a huge data frame with the selected clusters of each replicate of each treatment
# Takes some time, about a minute and a half
startTime <- Sys.time()
uberdf <- rbind(clusterall('N', 1, refafN, dpN, "NaCl_R1"),
                clusterall('N', 2, refafN, dpN, "NaCl_R2"),
                clusterall('N', 3, refafN, dpN, "NaCl_R3"),
                clusterall('N', 4, refafN, dpN, "NaCl_R4"),
                clusterall('LE', 1, refafE, dpE, "Ethanol_R1"),
                clusterall('LE', 2, refafE, dpE, "Ethanol_R2"),
                clusterall('LE', 3, refafE, dpE, "Ethanol_R3"),
                clusterall('LE', 4, refafE, dpE, "Ethanol_R4"),
                clusterall('LE', 1, refafL1, dpL1, "LiAc0.01_R1"),
                clusterall('LE', 2, refafL1, dpL1, "LiAc0.01_R2"),
                clusterall('LE', 4, refafL1, dpL1, "LiAc0.01_R4"),
                clusterall('LE', 5, refafL1, dpL1, "LiAc0.01_R5"),
                clusterall('LE', 1, refafL2, dpL2, "LiAc0.02_R1"),
                clusterall('LE', 2, refafL2, dpL2, "LiAc0.02_R2"),
                clusterall('LE', 3, refafL2, dpL2, "LiAc0.02_R3"),
                clusterall('LE', 4, refafL2, dpL2, "LiAc0.02_R4"),
                clusterall('LE', 5, refafL2, dpL2, "LiAc0.02_R5"))
endTime <- Sys.time()
endTime-startTime

## Save dataframe for plotting later
write.table(uberdf, file = clustersfile, sep = "\t", row.names=FALSE, quote=FALSE)

