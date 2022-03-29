#!/usr/bin/env Rscript

# vcf4adaptation_vcfR_plotter: Plot the coverage distribution from a vcf file
#############################################################################
# Part of the Snakemake pipeline vcf4adaptation.smk 
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021-09-24
# Version 1
#############################################################################

# ============================
# Check input
# ============================
vcffile <- snakemake@input$vcf
plotcov <- snakemake@output$plot
tablecov <- snakemake@output$table

lquantile = snakemake@params$lquantile
uquantile = snakemake@params$uquantile

# ============================
# Load the necessary libraries
# ============================
library(vcfR)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
library(tidyr) # For separate()
library(cowplot)
library(reshape2)

# ============================
# Reading the data
# ============================
vcf <- read.vcfR(vcffile, verbose = FALSE)

# ============================
# Prepare the functions
# ============================
# https://knausb.github.io/vcfR_documentation/sequence_coverage.html
coolviolins <- function(dp, color = "#C0C0C0"){
  p <- ggplot(dp, aes(x=Sample, y=Depth)) + geom_violin(fill=color, colour = color, adjust=1.0, scale = "count", trim=TRUE)
  p <- p + theme_bw()
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.text.x = element_text(angle = 60, hjust = 1))
  p <- p + stat_summary(fun.data=median_hilow, geom="pointrange", color="black") # by default the lower and upper quantiles computed are 0.025 and 0.975.
  p
  return(p)
}

# ============================
# Find multiallelic sites and recover only one representative
# ============================
# Make data frame of depth of coverage
# dp <- extract.gt(vcf, element="DP", as.numeric = TRUE)

# Because DP still includes bad quality reads in GATK, I should use the AD field
# Make it long so I can get the allele depth
dp_long <- extract.gt(vcf, element="AD") %>% data.frame() %>%
  tibble::rownames_to_column("POS") %>%
  gather("sample", "AD", -POS) %>% 
  separate(AD, c("ref", "alt"), sep = ",") 

dp_long$ref <- as.numeric(dp_long$ref)
dp_long$alt <- as.numeric(dp_long$alt)
# Finally get the total depth
dp_long <- mutate(dp_long, DP = ref + alt) %>% select(POS, sample, DP)

## Turn back into wide format
dp <- spread(dp_long, sample, DP )

## Find repeated sites
# How many sites are duplicated?
sites <- dp %>% separate(POS, c("contig", "POS", "n"), "_")
sites$Coordinate <- paste(sites$contig, sites$POS, sep = "_")

# Which sites are those?
# https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r
n_occur <- data.frame(table(sites$Coordinate))
# n_occur[n_occur$Freq > 1,] # Some sites are more than two times (3 and 4 times)
# Remove them
dp_nomulti <- dp[sites$Coordinate %in% n_occur$Var1[n_occur$Freq == 1],] 

# Now let's deal with the ones that are multiallelic
multisites <- sites[sites$Coordinate %in% n_occur$Var1[n_occur$Freq > 1],] # Get the multisites in data frame with site names
dp_multip <- dp[sites$Coordinate %in% n_occur$Var1[n_occur$Freq > 1],] # Get the multisites in original format
# What elements are repeated from those?
keep <- multisites$Coordinate %>% duplicated %>% rev
# Retrieve only one of them
dp_multiunique <- dp_multip[keep,]

#Put back together
dp <- rbind(dp_nomulti, dp_multiunique) # 41712 sites in total

# ============================
# Print violin plots of coverage
# ============================

# Make long format again
dpf <- gather(dp, "Sample", "Depth", -POS)
dpf <- separate(dpf, Sample, c("Environment", "Generation", "Replicate"), sep = "_", remove = FALSE) # Get details, but it takes a while


# Determine maximum coverage
rawmax <- dpf$Depth %>% max

# What is the color of this plot?
if ("NaCl" %in% unique(dpf$Environment)) {
  # thiscolor = "cadetblue3"
  thiscolor = "#4169E1"
} else if ("Ethanol" %in% unique(dpf$Environment)){
  # thiscolor = "coral4"
  thiscolor = "#AC1016"
} else if ("LiAc0.01" %in% unique(dpf$Environment)){
  # thiscolor = "darkgoldenrod3"
  thiscolor = "#800080"
} else if ("LiAc0.02" %in% unique(dpf$Environment)){
  # thiscolor = "darkseagreen4"
  thiscolor = "#228B22"
}

# Plot
envviolin <- coolviolins(dpf, color = thiscolor) + geom_hline(yintercept=50, linetype="dashed", color = "black")

ggsave(plotcov, plot = envviolin, width = 20, height = 10, units = "cm")

# Print a summary table of the coverage distribution per sample
dp_noPOS <- dp %>% select(-POS)
dpsum <- data.frame(sample = colnames(dp_noPOS), 
                    median = apply(dp_noPOS, 2, median, na.rm = TRUE),
                    mean = apply(dp_noPOS, 2, mean, na.rm = TRUE), 
                    lowquantile = apply(dp_noPOS, 2, quantile, probs = lquantile, na.rm = TRUE),
                    upquantile = apply(dp_noPOS, 2, quantile, probs = uquantile, na.rm = TRUE))

names(dpsum)[4] <- paste0("q",lquantile)
names(dpsum)[5] <- paste0("q",uquantile)

write.table(dpsum, file = tablecov, sep = "\t", row.names=FALSE, quote=FALSE)

sessionInfo()
