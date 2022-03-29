#!/usr/bin/env Rscript

# Y55_along_contigs: Plot the allele frequency of Y55 from the adaptation 
# experiment along the genome
#############################################################################
# Part of the Snakemake pipeline vcf4adaptation.smk 
# Here I spent a lot of time massaging the data frame to basically get the allele 
# frequency of the Y55 parental in all the samples
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021/10/15-19
# Version 2
#############################################################################
# Load the necessary libraries
# ============================
library(vcfR)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
library(tidyr) # For separate()
library(cowplot)
library(stringr) # For natural sorting using str_sort()

# ============================
# Check input
# ============================
vcffile <- snakemake@input$vcf
ploty55 <- snakemake@output$plot
tabley55win <- snakemake@output$tablewin
tabley55snp <- snakemake@output$tablesnp
tabley55_allsnp <- snakemake@output$tablesnpall
interestgen <- snakemake@params$gen

# ============================
# Reading the data
# ============================
cat("Reading vcf file ...\n")
vcf <- read.vcfR(vcffile, verbose = FALSE)

# Get genotypes
gt <- extract.gt(vcf, element='GT', return.alleles = TRUE) %>% data.frame()
gt <- tibble::rownames_to_column(gt, "POS") %>% data.frame() 

# ============================
# Get allele depth
# ============================
dp <- extract.gt(vcf, element="AD") %>% data.frame() %>%
  tibble::rownames_to_column("POS") %>% 
  gather("sample", "AD", -POS) %>% 
  separate(AD, c("refAD", "altAD"), sep = ",") 

dp$refAD <- as.numeric(dp$refAD)
dp$altAD <- as.numeric(dp$altAD)

# Get also the minor allele frequency
dp <- mutate(dp, DP = refAD + altAD, maf = pmin(refAD, altAD)/DP)

cat("Preparing data frame ...\n")
dp <- separate(dp, POS, c("Contig", "pos", "line"), sep = "_", remove = FALSE) # Get details of the site
dp <- separate(dp, sample, c("Environment", "Generation", "Replicate"), sep = "_", remove = FALSE) # Get details of the sample
dp$Coordinate <- paste(dp$Contig, dp$pos, sep = "_")

## Are there sites that are polymorphic in the parentals?
dp %>% filter(sample == "SK1") %>% filter(maf > 0) %>% nrow # 2144
dp %>% filter(sample == "Y55") %>% filter(maf > 0) %>% nrow # 1011
# A few
dp %>% filter(sample == "SK1") %>% filter(maf > 0.05) %>% nrow # 179
dp %>% filter(sample == "Y55") %>% filter(maf > 0.05) %>% nrow # 133

# ============================
# Get the alleles of each site
# ============================
vcffix <- vcf@fix %>% data.frame %>% select(-c(INFO, FILTER, QUAL, ID))
vcffix$POS <- as.numeric(vcffix$POS)
vcffix$Coordinate <- paste(vcffix$CHROM, vcffix$POS, sep = "_")

# Add the REF and ALT data to each variant in each sample
dpf <- merge(dp, vcffix %>% select(Coordinate, REF, ALT), by = "Coordinate", sort = FALSE)

# ============================
# Change the data frame a bit to have the Founder generation be equal to gen 0
# ============================
# First deal with the replicates and make the generation a number
dpf_R <- dpf %>% filter(!Generation %in% c("Founder"))
dpf_R$gen <- gsub("G", "", dpf_R$Generation) %>% as.numeric()

# Now deal with the founder
dpf_founder <- dpf %>% filter(Generation %in% c("Founder"))
dpf_founder$gen <- 0
dpf_founder$Replicate <- "Founder"

# Put back together
dpf2 <- rbind(dpf_R, dpf_founder)

# ============================
# Remove multiallelic sites
# ============================
cat("Removing multiallelic sites ...\n")
## Find repeated sites
# Get the actual coordinate of each site
# dpf2$Coordinate <- paste(dpf2$Contig, dpf2$pos, sep = "_")

dpf_biallelic <- data.frame()
for (samp in dpf2$sample %>% unique) {
  sites <- dpf2 %>% filter(sample == samp)
  # Which sites are those?
  # https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r
  n_occur <- data.frame(table(sites$Coordinate))
  # n_occur[n_occur$Freq > 1,] # Some sites are more than two times (3 and 4 times)
  # How many sites were multiallelic?
  
  # Remove them
  dpf_nomulti <- sites[sites$Coordinate %in% n_occur$Var1[n_occur$Freq == 1],]
  dpf_biallelic <- rbind(dpf_biallelic, dpf_nomulti)
}
# This is constant for each environment so I can just use the last sample
print(paste("Number of multiallelic sites:", length(n_occur$Var1[n_occur$Freq > 1])))
# There shouldn't be any left

# ============================
## Select for sites where the parental strains are homozygous.
# ============================
# A Few sites are not homozygous, but inspecting in IGV shows mostly bad mapping 
# areas or sequencing errors
# dpf_biallelic %>% filter(sample == "SK1") %>% filter(maf > 0) %>% dim
# dpf_biallelic %>% filter(sample == "SK1") %>% filter(maf > 0) %>% ggplot(aes(maf)) + geom_histogram()

HomoSK1 <- dpf_biallelic %>% filter(sample == "SK1") %>% filter(maf < 0.05) %>% .$Coordinate
HomoY55 <- dpf_biallelic %>% filter(sample == "Y55") %>% filter(maf < 0.05) %>% .$Coordinate
heteroFounder <- dpf_biallelic %>% filter(Generation == "Founder", maf > 0.3) %>% .$Coordinate # Truly polymorphic since the founder

# What sites are homozygous in both parents?
homosites <- intersect(HomoSK1, HomoY55)

# Keep only sites that are homozygous in both parents
dpf_homoparents <- dpf_biallelic %>% filter(Coordinate %in% homosites)

# ---
# Get the alleles (nucleotides) of the surviving SNPs for each parental, using 
# the fact that the parental allele should be the major allele
alleleSK1 <- dpf_homoparents %>% filter(sample == "SK1") %>% mutate(alleleSK1 = case_when(refAD > altAD ~ REF, refAD < altAD ~ ALT)) %>% .$alleleSK1
alleleY55 <- dpf_homoparents %>% filter(sample == "Y55") %>% mutate(alleleY55 = case_when(refAD > altAD ~ REF, refAD < altAD ~ ALT)) %>% .$alleleY55

# Because each sample has the same number of sites, the alleles should be repeated and attached correctly
dpf_homoparents <- cbind(dpf_homoparents, alleleSK1, alleleY55) 
dpf_homoparents$pos <- as.numeric(dpf_homoparents$pos) # Otherwise the filtering within winy55perwin() behaves weird and it's silent about it

# Infer the allele frequency of the Y55 allele. However, notice that sites where SK1 == Y55 will still look like Y55's frequency
# This data frame contains any potential de novo mutation
dpf_homoparentsy55_all <- dpf_homoparents %>% mutate(freqY55 = case_when(REF == alleleY55 ~ refAD/DP,
                                               ALT == alleleY55 ~ altAD/DP,
                                               TRUE ~ refAD))

# From those, which ones where polymorphic from the start (if I don't do this, a lot of the final Y55 sites will look "fixed" but that is because Y55=SK1; the disadvantage is that if I there was a mutation in a monomorphic area, I wouldn't see it)
goodsites <- intersect(homosites, heteroFounder) 
dpf_homoparentsy55 <- dpf_homoparentsy55_all %>% filter(Coordinate %in% goodsites)

## Save the files
write.table(dpf_homoparentsy55, file = tabley55snp, sep = "\t", row.names=FALSE, quote=FALSE)
write.table(dpf_homoparentsy55_all, file = tabley55_allsnp, sep = "\t", row.names=FALSE, quote=FALSE)

cat(paste("In the end", dpf_homoparentsy55$Coordinate %>% unique %>% length(), "variants survived."))

# ============================
# Calculate the frequency of one parental in windows (not used in the paper, just exploratory)
# ============================
cat("Calculating MAF in windows ...\n")
chromosomes <- data.frame(Chr = c("ChrI", "ChrII", "ChrIII", "ChrIV", "ChrV", "ChrVI", "ChrVII", "ChrVIII", "ChrIX", "ChrX", "ChrXI", "ChrXII", "ChrXIII", "ChrXIV", "ChrXV", "ChrXVI"),
                          Contig = c("BK006935.2", "BK006936.2", "BK006937.2", "BK006938.2", "BK006939.2", "BK006940.2", "BK006941.2", "BK006934.2", "BK006942.2", "BK006943.2", "BK006944.2", "BK006945.2", "BK006946.2", "BK006947.3", "BK006948.2", "BK006949.2"),
                          len = c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066), 
                          DP = 0, ref = 0, alt = 0, y55 = 0)

### ---------------
# A function to create windows of MAF per windows using a non-overlapping window size
winy55perwin <- function(dfmaf, chrlines, windowsize = 20000, minsnps = 0){
  ## Notice that chrlines is a data frame with two columns, Contig and len, like such:
  #             Contig len
  # 1 chromosome_1   8813524
  # 2 chromosome_2   5165621
  
  # Make empty data frame
  winmaf <- data.frame()
  
  # Do this per sample
  for (samp in dfmaf$sample %>% unique()){
    print(samp)
    samplemaf <- dfmaf %>% filter(sample == samp) # Get the data for only that sample
    
    # Get the metadata (should be the same for all the lines in the dataframe)
    gen = samplemaf[1,]$Generation
    rep = samplemaf[1,]$Replicate
    
    # Calculate the window cov_allele per chromosome
    for (chromi in chrlines$Contig){
      chrsamplemaf <- samplemaf %>% filter(Contig == chromi)
      
      # What's the length of the chromosome?
      chrlen <- chrlines %>% filter(Contig == chromi) %>% .$len
      
      starts <- seq(1, chrlen-windowsize, by = windowsize)
      n <- length(starts)    # Find the length of the vector "starts"
      
      # Loop through each window and get the MAF of SNPs in it
      for (i in 2:n-1){
        currentwin <- chrsamplemaf %>% filter(pos >= starts[i] & pos < starts[i+1])
        chunk <- currentwin$freqY55 %>% median()
        nsnps_win <- currentwin$freqY55 %>% length()
        
        if (nsnps_win >= minsnps) {
          winmaf <- rbind(winmaf, data.frame(Contig = chromi, sample = samp, pos = starts[i], y55 = chunk, Generation = gen, Replicate = rep))
        } else {
          winmaf <- rbind(winmaf, data.frame(Contig = chromi, sample = samp, pos = starts[i], y55 = NA, Generation = gen, Replicate = rep))
        }
        
      }
      
    }
  }
  return(winmaf)
}
### ---------------

## Calculate the MAF in windows and take the median (not used in the paper)
win_dpf_y55 <- winy55perwin(dpf_homoparentsy55, windowsize = 10000, chrlines = chromosomes, minsnps = 2)

# Change names of chromosomes
win_dpf_y55$Contig <- as.factor(win_dpf_y55$Contig)
win_dpf_y55$Contig <- recode_factor(win_dpf_y55$Contig, BK006935.2 = "ChrI", BK006936.2 = "ChrII", BK006937.2 = "ChrIII", BK006938.2 = "ChrIV", BK006939.2 = "ChrV", BK006940.2 = "ChrVI", BK006941.2 = "ChrVII", BK006934.2 = "ChrVIII", BK006942.2 = "ChrIX", BK006943.2 = "ChrX", BK006944.2 = "ChrXI", BK006945.2 = "ChrXII", BK006946.2 = "ChrXIII", BK006947.3 = "ChrXIV", BK006948.2 = "ChrXV", BK006949.2 = "ChrXVI")
chromosomes$Contig <- recode_factor(chromosomes$Contig, BK006935.2 = "ChrI", BK006936.2 = "ChrII", BK006937.2 = "ChrIII", BK006938.2 = "ChrIV", BK006939.2 = "ChrV", BK006940.2 = "ChrVI", BK006941.2 = "ChrVII", BK006934.2 = "ChrVIII", BK006942.2 = "ChrIX", BK006943.2 = "ChrX", BK006944.2 = "ChrXI", BK006945.2 = "ChrXII", BK006946.2 = "ChrXIII", BK006947.3 = "ChrXIV", BK006948.2 = "ChrXV", BK006949.2 = "ChrXVI")

# Get some info for the plot
if ("NaCl" %in% unique(dpf_biallelic$Environment)) {
  thisenv = "NaCl"
} else if ("Ethanol" %in% unique(dpf_biallelic$Environment)){
  thisenv = "Ethanol"
} else if ("LiAc0.01" %in% unique(dpf_biallelic$Environment)){
  thisenv = "LiAc0.01"
} else if ("LiAc0.02" %in% unique(dpf_biallelic$Environment)){
  thisenv = "LiAc0.02"
}

cat("Plotting ...\n")
## Decorate the plot
centromeres <- data.frame(cen = c("CEN1", "CEN2", "CEN3", "CEN4", "CEN5", "CEN6", "CEN7", "CEN8", "CEN9", "CEN10", "CEN11", "CEN12", "CEN13", "CEN14", "CEN15", "CEN16"),
                          start = c(151465, 238207, 114385, 449711, 151987, 148510, 496920, 105586, 355629, 436307, 440129, 150828, 268031, 628758, 326584, 555957),
                          end = c(151582, 238323, 114501, 449829, 152104, 148627, 497038, 105703, 355745, 436425, 440246, 150947, 268149, 628875, 326702, 556073),
                          Contig = c("ChrI", "ChrII", "ChrIII", "ChrIV", "ChrV", "ChrVI", "ChrVII", "ChrVIII", "ChrIX", "ChrX", "ChrXI", "ChrXII", "ChrXIII", "ChrXIV", "ChrXV", "ChrXVI"),
                          y55 = 0)
centromeres <- centromeres %>% mutate(pos = (end - start)/2 + start)

chrlines <- data.frame(pos =c(0), endchr = c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066),
                       Contig = c("ChrI", "ChrII", "ChrIII", "ChrIV", "ChrV", "ChrVI", "ChrVII", "ChrVIII", "ChrIX", "ChrX", "ChrXI", "ChrXII", "ChrXIII", "ChrXIV", "ChrXV", "ChrXVI"), 
                       Strain = character(1))

# To select from a constant color palette
colorreplicates <- c("Founder" = "#000000", "R1" = "#009E73", "R2" = "#0072B2", "R3" = "#E69F00", "R4" = "#D55E00", "R5" = "#CC79A7")

# interestgen <- "G700"
win_dpf_y55gen <- win_dpf_y55 %>% filter(Generation %in% c(interestgen, "Founder")) #%>% filter(Replicate == "R2")
freqplot <- ggplot(win_dpf_y55gen, aes( x = pos, y = y55, colour = Replicate)) +
  facet_grid(Contig ~ .) +
  geom_segment(data = chrlines, # ## Add lines to locate the chromosomes and the centromeres
               aes(x = pos, y = 0, xend = endchr, yend = 0),
               # aes(x = pos, y = Strain, xend = endchr, yend = 1), 
               colour = "gray",
               size = 1, alpha = 0.7) +
  geom_point(data = centromeres, aes(x = pos), colour = "black", size = 1, shape = 19) +
  # geom_line(data = win_dpf_y55 %>% filter(Generation %in% c("Founder")), aes(x = pos, y = maf), colour = "black") +
  geom_line(alpha = 0.8) +
  xlab(expression(paste("Chromosome position (bp)"))) + 
  ylab(expression(paste('Frequency of the Y55 allele'))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(title = paste(interestgen, "in", thisenv)) +
  guides(colour=guide_legend(title=NULL)) + # Remove the legend title
  scale_y_continuous(limits = c(0, 1), breaks=c(0, 0.5, 1.0)) +
  scale_color_manual(values = colorreplicates[win_dpf_y55gen$Replicate %>% unique()] %>% str_sort()) #  

## Save the files
ggsave(ploty55, plot = freqplot, width = 16, height = 22, units = "cm")
write.table(win_dpf_y55, file = tabley55win, sep = "\t", row.names=FALSE, quote=FALSE)
