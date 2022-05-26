#!/usr/bin/env Rscript

# ParallelFixation: How much fixation is parallel in the Adaptation dynamics experiment
#############################################################################
# Part of the Snakemake pipeline vcf4adaptation.smk
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022/04/12 - 2022/05/09
#############################################################################
# Load the necessary libraries
# ============================

library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
library(tidyr) # For separate()

# ============================
# Check input
# ============================
# Input
trajectos_file <- snakemake@input$table
denovomutations_file <- snakemake@input$mutations
# Output
windowsfile <- snakemake@output$wins
G100vsG700file <- snakemake@output$G100vsG700
fixwinplotfile <- snakemake@output$fixwins
G30file <- snakemake@output$G30

# ============================
# Reading the data
# ============================
allfreq_set <- read.table(trajectos_file, header=TRUE)

denovomutations <- read.table(denovomutations_file, header=TRUE) %>% 
  mutate(Coordinate = paste0(Contig, '_', pos))

# ============================
# Calculate the frequency of one parental in windows
# ============================
cat("Calculating parental allele frequency in windows ...\n")

## Can I add the chromosomes?
chromosomes <- data.frame(Contig = c("ChrI", "ChrII", "ChrIII", "ChrIV", "ChrV", "ChrVI", "ChrVII", "ChrVIII", "ChrIX", "ChrX", "ChrXI", "ChrXII", "ChrXIII", "ChrXIV", "ChrXV", "ChrXVI"),
                          # Contig = c("BK006935.2", "BK006936.2", "BK006937.2", "BK006938.2", "BK006939.2", "BK006940.2", "BK006941.2", "BK006934.2", "BK006942.2", "BK006943.2", "BK006944.2", "BK006945.2", "BK006946.2", "BK006947.3", "BK006948.2", "BK006949.2"),
                          len = c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066), 
                          pos = 0, Environment = character(1), parallel = 0)

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
    gen <- samplemaf$Generation[1]
    rep <- samplemaf$Replicate[1]
    
    # Calculate the window median allele frequency per chromosome
    for (chromi in chrlines$Contig){
      chrsamplemaf <- samplemaf %>% filter(Contig == chromi)
      
      # What's the length of the chromosome?
      chrlen <- chrlines %>% filter(Contig == chromi) %>% .$len
      
      starts <- seq(1, chrlen-windowsize, by = windowsize)
      n <- length(starts)    # Find the length of the vector "starts"
      
      # Loop through each window and get the parentall allele frequency of SNPs in it
      for (i in 2:n-1){
        currentwin <- chrsamplemaf %>% filter(pos >= starts[i] & pos < starts[i+1])
        chunk <- currentwin$freqY55 %>% median()
        nsnps_win <- currentwin$freqY55 %>% length()
        meanpos <- windowsize/2 + starts[i]
        
        if (nsnps_win >= minsnps) {
          winmaf <- rbind(winmaf, data.frame(Contig = chromi, sample = samp, start = starts[i], end = starts[i+1], pos = meanpos, Y55 = chunk, Generation = gen, Replicate = rep))
        } else {
          winmaf <- rbind(winmaf, data.frame(Contig = chromi, sample = samp, start = starts[i], end = starts[i+1], pos = meanpos, Y55 = NA, Generation = gen, Replicate = rep))
        }
        
      }
      
    }
  }
  return(winmaf)
}
### ---------------

cat("Calculating the parental allele frequency in windows and take the median (it'll take a long time!!)\n")
allfreqWinAll <- winy55perwin(allfreq_set, windowsize = 10000, chrlines = chromosomes, minsnps = 4)
# Save this precious file
write.table(allfreqWinAll, file = windowsfile, sep = "\t", row.names=FALSE, quote=FALSE)
