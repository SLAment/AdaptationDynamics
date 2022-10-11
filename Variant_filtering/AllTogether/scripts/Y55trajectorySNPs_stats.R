#!/usr/bin/env Rscript

# Y55trajectorySNPs_stats: script to characterize the loss of genetic diversity 
# in the experiment
#############################################################################
# Part of the Snakemake pipeline vcf4adaptation.smk 
# The name Y55 just makes a reference to the fact that I track one of the parents
# but in the end I plot the SK1 allele in the paper
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022-01-06 - 2022-05-26
# Version 1
#############################################################################
# Load the necessary libraries
# ============================
library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
library(psych) # for harmonic and geometric means
library(spatstat) # for weighted median
library(sicegar) # to fit a sigmoidal model
library(ggpubr)
library(rstatix) # for wilcox_test
library(cowplot)

# ============================
# Get files
# ============================

# Snakemake
tablemaf_file <- snakemake@input$table
# De novo mutations
denovomutations_file <- snakemake@input$mutations
# Read coverage data
coveragesnpsfile <- snakemake@input$coverage

# ============================
# Reading the data
# ============================
allfreq_all <- read.table(tablemaf_file, header=TRUE)

allfreq_NaCl <- filter(allfreq_all, Environment == "NaCl")
allfreq_EtOH <- filter(allfreq_all, Environment == "Ethanol")
allfreq_LiAc01 <- filter(allfreq_all, Environment == "LiAc0.01")
allfreq_LiAc02 <- filter(allfreq_all, Environment == "LiAc0.02")

denovomutations <- read.table(denovomutations_file, header=TRUE)

coveragesnps <- read.table(coveragesnpsfile, header=TRUE)

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
    new_allfreq_env <- rbind(new_allfreq_env, allfreq_env_F)
  }
  
  return(new_allfreq_env)
}

# put them together
allfreq <- rbind(fixfounder(allfreq_EtOH),
                 fixfounder(allfreq_NaCl),
                 fixfounder(allfreq_LiAc01),
                 fixfounder(allfreq_LiAc02)) 

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

commonset <- intersect(LEfounder$Coordinate, Nfounder$Coordinate)
allfreq_set <- filter(allfreq, Coordinate %in% commonset) # %>% dplyr::select(-POS)

# ============================
# Putting together the trajectories of the "all-samples" set with the de novo mutations
# ============================
# These are included as supplementary figures in the paper
allfreq_set_novo  <- rbind(cbind(allfreq_set, class = "standing_variation"), 
                           cbind(select(denovomutations, names(allfreq_set)), class = "de_novo" ) )

allelesInTime_env_mut <- function(env = "NaCl"){
  p <- ggplot(allfreq_set_novo %>% filter(!sample %in% c("SK1", "Y55")) %>% 
                filter(Environment == env, class == "standing_variation"), 
              aes(x = gen, y = 1 - freqY55, colour = Coordinate)) + 
    geom_line(alpha = 0.1) +
    geom_line(data = allfreq_set_novo %>% filter(Environment == env, class == "de_novo"), aes(x = gen, y = 1 - freqY55, fill = Coordinate), alpha = 1, colour = "black") +
    #  geom_line(alpha = 0.1, colour = "black") + # Have the raw SNPs in black and mutation in colour
    # geom_line(data = allfreq_set_novo %>% filter(Environment == env, class == "de_novo"), aes(x = gen, y = 1 - freqY55, colour = Coordinate), alpha = 1, size = 0.8) +
    # labs(title = env) +
    facet_grid(Contig~Replicate) +
    xlab("Generation") + ylab(paste0('Allele Frequency in ', env)) +
    scale_x_continuous(breaks = allfreq_set$gen %>% unique %>% sort,
                       labels = c(0, "", 60, "", 200, 300, 400, 500, 700, 1000)) +
    theme_bw() + 
    theme(legend.position = "none",
          axis.title = element_text(size=16),
          strip.text.x = element_text(face = "bold", size = 12),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}

# Snakemake
ggsave(plot = allelesInTime_env_mut("NaCl"), file = snakemake@output$nacl , width = 7, height = 11)
ggsave(plot = allelesInTime_env_mut("Ethanol"), file = snakemake@output$ethanol , width = 7, height = 11)
ggsave(plot = allelesInTime_env_mut("LiAc0.01"), file = snakemake@output$li1 , width = 7, height = 11)
ggsave(plot = allelesInTime_env_mut("LiAc0.02"), file = snakemake@output$li2 , width = 7, height = 11)

# ============================
# What percentage of sites are fixed at the end of the experiment?
# ============================
cat("No. of SNPs per sample:", length(commonset), "\n")

# For this I don't need the founders, since they are all supposed to be 100% heterozygous
# Let's double check that
# filter(LEfounder, freqY55 >= 0.9 | freqY55 <= 0.1) # 0
# filter(Nfounder, freqY55 >= 0.9 | freqY55 <= 0.1) # 0

propfixdf <- data.frame()
for (env in unique(allfreq_set$Environment)) {
  thisenv <- allfreq_set %>% filter(Environment == env)
  for (r in unique(thisenv$Replicate)) {
    for (g in unique(thisenv$gen)) {
      thissample <- thisenv %>% filter(gen == g, Replicate == r)
      
      if(nrow(thissample) != 0){ # Does this exist? some replicates lost G1000
        nearfix_genome <- filter(thissample, freqY55 >= 0.9 | freqY55 <= 0.1)
        propfixed_genome <- length(nearfix_genome$Coordinate)/length(thissample$Coordinate)
      
        numY55 <- filter(thissample, freqY55 >= 0.9) %>% nrow # Of the nearly fixed ones, how many are Y55?
        
        # Save the genomewide proportion
        propfixdf <- rbind(propfixdf, data.frame(Contig = "genome", thissample[1,] %>% dplyr::select(sample, Environment, Generation, gen, Replicate), ntotal = nrow(thissample), numfix = nrow(nearfix_genome), propfix = propfixed_genome, numY55 = numY55, propY55 = numY55/nrow(nearfix_genome)))
        
         # But we can get it per chromosome too
        for (cont in unique(thissample$Contig) %>% sort) {
          thissample_chr <- thissample %>% filter(Contig == cont)
          nearfix_cont <- filter(thissample_chr, freqY55 >= 0.9 | freqY55 <= 0.1)
          propfixed_cont <- length(nearfix_cont$Coordinate)/length(thissample_chr$Coordinate)
          
          numY55_chr <- filter(thissample_chr, freqY55 >= 0.9) %>% nrow # Of the nearly fixed ones, how many are Y55?
          
          propfixdf <- rbind(propfixdf, data.frame(thissample_chr[1,] %>% dplyr::select(Contig, sample, Environment, Generation, gen, Replicate), ntotal = nrow(thissample_chr), numfix = nrow(nearfix_cont), propfix = propfixed_cont, numY55 = numY55_chr, propY55 = numY55_chr/nrow(nearfix_cont)))
        }
      }
    }
  }
}

# Plot with the colors of Ciaran's palette
CiaransPalette <- c("NaCl_R1" = "#87CEFA", "NaCl_R2" = "#5AA3CF", "NaCl_R3" = "#1E90FF", "NaCl_R4" = "#4169E1",
                    "Ethanol_R1" = "#CD5C5C", "Ethanol_R2" = "#FF0000", "Ethanol_R3" = "#AC1016", "Ethanol_R4" = "#67000C",
                    "LiAc0.01_R1" = "#D8BFD8", "LiAc0.01_R2" = "#DA70D6", "LiAc0.01_R3" = "#9400D3", "LiAc0.01_R4" = "#800080", "LiAc0.01_R5" = "#4B0082",
                    "LiAc0.02_R1" = "#98D493", "LiAc0.02_R2" = "#9ACD32", "LiAc0.02_R3" = "#228B22", "LiAc0.02_R4" = "#6B8E23", "LiAc0.02_R5" = "#00441B",
                    "Ethanol" = "#AC1016", "NaCl" = "#4169E1", "LiAc0.01" = "#800080", "LiAc0.02" = "#228B22" )


# So I can tell them apart more easily
propfixdf <- propfixdf %>% mutate(Condition = paste0(Environment, "_", Replicate))

# ============================
# Getting some descriptive numbers - sigmoidal curve 
# ============================
## Make a function to extract the key values of a sigmoidal fit
# https://cran.r-project.org/web/packages/sicegar/vignettes/fitting_individual_models.html
getsigmoidal <- function(df, condition){
  subdf <- df %>% filter(Condition == condition, Contig == "genome")
  rawdata <- data.frame(time = subdf$gen, intensity = subdf$propfix) # The dataframe has to have those names
  
  # Normalize data
  normalizedInput <- normalizeData(dataInput = rawdata, 
                                   dataInputName = condition)
  # Do the sigmoidal fit
  sigmoidalModel <- multipleFitFunction(dataInput=normalizedInput,
                                        model="sigmoidal")
  # Get additional parameters
  sigmoidalModelAugmented <- parameterCalculation(sigmoidalModel)
  
  # --- Get the fitted points ---
  fittedXmax_sigmoidal <- sigmoidalModelAugmented$dataScalingParameters.timeRange
  time <- seq(0, fittedXmax_sigmoidal, fittedXmax_sigmoidal/1000)
  intensityTheoreticalSigmoidal <- sicegar::sigmoidalFitFormula(time, 
                                                                maximum = sigmoidalModelAugmented$maximum_y, 
                                                                slopeParam = sigmoidalModelAugmented$slopeParam_Estimate, 
                                                                midPoint = sigmoidalModelAugmented$midPoint_x)
  intensityTheoreticalSigmoidalDf <- data.frame(gen = time, propfix = intensityTheoreticalSigmoidal, Condition = condition, Environment = unique(subdf$Environment))
  # ----------
  report <- data.frame(
    Condition = condition, 
    Environment = subdf$Environment[1],
    Replicate = subdf$Replicate[1],
    maximum_y = sigmoidalModelAugmented$maximum_y, 
    midPoint_x =sigmoidalModelAugmented$midPoint_Estimate,
    slope =sigmoidalModelAugmented$slope,
    reachMaximum_x =sigmoidalModelAugmented$reachMaximum_x
  )
  # maximum_y: The maximum intensity the fitted curve reaches at infinity. The value is equal to maximum_Estimate. 
  # midPoint_x: The x value (i.e., time) at which the fitted curve reaches the midpoint. The value is equal to midPoint_Estimate. # For some reason comes NA so I just used midPoint_Estimate
  # slope: The maximum slope of the fitted curve. This is the slope at the midpoint. The value is equal to slopeParam_Estimate * maximum_y / 4.
  # reachMaximum_x: The x value (i.e., time) of the reach maximum point. The reach maximum point is defined as the point where the slope tangent intersects with y = maximum_y. It approximately represents the moment in time when the intensity signal reaches its maximum. Its value is equal to midPoint_x + (incrementTime/2).
  
  ## -- output
  out <- list()
  
  out$report <- report
  out$curve <- intensityTheoreticalSigmoidalDf
  return(out)
}

## Make a dataframe with all conditions with fit parameters
propfixsig <- data.frame()
sigmoidcurve <- data.frame()
for (con in unique(propfixdf$Condition)) {
  sigmoidalfit <- getsigmoidal(propfixdf, con) # calculate the curve
  propfixsig <- rbind(propfixsig, sigmoidalfit$report) # for comparing the main parameters
  sigmoidcurve <- rbind(sigmoidcurve, sigmoidalfit$curve) # for plotting
}

## Let's plot the fitted curve and raw values as points
propfixp <- ggplot(propfixdf %>% filter(Contig == "genome"),
                    aes(x = gen, y = propfix, colour = Condition)) + 
    geom_line(data = sigmoidcurve, size = 1) +
    geom_point(alpha = 0.7) +
    # geom_point(aes(shape = Replicate), size = 2) +
    facet_grid(Environment~.) +
    xlab("Generation") + ylab('Proportion of sites near fixation') +
    scale_x_continuous(breaks = allfreq_set$gen %>% unique %>% sort,
                       labels = c(0, "", 60, "", 200, 300, 400, 500, 700, 1000)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          # strip.background = element_blank()) +
          strip.background = element_rect(fill="gray94")) +
    scale_colour_manual(values = CiaransPalette) + guides(colour = "none")

#Are the differences significant?
kruskal.test(maximum_y ~ Environment, data = propfixsig)
# Kruskal-Wallis chi-squared = 10.549, df = 3, p-value = 0.01443
# we reject the hypothesis that all means are equal

kruskal.test(reachMaximum_x ~ Environment, data = propfixsig)
# Kruskal-Wallis chi-squared = 11.325, df = 3, p-value = 0.01009

# How is different?
propfixsig %>% pairwise_wilcox_test(maximum_y ~ Environment, p.adjust.method = "bonferroni")
propfixsig %>% pairwise_wilcox_test(reachMaximum_x ~ Environment, p.adjust.method = "bonferroni")
# Nothing is significant after correction :(

propfixsig %>% group_by(Environment) %>% summarise(var_maxy = var(maximum_y), var_reach = var(reachMaximum_x))

# Time to max proportion of fixed sites for each environment
timetomax <- ggplot(propfixsig, aes(x = Environment, y = reachMaximum_x, colour = Environment)) + 
    geom_boxplot() + 
    geom_jitter(size = 3, aes(colour = Condition), width = 0.05) +
    ylab("Time to maximum proportion \nof fixed sites (generations)") +
    theme_bw() +
    ylim(0,400) +
    theme(legend.position = "none") +
    scale_colour_manual(values = CiaransPalette)

# Slope
slopep <- ggplot(propfixsig, aes(x = Environment, y = slope, colour = Environment)) + 
  geom_boxplot() + 
  geom_jitter(size = 3, aes(colour = Condition), width = 0.05) +
  ylab("Slope of the fitted \nsigmoid curve") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_manual(values = CiaransPalette)

# maximum_y
maxyp <- ggplot(propfixsig, aes(x = Environment, y = maximum_y, colour = Environment)) + 
    geom_boxplot() + 
    geom_jitter(size = 3, aes(colour = Condition), width = 0.05) +
    ylab("Maximum proportion of \nnearly fixed sites") +
    xlab("") +
    theme_bw() +
    expand_limits(y=0) + # Only affect lower bound, otherwise the points at 1 are lost
    theme(legend.position = "none") +
    scale_colour_manual(values = CiaransPalette)

rightside <- plot_grid(maxyp, timetomax, nrow = 2, labels = c('B', 'C'))
(sigmoidallplots <- plot_grid(propfixp, rightside, ncol = 2, rel_widths = c(2.5, 1.5), labels = c("A")))

ggsave(plot = sigmoidallplots, file = snakemake@output$sigmoid, width = 8.5, height = 6)

# ============================
# Is there a correlation with coverage and proportion of fixed sites?
# ============================
propfixdfcov <- merge(propfixdf, coveragesnps, by = "sample")

# Milo Johnson asked: Do time points where previously fixed SNPs go down generally have low coverage? 

# Sort by generation so I can use a loop
propfixdfcovgenome <- propfixdfcov %>% filter(Contig == "genome") %>% arrange(gen)

propfixdf_delta <- data.frame()
for (con in unique(propfixdfcovgenome$Condition)) {
  thiscon <- propfixdfcovgenome %>% filter(Condition == con)
  gens <- unique(thiscon$gen)
  for (sam in thiscon$sample) {
    thisample <- thiscon %>% filter(sample == sam)
    indexgen <- which(gens == thisample$gen)
    if (indexgen == 1) {
      deltapropfix <- 0
      deltacov <- 0
    } else {
      previousample <- thiscon %>% filter(gen == gens[indexgen - 1])
      deltapropfix <- thisample$propfix - previousample$propfix
      deltacov <- thisample$median - previousample$median
    }
    propfixdf_delta <- rbind(propfixdf_delta, cbind(thisample, deltapropfix = deltapropfix, deltacov = deltacov))
  }
}

# Samples where the coverage was smaller and the proportion of nearly fixed sites was bigger in t-1 
worrisomesamples <- propfixdf_delta %>% filter(deltapropfix < 0) %>% .$sample
badsamples <- propfixdf_delta %>% filter(deltacov > 0 & deltapropfix < 0) %>% .$sample

propfixdf_delta <- propfixdf_delta %>% 
  mutate(pointcolor = case_when(deltacov > 0 & deltapropfix < 0 ~ "bad", deltapropfix < 0 ~ "worry", TRUE ~ "ok"))

# Now correlate the change in fixation proportion and coverage at time t-1
deltapropfixVscoverage <- ggplot(propfixdf_delta %>% filter(Contig %in% c("genome"), Generation != "Founder"), aes(x = median - deltacov, y = deltapropfix)) +
  geom_point(size = 3, alpha = 0.7, aes(colour = pointcolor)) + theme_bw() +
  geom_hline(yintercept=0, colour = "gray") +
  scale_color_manual(values=c("ok" = "black", "bad" = "red", "worry" = "darkcyan")) +
  xlab("Median depth of coverage at time t-1") + ylab(bquote(Delta~'Fixation proportion = '~Fix[t]~'-'~Fix[t-1])) +
  facet_grid(Environment ~ ., scales = "free") +
  stat_cor(method="pearson") + guides(colour = "none")

# And what about the difference in coverage between time points?
deltapropfixVsdeltacoverage <- ggplot(propfixdf_delta %>% filter(Contig %in% c("genome"), Generation != "Founder"), aes(x = deltacov, y = deltapropfix)) +
  geom_point(size = 3, alpha = 0.7, aes(colour = pointcolor)) + theme_bw() + 
  geom_hline(yintercept=0, colour = "gray") +
  geom_vline(xintercept=0, colour = "gray") +
  scale_color_manual(values=c("ok" = "black", "bad" = "red", "worry" = "darkcyan")) +
  xlab(bquote(Delta~'depth of coverage = '~Cov[t]~'-'~Cov[t-1])) +
  ylab("") +
  facet_grid(Environment ~ ., scales = "free") + 
  stat_cor(method="pearson") + guides(colour = "none")

## Figure in the report
ggsave(plot_grid(deltapropfixVscoverage, deltapropfixVsdeltacoverage, ncol = 2, labels = "AUTO", align = "h"),
       file = snakemake@output$covcorr, width = 8, height = 8)

