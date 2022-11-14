#!/usr/bin/env Rscript

# ParallelFixationPlot: Plot much fixation is parallel in the Adaptation dynamics experiment
#############################################################################
# Part of the Snakemake pipeline vcf4adaptation.smk
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022/04/12 - 2022/05/09
#############################################################################
# Load the necessary libraries
# ============================

library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
library(tidyr) # For separate()
library(gtools) # for mixedsort()
library(cowplot)

# ============================
# Check input
# ============================
# Input
windowsfile <- snakemake@input$wins
# Output
G100vsG700file <- snakemake@output$G100vsG700
fixwinplotfile <- snakemake@output$fixwins
G30file <- snakemake@output$G30
NaClfile <- snakemake@output$NaCl

# ============================
# Reading the data
# ============================
allfreqWinAll <- read.table(windowsfile, header=TRUE)

# Let's remove the windows that have NAs
allfreqWinAllclean <- allfreqWinAll[!is.na(allfreqWinAll$Y55),] %>% 
  mutate(Coordinate = paste0(Contig, "_", pos), SK1 = 1-Y55) %>% 
  separate(sample, c("Environment", NA, NA), sep = "_", remove = FALSE) # Recover the environment column

# Sorted sample names
sortedsamples <- allfreqWinAllclean$sample %>% as.factor() %>% levels() %>% mixedsort()
allfreqWinAllclean$sample <- factor(allfreqWinAllclean$sample, levels = sortedsamples)

cat("How many windows survived?\n")
(allfreqWinAllclean %>% filter(sample == "Ethanol_G1000_R1") %>% nrow)*100/(allfreqWinAll %>% filter(sample == "Ethanol_G1000_R1") %>% nrow) 
# 42.56757

# Change the names so they are fancier
allfreqWinAllclean$Environment[allfreqWinAllclean$Environment == "LiAc0.01"] <- "LiAc 0.01M"
allfreqWinAllclean$Environment[allfreqWinAllclean$Environment == "LiAc0.02"] <- "LiAc 0.02M"
# ============================
# Make heatmap of parallelism
# ============================

# Count how many times there is fixation or extinction for the SK1 allele
# Ideally it should be ran for a given environment and a given generation, 
# but one can play with that
countfixation <- function(afdf = allfreqWinAllclean %>% filter(Environment == "NaCl", Generation == "G700")){
  countwin <- data.frame()
  maxnumsamples <- length(unique(afdf$sample))
  
  for (coord in unique(afdf$Coordinate)) { # For every window
    thiswindowdf <- afdf %>% filter(Coordinate == coord)
    count <- 0
    for (pop in thiswindowdf$SK1) {
      if (pop >= 0.9) { count <- count + 1} 
      else if (pop <= 0.1) {count <- count - 1}
      else {count <- count + 0}
    }
    thisline <- cbind(thiswindowdf[1,] %>% select(Contig, Environment, start, end, pos, Generation, Coordinate), parallelraw = count)
    countwin <- rbind(countwin, thisline)
  }
  # Make the counts relative
  countwin$parallel <- countwin$parallelraw/maxnumsamples
  
  # Add an extra column for strictly fixed things
  countwin <- mutate(countwin, fix = case_when(parallel == 1 ~ 1, parallel == -1 ~ -1, TRUE ~ 0))
  # Change names of chromosomes
  countwin$Chromosome <- countwin$Contig
  countwin$Contig <- as.factor(countwin$Contig)
  countwin$Contig <- recode_factor(countwin$Contig, ChrI = "BK006935.2", ChrII = "BK006936.2", ChrIII = "BK006937.2", ChrIV = "BK006938.2", ChrV = "BK006939.2", ChrVI = "BK006940.2", ChrVII = "BK006941.2", ChrVIII = "BK006934.2", ChrIX = "BK006942.2", ChrX = "BK006943.2", ChrXI = "BK006944.2", ChrXII = "BK006945.2", ChrXIII = "BK006946.2", ChrXIV = "BK006947.3", ChrXV = "BK006948.2", ChrXVI = "BK006949.2")
  # For better plotting
  levels(countwin$Chromosome) <- unique(afdf$Contig)
  levels(countwin$Coordinate) <- unique(countwin$Coordinate)
  
  return(countwin)
}

### ----- Parallelism at the environment level -----
relevantenvs_G700 <- filter(allfreqWinAllclean, Generation == "G700") %>% .$Environment %>% unique()

# Produce a dataframe of counts for each environment
parallelenvs_G700 <- data.frame()
for (env in relevantenvs_G700){
  thisenv <- countfixation(afdf = allfreqWinAllclean %>% filter(Environment == env, Generation == "G700"))
  parallelenvs_G700 <- rbind(parallelenvs_G700, thisenv)
}

### ---- Parallelism in LiAc ----
parallelLiAC <- countfixation(afdf = allfreqWinAllclean %>% filter(Generation == "G700", Environment %in% c("LiAc 0.01M", "LiAc 0.02M")))
# The "Environment" is an artifact from the function
parallelLiAC$Environment <- "LiAC"
          
# I need an extra data frame to plot interesting regions in addition to the full plot
parallelLiAC_fix <- rbind(parallelLiAC %>% filter(fix %in% c(-1,1)) %>% mutate(Environment = "LiAc 0.01M"),
                          parallelLiAC %>% filter(fix %in% c(-1,1)) %>% mutate(Environment = "LiAc 0.02M"))

##### ---- A summary plot -----
# I need the windows for marking the chromosomes at the base of the plot
concatchrs <- parallelenvs_G700 %>% filter(Environment == "NaCl") %>% 
  mutate(colchr = case_when(Chromosome %in% unique(parallelenvs_G700$Chromosome)[seq(1,16,2)] ~ "odd",
                            TRUE ~ "pair"))
# concatchrs %>% count(Chromosome, colchr) # How many windows per chromosome?

# Plot the counts as a barcode
plotParaEnvs_G700 <- ggplot(parallelenvs_G700, aes(x = Coordinate, Environment, fill = parallel)) + 
  geom_tile() +
  theme(axis.text.x = element_blank(),
        axis.title.y= element_blank(),
        axis.title.x = element_text(margin = margin(t = -1), size=9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=11),
        axis.ticks.x = element_blank()) +
  geom_point(data = parallelenvs_G700 %>% filter(abs(parallelraw) >= 4), aes(x = Coordinate), color = "black", size = 0.5) +
  geom_point(data = parallelLiAC_fix, aes(x = Coordinate), color = "hotpink1", shape = 3, size = 1) +
  scale_fill_distiller(palette = "BrBG") + # BrBG, PRGn, PuOr, RdBu
  labs(fill='Parallel\n index') + 
  geom_segment(data = concatchrs, 
               aes(x = Coordinate, xend = Coordinate, y = 0.4, yend = 0.5, colour = colchr),
               size = 0.5) +
  scale_colour_manual(values = c("gray", "azure4")) + guides(colour = "none") +
  xlab("Genome coordinate (non-overlapping 10 kb windows)") +
  ggtitle("Generation 700") +
  scale_y_discrete(limits=rev) # To keep Ethanol on top, like in the rest of the paper's figures

#### How about G100, when the first fitness measure happened in the experiment?
### ----- Parallelism at the environment level -----
relevantenvs_G100 <- filter(allfreqWinAllclean, Generation == "G100") %>% .$Environment %>% unique()

# Produce a dataframe of counts for each environment
parallelenvs_G100 <- data.frame()
for (env in relevantenvs_G100){
  thisenv <- countfixation(afdf = allfreqWinAllclean %>% filter(Environment == env, Generation == "G100"))
  parallelenvs_G100 <- rbind(parallelenvs_G100, thisenv)
}

parallelLiAC_G100 <- countfixation(afdf = allfreqWinAllclean %>% filter(Generation == "G100", Environment %in% c("LiAc 0.01M", "LiAc 0.02M")))
# I need an extra data frame to plot interesting regions in addition to the full plot
parallelLiAC_fix_G100 <- rbind(parallelLiAC_G100 %>% filter(fix %in% c(-1,1)) %>% mutate(Environment = "LiAc 0.01M"),
                               parallelLiAC_G100 %>% filter(fix %in% c(-1,1)) %>% mutate(Environment = "LiAc 0.02M"))

plotParaEnvs_G100 <- ggplot(parallelenvs_G100, aes(x = Coordinate, Environment, fill= parallel)) + 
  geom_tile() +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = -1), size=9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=11),
        axis.ticks.x = element_blank()) +
  geom_point(data = parallelenvs_G100 %>% filter(abs(parallelraw) >= 4), aes(x = Coordinate), color = "black", size = 0.5) +
  geom_point(data = parallelLiAC_fix_G100, aes(x = Coordinate), color = "hotpink1", shape = 3, size = 1) +
  scale_fill_distiller(palette = "BrBG") + # BrBG, PRGn, PuOr, RdBu
  geom_segment(data = concatchrs, 
                 aes(x = Coordinate, xend = Coordinate, y = 0.4, yend = 0.5, colour = colchr),
                 size = 0.5) +
  scale_colour_manual(values = c("gray", "azure4")) + guides(colour = "none") +
  labs(fill='Parallel\n index') + 
  ggtitle("Generation 100") +
  xlab("Genome coordinate (non-overlapping 10 kb windows)") +
  scale_y_discrete(limits=rev)

legend <- get_legend( plotParaEnvs_G100 + theme(legend.title = element_text(size=8), 
                                                legend.text = element_text(size=8)))

# Put together
G100vsG700_plots <- plot_grid(plotParaEnvs_G100 + theme(axis.title.x = element_blank(), legend.position = "none"), 
                        plotParaEnvs_G700 + theme(legend.position = "none"), 
                        nrow = 2, align = "v", rel_heights = c(1, 1.05))

G100vsG700 <- plot_grid(G100vsG700_plots,legend, nrow = 1, rel_widths = c(1, 0.1)) # Paper main figure of parallelism
  
ggsave(plot = G100vsG700, 
       filename = G100vsG700file,
       width = 7, height = 4)

######## Distribution of fixed windows
fixedwinscount <- rbind(data.frame(parallelenvs_G100 %>% filter(abs(parallelraw) >= 4) %>% count(Environment), Gen = 100), 
                        data.frame(Environment = "NaCl", n = 0, Gen = 100),
                        data.frame(parallelenvs_G700 %>% filter(abs(parallelraw) >= 4) %>% count(Environment), Gen = 700))

# The probability of an F2 being homozygous for any of the two parental alleles is (0.5*0.5)+(0.5*0.5)=0.5
# The probability of the two F2s being homozygous for the same allele would then be
# (Prob of being homozygous F2_1)*(Prob of being homozygous F2_2)*(Prob_of_being_same_allele)
# = (0.5)*(0.5)*(0.5)
# Equivalently, tracing each of the two individual alleles from the parents, the prob of
# being homozygous for a given allele is 0.25, so the probability of sharing the same allele
# in the two F2s is (0.25*0.25)+(0.25*0.25)=0.125
# So if different replicates fix for a random single genotype, they should share in 
# average 12.5% of their loci. Likewise, the expected shared amount of loci of 4 
# replicates is:
sharedloci4 <- 2*(0.25^4)  # 0.0078125
# And how many windows are there?
nwins <- parallelenvs_G100 %>% filter(Environment == "Ethanol") %>% nrow()
# So the expected number of shared windows is
nwins*sharedloci4 # 3.9375

fixwinplot <- ggplot(fixedwinscount, aes(x = Environment, y = n, fill = factor(Gen))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="bottom") +
  geom_bar(stat = "identity", position="dodge") +
  labs(fill = "Generation") +
  scale_fill_manual(values = c("gray", "azure4")) +
  ylab("Number of windows fixed in four or more replicates") +
  geom_hline(yintercept= round(nwins*sharedloci4), linetype='dashed')

ggsave(plot = fixwinplot, 
       filename = fixwinplotfile,
       width = 3, height = 5)

cat("Coverage of fixed windows in four or more replicates\n")
fixedwinscount %>% mutate(kbs = n * 10)

######## How fast is the red LiAc region fixing?
relevantenvs_G30 <- filter(allfreqWinAllclean, Generation == "G30") %>% .$Environment %>% unique()

# Produce a dataframe of counts for each environment
parallelenvs_G30 <- data.frame()
for (env in relevantenvs_G30){
  thisenv <- countfixation(afdf = allfreqWinAllclean %>% filter(Environment == env, Generation == "G30"))
  parallelenvs_G30 <- rbind(parallelenvs_G30, thisenv)
}
parallelLiAC_G30 <- countfixation(afdf = allfreqWinAllclean %>% filter(Generation == "G30", Environment %in% c("LiAc 0.01M", "LiAc 0.02M")))
# I need an extra data frame to plot interesting regions in addition to the full plot
parallelLiAC_fix_G30 <- rbind(parallelLiAC_G30 %>% filter(abs(parallelraw) >= 8) %>% mutate(Environment = "LiAc 0.01M"),
                              parallelLiAC_G30 %>% filter(abs(parallelraw) >= 8) %>% mutate(Environment = "LiAc 0.02M"))


plotParaEnvs_G30 <- ggplot(parallelenvs_G30, aes(x = Coordinate, Environment, fill= parallel)) + 
    geom_tile() +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(margin = margin(t = -1), size=9),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size=11),
          axis.ticks.x = element_blank()) +
    geom_point(data = parallelenvs_G30 %>% filter(abs(parallelraw) >= 4), aes(x = Coordinate), color = "black", size = 0.5) +
    geom_point(data = parallelLiAC_fix_G30, aes(x = Coordinate), color = "orange", size = 0.5) +
    geom_point(data = parallelLiAC_fix_G100, aes(x = Coordinate), color = "red", size = 0.5) +
    scale_fill_distiller(palette = "BrBG") + # BrBG, PRGn, PuOr, RdBu
    geom_segment(data = concatchrs, 
                 aes(x = Coordinate, xend = Coordinate, y = 0.4, yend = 0.5, colour = colchr),
                 size = 0.5) +
    scale_colour_manual(values = c("gray", "azure4")) + guides(colour = "none") +
    labs(fill='Parallel\n index') + 
    ggtitle("Generation 30") +
    xlab("Genome coordinate (non-overlapping 10 kb windows)") +
    scale_y_discrete(limits=rev)

ggsave(plot = plotParaEnvs_G30, 
       filename = G30file,
       width = 7, height = 2)

######## How much parallelism is there between the two types of replicates in NaCl, 
# the ones that fixed a single genotype, and those that didn't?

# Produce a dataframe of counts for each environment
parallel_NaClG700_MG <- countfixation(afdf = allfreqWinAllclean %>% filter(Environment == "NaCl", 
                                                                           Generation == "G700", 
                                                                           Replicate %in% c("R1", "R4"))) %>% mutate(type = "R1 & R4")

parallel_NaClG700_OG <- countfixation(afdf = allfreqWinAllclean %>% filter(Environment == "NaCl", 
                                                                           Generation == "G700", 
                                                                           Replicate %in% c("R2", "R3"))) %>% mutate(type = "R2 & R3")

parallel_NaClG700_all <- countfixation(afdf = allfreqWinAllclean %>% filter(Environment == "NaCl", 
                                                                            Generation == "G700")) %>% mutate(type = "All NaCl")
# Put them together
parallel_NaClG700 <- rbind(parallel_NaClG700_MG, parallel_NaClG700_OG, parallel_NaClG700_all)


# I need an extra data frame to plot interesting regions in addition to the full plot
parallel_NaClG700_fix <- rbind(parallel_NaClG700_MG %>% filter(abs(parallelraw) >= 8) %>% mutate(Environment = "R1 & R4"),
                               parallel_NaClG700_OG %>% filter(abs(parallelraw) >= 8) %>% mutate(Environment = "R2 & R3"),
                               parallel_NaClG700_all %>% filter(abs(parallelraw) >= 8) %>% mutate(Environment = "All NaCl"))



plotParaNaCl_G700 <- ggplot(parallel_NaClG700, aes(x = Coordinate, type, fill= parallel)) + 
  geom_tile() +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = -1), size=9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=11),
        axis.ticks.x = element_blank()) +
  scale_fill_distiller(palette = "BrBG") + # BrBG, PRGn, PuOr, RdBu
  geom_segment(data = concatchrs, 
               aes(x = Coordinate, xend = Coordinate, y = 0.4, yend = 0.5, colour = colchr),
               size = 0.5) +
  geom_point(data = parallel_NaClG700 %>% filter(abs(parallelraw) >= 4), aes(x = Coordinate), color = "red", size = 0.5) +
  geom_point(data = parallel_NaClG700 %>% filter(type != "All NaCl", abs(parallelraw) >= 2), aes(x = Coordinate), color = "black", size = 0.5) +
  scale_colour_manual(values = c("gray", "azure4")) + guides(colour = "none") +
  labs(fill='Parallel\n index') + 
  ggtitle("Generation 700 for NaCl populations") +
  xlab("Genome coordinate (non-overlapping 10 kb windows)") +
  scale_y_discrete(limits=rev)

ggsave(plot = plotParaNaCl_G700, 
       filename = NaClfile,
       width = 7, height = 2)

cat("Coverage of fixed windows in 2 replicates\n")
parallel_NaClG700 %>% filter(abs(parallelraw) == 2) %>% count(type) %>% mutate(kbs = n * 10, per = (n/nwins)*100)

### Number of expected shared windows by chance for two F2s
# For a given locus with two alleles, A and B, the probability of an F2 to be 
# homozygous for one allele is 0.5*0.5 = 0.25, and for the two possible
# homozygous cases (AA and BB) is 0.25*2 = 0.5. The probability of two F2s to be 
# homozygous for the same allele would then be: 
# (Prob of homozygous first F2)*(Prob of homozygous second F2)*(prob of same allele)
#  = 0.5*0.5*0.5 = 0.125
0.125*nwins # 63

# How many windows are overlapping in the two types of populations?
setR1R4 <- parallel_NaClG700 %>% filter(abs(parallelraw) == 2, type == "R1 & R4") %>% .$Coordinate
setR2R3 <- parallel_NaClG700 %>% filter(abs(parallelraw) == 2, type == "R2 & R3") %>% .$Coordinate

intersect(setR1R4, setR2R3) %>% length # 16
# So both conditions fix 129 windows but 16 of those are shared!
