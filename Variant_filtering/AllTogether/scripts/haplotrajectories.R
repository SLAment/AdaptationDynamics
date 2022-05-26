#!/usr/bin/env Rscript

# haplotrajectories: plot the standing genetic variation clusters and mutations
#############################################################################
# Part of the Snakemake pipeline vcf4adaptation.smk 
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022-01-21 to 2022-04-25
# Version 1
#############################################################################
# Load the necessary libraries
# ============================
library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
library(vcfR)
library(tidyr) # For separate()

# ============================
# Set paths
# ============================
haplos_file <- snakemake@input$clusters
trajectos_file <-  snakemake@input$rawsnps
denovomutations_file <- snakemake@input$mutations

# ============================
# Reading the data
# ============================
haplos <- read.table(haplos_file, header=TRUE) %>% 
  mutate(Coordicon = paste0(Coordinate, "_", Condition)) %>% 
  mutate(clust_id = paste0(Condition, "_", clust))

allfreq_set <- read.table(trajectos_file, header=TRUE) %>% 
  mutate(Condition = paste0(Environment, "_", Replicate)) %>% 
  mutate(Coordicon = paste0(Coordinate, "_", Condition)) 

denovomutations <- read.table(denovomutations_file, header=TRUE) %>% 
  mutate(Condition = paste0(Environment, "_", Replicate)) %>% 
  mutate(Coordicon = paste0(Coordinate, "_", Condition)) %>% 
  select(names(allfreq_set))

# Attach the missing info to the selected sites
selectsites <- merge(haplos, allfreq_set, by = c("Coordinate", "Condition", "Coordicon"))
# Reorder columns
selectsites <- select(selectsites, c(colnames(allfreq_set), "clust", "N_var", "clust_id"))

# But then there are haplotypes that are larger than others. We used a threshold
# of 10 SNPs to find the "big boi" haplotypes
bigbois <- filter(selectsites, N_var >= 10)
notbigbois <- filter(selectsites, N_var < 10)

# Put them back together
trajectories <- rbind(bigbois, notbigbois)

# -----
# Plot with the colors of Ciaran's palette
CiaransPalette <- c("NaCl_R1" = "#87CEFA", "NaCl_R2" = "#5AA3CF", "NaCl_R3" = "#1E90FF", "NaCl_R4" = "#4169E1",
                    "Ethanol_R1" = "#CD5C5C", "Ethanol_R2" = "#FF0000", "Ethanol_R3" = "#AC1016", "Ethanol_R4" = "#67000C",
                    "LiAc0.01_R1" = "#D8BFD8", "LiAc0.01_R2" = "#DA70D6", "LiAc0.01_R3" = "#9400D3", "LiAc0.01_R4" = "#800080", "LiAc0.01_R5" = "#4B0082",
                    "LiAc0.02_R1" = "#98D493", "LiAc0.02_R2" = "#9ACD32", "LiAc0.02_R3" = "#228B22", "LiAc0.02_R4" = "#6B8E23", "LiAc0.02_R5" = "#00441B",
                    "Ethanol" = "#AC1016", "NaCl" = "#4169E1", "LiAc0.01" = "#800080", "LiAc0.02" = "#228B22" )

### Sample a random site per cluster per replicate
sampleEnv <- function(env){
  envsiteclust <- bigbois %>% filter(Environment == env) %>% select(Coordinate, clust_id) %>% distinct()
  # Split by Replicate and cluster, and then sample one site
  set.seed(123)
  envsiteclust_sample <- do.call(rbind, lapply( split(envsiteclust, envsiteclust$clust_id), function(envsiteclust) envsiteclust[sample(nrow(envsiteclust), 1) , ] ) )
  bigbois_envsamp <- bigbois %>% filter(Environment == env) %>% filter(Coordinate %in% envsiteclust_sample$Coordinate & clust_id %in% envsiteclust_sample$clust_id)
  return(bigbois_envsamp)
}

bigbois_sampled <- rbind(sampleEnv("NaCl"),
      sampleEnv("Ethanol"),
      sampleEnv("LiAc0.01"),
      sampleEnv("LiAc0.02") )

## The full experiment with the de novo mutations # Paper figure!!
sampledhaplo <- ggplot(bigbois_sampled %>% filter(!sample %in% c("SK1", "Y55")), 
       aes(x = gen, y = 1 - freqY55, colour = as.factor(clust), fill = Coordicon)) + 
  geom_line(alpha = 0.5, size = 0.8) +
  geom_line(data = denovomutations, colour = "black", alpha = 0.8) +
  facet_grid(Environment~Replicate) +
  xlab("Generation") + ylab('Frequency') + # the frequency of the SK1 allele or the de novo mutation
  scale_x_continuous(breaks = trajectories$gen %>% unique %>% sort,
                     labels = c(0, "", 60, "", 200, 300, 400, 500, 700, 1000)) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title = element_text(size=16),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Save the plot
ggsave(plot = sampledhaplo, 
       file = snakemake@output$sgvplot,
       width = 8.5, height = 6)

### How many bigboi clusters are there per environment? (not in the paper)
bigboistable <- bigbois_sampled %>% select(Condition, clust) %>% distinct() %>% 
  group_by(Condition) %>% 
  summarise(NoBigBois = length(clust))

propbiboi <- ggplot(bigboistable, aes(x = NoBigBois, y = Condition, fill = Condition)) + 
  geom_bar(stat="identity") + 
  xlab("Number of haplotypes with > 10 sites") +
  ylab("Replicate") +
  theme_bw() +
  scale_fill_manual(values = CiaransPalette) + guides(fill = "none") +
  scale_y_discrete(limits = rev)

ggsave(plot = propbiboi, 
       file = snakemake@output$bigbois,
       width = 4, height = 6)

