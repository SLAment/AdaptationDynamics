#!/usr/bin/env Rscript

# DeNovoMutations: Trajectory of de novo mutations in the adaptation experiment
#############################################################################
# Part of the Snakemake pipeline vcf4adaptation.smk 
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2021-11-01 -- 2022-03-28
# Version 1
#############################################################################
# Load the necessary libraries
# ============================
library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff
library(ggpubr) # for get_palette
# library(cowplot)

# ============================
# Read file names
# ============================
## Snakemake
# SNPs
tablemaf_EtOH_snp <- snakemake@input$EtOH_snp
tablemaf_NaCl_snp <- snakemake@input$NaCl_snp
tablemaf_LiAc01_snp <- snakemake@input$LiAc01_snp
tablemaf_LiAc02_snp <- snakemake@input$LiAc02_snp

# # INDELs
tablemaf_EtOH_indel <- snakemake@input$EtOH_indel
tablemaf_NaCl_indel <- snakemake@input$NaCl_indel
tablemaf_LiAc01_indel <- snakemake@input$LiAc01_indel
tablemaf_LiAc02_indel <- snakemake@input$LiAc02_indel

# SnpEff data
snpefffile_EtOH_snp <- snakemake@input$EtOH_snp_snpeff
snpefffile_NaCl_snp <- snakemake@input$NaCl_snp_snpeff
snpefffile_LiAc01_snp <- snakemake@input$LiAc01_snp_snpeff
snpefffile_LiAc02_snp <- snakemake@input$LiAc02_snp_snpeff

snpefffile_EtOH_indel <- snakemake@input$EtOH_indel_snpeff
snpefffile_NaCl_indel <- snakemake@input$NaCl_indel_snpeff
snpefffile_LiAc01_indel <- snakemake@input$LiAc01_indel_snpeff
snpefffile_LiAc02_indel <- snakemake@input$LiAc02_indel_snpeff

# Bad mutations
bad_mutations_file <- snakemake@input$bad_mutations

# ============================
# Reading the data
# ============================

## Make them tables
## SNPs
# Read data of allele frequencies
allfreq_EtOH_snps <- read.table(tablemaf_EtOH_snp, header=TRUE)
allfreq_NaCl_snps <- read.table(tablemaf_NaCl_snp, header=TRUE)
allfreq_LiAc01_snps <- read.table(tablemaf_LiAc01_snp, header=TRUE)
allfreq_LiAc02_snps <- read.table(tablemaf_LiAc02_snp, header=TRUE)

# Read data of snpEff
snpeff_EtOH_snps <- read.table(snpefffile_EtOH_snp, header=TRUE) %>% mutate(Environment = "Ethanol", Coordinate = paste0(Contig, '_', position))
snpeff_NaCl_snps <- read.table(snpefffile_NaCl_snp, header=TRUE) %>% mutate(Environment = "NaCl", Coordinate = paste0(Contig, '_', position))
snpeff_LiAc01_snps <- read.table(snpefffile_LiAc01_snp, header=TRUE) %>% mutate(Environment = "LiAc0.01", Coordinate = paste0(Contig, '_', position))
snpeff_LiAc02_snps <- read.table(snpefffile_LiAc02_snp, header=TRUE) %>% mutate(Environment = "LiAc0.02", Coordinate = paste0(Contig, '_', position))

# Put together the snpeff data
allsnpeff_snps <- rbind(snpeff_EtOH_snps, snpeff_NaCl_snps, snpeff_LiAc01_snps, snpeff_LiAc02_snps)

## INDELs
# Read data of allele frequencies
allfreq_EtOH_indels <- read.table(tablemaf_EtOH_indel, header=TRUE)
allfreq_NaCl_indels <- read.table(tablemaf_NaCl_indel, header=TRUE)
allfreq_LiAc01_indels <- read.table(tablemaf_LiAc01_indel, header=TRUE)
allfreq_LiAc02_indels <- read.table(tablemaf_LiAc02_indel, header=TRUE)

# Put together the snpeff data
snpeff_EtOH_indels <- read.table(snpefffile_EtOH_indel, header=TRUE) %>% mutate(Environment = "Ethanol", Coordinate = paste0(Contig, '_', position))
snpeff_NaCl_indels <- read.table(snpefffile_NaCl_indel, header=TRUE) %>% mutate(Environment = "NaCl", Coordinate = paste0(Contig, '_', position))
snpeff_LiAc01_indels <- read.table(snpefffile_LiAc01_indel, header=TRUE) %>% mutate(Environment = "LiAc0.01", Coordinate = paste0(Contig, '_', position))
snpeff_LiAc02_indels <- read.table(snpefffile_LiAc02_indel, header=TRUE) %>% mutate(Environment = "LiAc0.02", Coordinate = paste0(Contig, '_', position))

# Now the snpwff data
allsnpeff_indels <- rbind(snpeff_EtOH_indels, snpeff_NaCl_indels, snpeff_LiAc01_indels, snpeff_LiAc02_indels)

# Absolutely all snpeff together
allsnpeff <- rbind(data.frame(allsnpeff_snps, class = "SNP"),
                   data.frame(allsnpeff_indels, class = "INDEL"))


bad_mutations <- read.table(bad_mutations_file, header = FALSE)$V1 # A list of mutations to ignore because manual curation showed they are not real/reliable

# ============================
# Define the thresholds for the filtering and plot
# ============================
#############
minfreqmutation <- 0.10
selectedfreq <- 0.65
#############
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

# put SNPs together
allfreq_snps <- rbind(fixfounder(allfreq_EtOH_snps),
                      fixfounder(allfreq_NaCl_snps),
                      fixfounder(allfreq_LiAc01_snps),
                      fixfounder(allfreq_LiAc02_snps))

# put INDELs together
allfreq_indels <- rbind(fixfounder(allfreq_EtOH_indels),
                        fixfounder(allfreq_NaCl_indels),
                        fixfounder(allfreq_LiAc01_indels),
                        fixfounder(allfreq_LiAc02_indels))

# Now together
allfreq <- rbind(data.frame(allfreq_snps, class = "SNP"),
                 data.frame(allfreq_indels, class = "INDEL"))

# ============================
# Explore the Y55 allele frequency changes in time
# ============================
# Change names of chromosomes
allfreq$Contig <- as.factor(allfreq$Contig)
allfreq$Contig <- recode_factor(allfreq$Contig, BK006935.2 = "ChrI", BK006936.2 = "ChrII", BK006937.2 = "ChrIII", BK006938.2 = "ChrIV", BK006939.2 = "ChrV", BK006940.2 = "ChrVI", BK006941.2 = "ChrVII", BK006934.2 = "ChrVIII", BK006942.2 = "ChrIX", BK006943.2 = "ChrX", BK006944.2 = "ChrXI", BK006945.2 = "ChrXII", BK006946.2 = "ChrXIII", BK006947.3 = "ChrXIV", BK006948.2 = "ChrXV", BK006949.2 = "ChrXVI")

# Find the de novo mutations
allfreq <- allfreq %>% mutate(altY55 = case_when(REF == alleleSK1 & ALT == alleleY55 ~ "SK1",
                                                 ALT == alleleSK1 & REF == alleleY55 ~ "SK1",
                                                 TRUE ~ "novo"))
# Some sites are the same in SK1 and Y55 and they will look nearly fixed or fixed at the founder,
# but that is because freqY55 then reflects the ancestral state between them, not Y55 allele only.

# There are a few things for sure at intermediate frequencies in the founder
allfreq %>% filter(Generation == "Founder", altY55 == "novo", freqY55 < 0.7) %>% .$Coordinate %>% unique %>% length() # 31
# For the file that has only "nice" biallelic sites I used a maf > 0.3 to filter. 

# How many de novo mutations (plus ancestral allele of SK1 and Y55) per environment?
allfreq %>% filter(altY55 == "novo", Replicate == "R1", Generation == "G30") %>% count(Environment)

# Extract the SNPs that had intermediate freqs in the Founder but don't look like the parentals
# I have to do it like this because the SNPs present in each environment are different

allfreq_denovo <- data.frame()
for (env in c("NaCl", "LiAc0.01", "LiAc0.02", "Ethanol")) {
  allfreq_env <- allfreq %>% filter(Environment == env, altY55 == "novo") 
  
  # Remove sites that were at intermediate frequencies in the founder
  listsnps <- allfreq_env %>% filter(Generation == "Founder", freqY55 < 0.7) %>% .$Coordinate
  allfreq_env_keep <- allfreq_env[!allfreq_env$Coordinate %in% listsnps,]
  
  # Remove sites that never reach at least a certain frequency at a given time point
  gooddenovo <- allfreq_env_keep %>% filter(freqY55 <= 1 - minfreqmutation) %>% .$Coordinate
  allfreq_env_keep <- allfreq_env_keep[allfreq_env_keep$Coordinate %in% gooddenovo,]
  
  allfreq_denovo <- rbind(allfreq_denovo, allfreq_env_keep)
  # print(allfreq_env_keep$Coordinate %>% unique %>% length)
}

###  ----- Some numbers -----
# How many de novo mutations are left in the whole experiment?
allfreq_denovo$Coordinate %>% unique %>% length # 310

# How many SNPs and INDELs?
filter(allfreq_denovo, class == "SNP")$Coordinate %>% unique %>% length # 183 SNPs
filter(allfreq_denovo, class == "INDEL")$Coordinate %>% unique %>% length # 127 INDELs

# How many de novo mutations per environment?
allfreq_denovo %>% filter(Replicate == "R1", Generation == "G30") %>% dplyr::count(Environment)
# What if we sum them, is it the same as the total?
allfreq_denovo %>% filter(Replicate == "R1", Generation == "G30") %>% dplyr::count(Environment) %>% .$n %>% sum # 392
# That is considerably more than the actual total, so some sites are repeated
# Sites that are repeated
allfreq_denovo %>% filter(Replicate == "R1", Generation == "G30") %>% dplyr::count(Coordinate) %>% filter(n > 1) %>% nrow()
# 54 sites appear more than once!
# ---------------------------

# ============================
### Find the mutations that rose to fixation
# ============================
deNovofixers <- data.frame()
for (env in c("NaCl", "LiAc0.01", "LiAc0.02", "Ethanol")) {
  allfreq_denovo_env <- allfreq_denovo %>% filter(Environment == env)
  listsnps <- c()
  for (repi in allfreq_denovo_env$Replicate %>% unique) {
    allfreq_denovo_env_rep <- allfreq_denovo_env %>% filter(Replicate == repi)
    
    # Anything that reached at least a min freq (1 - selectedfreq) at some point in the experiment
    listsnps <- c(listsnps, allfreq_denovo_env_rep %>% filter(freqY55 < selectedfreq) %>% .$Coordinate)
    
    # # Do they go to fixation at the last time point?
    # if ("G1000" %in% allfreq_denovo_env_rep$Generation) { # If I do it like this, the criteria is that the mutation stays high till the end of the experiment
    #   listsnps <- c(listsnps, allfreq_denovo_env_rep %>% filter(Generation == "G1000") %>% filter(freqY55 < 0.25) %>% .$Coordinate)
    # } else {
    #   listsnps <- c(listsnps, allfreq_denovo_env_rep %>% filter(Generation == "G700") %>% filter(freqY55 < 0.25) %>% .$Coordinate)
    # }
  }
  
  # listsnps <- allfreq_denovo_env %>% filter(Generation == "G700" | Generation == "G1000") %>% filter(freqY55 < 0.25) %>% .$Coordinate
  allfreq_env_keep <- allfreq_denovo_env[allfreq_denovo_env$Coordinate %in% listsnps,]
  
  # allfreq_env_keep[allfreq_env_keep$Coordinate %in% fixedsnps,]
  deNovofixers <- rbind(deNovofixers, allfreq_env_keep)
}

### Get the list of mutations for manual curation
# Add effects of mutations (SNPeff)
deNovofixers_snpeff_raw <- merge(deNovofixers, allsnpeff %>% select(-c(Contig, position)), by = c("Coordinate", "Environment", "class"))

#################### Make a table with the relevant info to report -----
infogenes <- deNovofixers_snpeff_raw %>% filter(freqY55 <= 1 - minfreqmutation) %>% # The frequency requirements helps me focus only on the relevant time points (but before I already had asked for mutations that reach more than 0.35)
  select(Coordinate, Contig, pos, gene_name, Gene_ID, type, putative_effect, Environment, Replicate, class) %>% 
  distinct(.keep_all = TRUE) %>% arrange(Contig, pos, Environment, Replicate)

infogenes$Coordinate %>% unique %>% length # 92 mutations that reach a frequency of >= (1 - selectedfreq) will be subject to manual curation

# This file was subject to manual curation to determine the "bad_mutations" list
write.table(infogenes, file = snakemake@output$infogenes, sep = "\t", row.names=FALSE, quote=FALSE)

# Make a table with the trajectories
write.table(deNovofixers_snpeff_raw, 
            file = snakemake@output$deNovofixers_trajectories,
            sep = "\t", row.names=FALSE, quote=FALSE)
####################

## After manual curation, some sites turned out to be bad SNP-calling
allfreq_denovo %>% filter(Coordinate %in% bad_mutations) %>% filter(class == "SNP") %>% .$Coordinate %>% unique %>% length # 20 bad SNPs
allfreq_denovo %>% filter(Coordinate %in% bad_mutations) %>% filter(class == "INDEL") %>% .$Coordinate %>% unique %>% length # 10 bad INDELs

##### So remove them!!!!
allfreq_denovo <- allfreq_denovo %>% filter(!Coordinate %in% bad_mutations)
deNovofixers <- deNovofixers %>% filter(!Coordinate %in% bad_mutations)
#####

# How many mutations >= (1 - selectedfreq) freq survived?
deNovofixers$Coordinate %>% unique %>% length # 62

# How many SNPs and INDELs?
filter(deNovofixers, class == "SNP")$Coordinate %>% unique %>% length # 56 SNPs
filter(deNovofixers, class == "INDEL")$Coordinate %>% unique %>% length # 6 INDELs

# Add effects to the fixers after removing bad mutations
allfreq_denovo_snpeff <- merge(allfreq_denovo, allsnpeff %>% select(-c(Contig, position)), by = c("Coordinate", "Environment", "class"))
deNovofixers_snpeff <- merge(deNovofixers, allsnpeff %>% select(-c(Contig, position)), by = c("Coordinate", "Environment", "class"))

# Trajectories of the good mutations
# Make a table with the trajectories
write.table(deNovofixers_snpeff, 
            file = snakemake@output$deNovofixers_trajectories_curated,
            sep = "\t", row.names=FALSE, quote=FALSE)

# ============================
cat(paste0("Plot interesting mutations", "\n"))
# ============================

plotgenes_experiment <- function(coolgenes, subtitle = "", deNovofixersdf = deNovofixers_snpeff){
  # mypalette <- get_palette(c("#FC4E07", "#E7B800","#00AFBB"), length(deNovofixersdf %>% filter(gene_name %in% coolgenes) %>% .$gene_name %>% unique()))
  mypalette <- ggpubr::get_palette(c("#FC4E07", "#ffdd55ff","#00AFBB"), length(deNovofixersdf %>% filter(gene_name %in% coolgenes) %>% .$gene_name %>% unique()))

  names(mypalette) <- deNovofixersdf %>% filter(gene_name %in% coolgenes) %>% .$gene_name %>% unique()
  mutcolors <- data.frame(gene_name = allfreq_denovo_snpeff$gene_name %>% unique) %>% mutate(col = case_when(
    gene_name %in% coolgenes ~ mypalette[gene_name],
    gene_name %in% deNovofixersdf$gene_name ~ "black", 
    TRUE ~ "gray"))
  
  # Make a list to just highlight the interesting mutations
  mutcolors_list <- mutcolors %>% filter(col != "gray") %>% .$col
  names(mutcolors_list) <- mutcolors %>% filter(col != "gray") %>% .$gene_name
  
  p <- ggplot(deNovofixersdf, aes(x = gen, y = 1 - freqY55, colour = gene_name, fill = Coordinate)) + 
    geom_line(alpha = 1, size = 0.4) +
    # geom_line(data = deNovofixersdf %>% filter(gene_name %in% coolgenes), aes(x = gen, y = 1 - freqY55, colour = gene_name, linetype = gene_name), alpha = 1, size = 0.8) +
    geom_line(data = deNovofixersdf %>% filter(gene_name %in% coolgenes), aes(x = gen, y = 1 - freqY55, colour = gene_name), alpha = 1, size = 0.8) +
    geom_point(data = deNovofixersdf %>% filter(gene_name %in% coolgenes), aes(x = gen, y = 1 - freqY55, colour = gene_name, shape = gene_name)) +
    geom_line(data = anti_join(allfreq_denovo_snpeff, deNovofixersdf), aes(x = gen, y = 1 - freqY55, colour = gene_name), alpha = 0.25) + # The grays
    facet_grid(Environment~Replicate) +
    xlab("Generation") + ylab('Frequency of the de novo mutation') +
    theme_bw() + 
    scale_x_continuous(breaks = allfreq_denovo_snpeff$gen %>% unique %>% sort,
                       labels = c(0, "", 60, "", 200, 300, 400, 500, 700, 1000)) +
    scale_shape_manual(name = "Gene", breaks = coolgenes, values = seq(1, length(coolgenes))) + # Having the same name and breaks leads to a combined legend
    scale_fill_manual(name = "Gene", breaks = coolgenes, values = mutcolors_list) + # Having the same name and breaks leads to a combined legend # I need fill or the points extend to all genes, not just coolgenes
    scale_color_manual(name = "Gene", breaks = coolgenes, values = mutcolors_list) + # Having the same name and breaks leads to a combined legend
    # # scale_linetype_manual(name = "Gene", breaks = coolgenes, values = seq(1, length(coolgenes))) + # Having the same name and breaks leads to a combined legend
    labs(title = paste0("De novo mutations that reach a frequency of at least ", minfreqmutation*100, "% at some point\n", subtitle)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), 
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(size = 0.25),
          panel.grid.major.y = element_line(size = 0.25),
          axis.title = element_text(size=14),
          strip.text = element_text(face = "bold", size = 12),
          # strip.background = element_rect(fill="gray94"),
          strip.background = element_blank(),
          # legend.position = "none")
          # legend.position="right")
          legend.position="bottom")
  return(p)
}

# Change the names so they match the main text of the manuscript
deNovofixers_snpeff$Environment[deNovofixers_snpeff$Environment == "LiAc0.01"] <- "LiAc 0.01M"
deNovofixers_snpeff$Environment[deNovofixers_snpeff$Environment == "LiAc0.02"] <- "LiAc 0.02M"
allfreq_denovo_snpeff$Environment[allfreq_denovo_snpeff$Environment == "LiAc0.01"] <- "LiAc 0.01M"
allfreq_denovo_snpeff$Environment[allfreq_denovo_snpeff$Environment == "LiAc0.02"] <- "LiAc 0.02M"

# BioRxiv
biorxiv <- plotgenes_experiment(c("ISW2", "CYC8", "CMK2", "CNB1", "SNF4", "SNF3", "WIP1"), "Genes with multiple independent mutations") +
  labs(title = NULL)
ggsave(plot = biorxiv, 
       filename = snakemake@output$plotgenes, 
       width = 8.5, height = 6.5)

# ============================
### What is the effect of the mutations?
# ============================
### Save the genes for https://string-db.org/
# deNovofixers_snpeff %>% filter(gene_name %in% c("ISW2", "CYC8", "CMK2", "CNB1", "SNF4", "SNF3", "WIP1")) %>% .$Coordinate %>% unique

write.table(deNovofixers_snpeff %>% filter(putative_effect != "MODIFIER") %>% select(gene_name) %>% unique %>% .$gene_name %>% sort, 
            file = snakemake@output$genes4stringdb,
            sep = "\t", row.names=FALSE, quote=FALSE, col.names = FALSE)

## I think the sites annotated as "MODIFIER" fall in non-coding regions
# How many sites are non-coding?
allfreq_denovo_snpeff %>% filter(putative_effect == "MODIFIER") %>% .$Coordinate %>% unique() %>% length() # 156

allfreq_denovo_snpeff %>% filter(putative_effect == "LOW") %>% .$Coordinate %>% unique() %>% length() # 20
allfreq_denovo_snpeff %>% filter(putative_effect == "MODERATE") %>% .$Coordinate %>% unique() %>% length() # 76
allfreq_denovo_snpeff %>% filter(putative_effect == "HIGH") %>% .$Coordinate %>% unique() %>% length() # 28

allfreq_denovo_snpeff %>% filter(type == "missense_variant") %>% .$Coordinate %>% unique() %>% length() # 68
allfreq_denovo_snpeff %>% filter(type == "synonymous_variant") %>% .$Coordinate %>% unique() %>% length() # 20
allfreq_denovo_snpeff %>% filter(type == "stop_gained") %>% .$Coordinate %>% unique() %>% length() # 18
allfreq_denovo_snpeff %>% filter(type == "start_lost") %>% .$Coordinate %>% unique() %>% length() # 1

# How many unique genes?
allfreq_denovo_snpeff %>% filter(putative_effect != "MODIFIER") %>% .$gene_name %>% unique %>% length() # 103

## How about the chosen mutations with higher frequency?
deNovofixers_snpeff %>% select(Coordinate, class) %>% count(Coordinate, class) %>% filter(class == "SNP") %>% nrow() #56
deNovofixers_snpeff %>% select(Coordinate, class) %>% count(Coordinate, class) %>% filter(class == "INDEL") %>% nrow() #6

deNovofixers_snpeff %>% .$Coordinate %>% unique() %>% length() # 62
deNovofixers_snpeff %>% filter(putative_effect == "MODIFIER") %>% .$Coordinate %>% unique() %>% length() # 5

deNovofixers_snpeff %>% filter(putative_effect == "LOW") %>% .$Coordinate %>% unique() %>% length() # 6
deNovofixers_snpeff %>% filter(putative_effect == "MODERATE") %>% .$Coordinate %>% unique() %>% length() # 35
deNovofixers_snpeff %>% filter(putative_effect == "HIGH") %>% .$Coordinate %>% unique() %>% length() # 16

deNovofixers_snpeff %>% filter(type == "missense_variant") %>% .$Coordinate %>% unique() %>% length() # 35
deNovofixers_snpeff %>% filter(type == "synonymous_variant") %>% .$Coordinate %>% unique() %>% length() # 6
deNovofixers_snpeff %>% filter(type == "frameshift_variant") %>% .$Coordinate %>% unique() %>% length() # 6
deNovofixers_snpeff %>% filter(type == "stop_gained") %>% .$Coordinate %>% unique() %>% length() # 9
deNovofixers_snpeff %>% filter(type == "start_lost") %>% .$Coordinate %>% unique() %>% length() # 1

## Plot of the mutation effects
deNovofixers_snpeff_effects <- allsnpeff %>% filter(Coordinate %in% deNovofixers_snpeff$Coordinate)

metamutations <- rbind(data.frame(allsnpeff %>% filter(Coordinate %in% allfreq_denovo_snpeff$Coordinate), dataset = "All"),
                       data.frame(deNovofixers_snpeff_effects, dataset = "Selected"))

# Fix plotting order
metamutations$putative_effect <- factor(metamutations$putative_effect, levels = c("MODIFIER", "LOW", "MODERATE", "HIGH"))

denovoeffectsp <- ggplot(metamutations, aes(putative_effect, fill = class)) + 
  geom_bar() + facet_grid(dataset ~ Environment, scales = "free_y") + 
  scale_fill_manual(values = c("dodgerblue3", "firebrick4")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="gray94")) +
  xlab("Putative effect of mutation")

ggsave(plot = denovoeffectsp, 
       filename = snakemake@output$denovoeffectsp,
       width = 23, height = 10, units = "cm")

## What genes are in more than one sample?
# gene_counts <- deNovofixers_snpeff %>% filter(freqY55 < selectedfreq, putative_effect != "MODIFIER") %>% mutate(Condition = paste0(Environment, "_", Replicate)) %>% select(c(Condition, gene_name)) %>% distinct(.keep_all = TRUE) %>% count(gene_name)
gene_counts <- deNovofixers_snpeff %>% filter(freqY55 < 1 - minfreqmutation, putative_effect != "MODIFIER") %>% mutate(Condition = paste0(Environment, "_", Replicate)) %>% select(c(Condition, gene_name)) %>% distinct(.keep_all = TRUE) %>% dplyr::count(gene_name)
gene_counts %>% filter(n > 1)

## What about independent mutations?
# mutation_counts <- deNovofixers_snpeff %>% filter(freqY55 < selectedfreq, putative_effect != "MODIFIER") %>% mutate(Condition = paste0(Environment, "_", Replicate)) %>% select(c(Condition, gene_name, Coordinate)) %>% distinct(.keep_all = TRUE) %>% dplyr::count(gene_name)
mutation_counts <- deNovofixers_snpeff %>% filter(freqY55 < 1 - minfreqmutation, putative_effect != "MODIFIER") %>% mutate(Condition = paste0(Environment, "_", Replicate)) %>% select(c(Condition, gene_name, Coordinate)) %>% distinct(.keep_all = TRUE) %>% dplyr::count(gene_name)
mutation_counts %>% filter(n > 1)

# Or should I count only the ones that are not linked? 
# Manually looking at their trajectories, it looks like the mutations in SUR2 and two of the 
# SNF3 mutations are linked. The corrected list would be (with selectedfreq):
# CMK2 2
# CNB1 2
# CYC8 2
# ISW2 8
# SNF3 2 # Two mutations were linked (and SUR2 was removed)
# SNF4 2
# WIP1 2

# The corrected list would be (with 1 - minfreqmutation). But since I am reporting the 0.35 manually curated ones, I think the one above is better.
# CMK2 2
# CNB1 2
# CYC8 2
# ISW2 8
# PBP1 2 # This one is a bit weird, because it is the exact same mutation that just appears in the last time point of Ethanol_R1 and it doesn't reach more than 0.35 so maybe discard?
# SNF3 2
# SNF4 2
# WIP1 2 # It's the exact same substitution in two independent LiAc0.01 replicates (R3 and R5) ...

# ## -- For testing and manual curation --
# genetrajectory <- function(snpdf, geneid){
#   p <- ggplot(snpdf %>% filter(gene_name == geneid), aes(x = gen, y = 1-freqY55, colour = Coordinate)) +
#     geom_line(alpha = 0.9) +
#     facet_grid(Environment~Replicate) +
#     xlab("Generation") + ylab('Frequency of the de novo mutation') +
#     theme_bw() + 
#     ylim(0,1) +
#     labs(title = paste("Frequency of variant", geneid)) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.position="bottom")
#   return(p)
# }
# 
# genetrajectory(allfreq_denovo_snpeff, "CBK1")
# allfreq_denovo_snpeff %>% filter(gene_name == "ISW2", freqY55 <= 1 - minfreqmutation)
# 
# snptrajectory <- function(snpdf, snpid){
#   p <- ggplot(snpdf %>% filter(Coordinate == snpid), aes(x = gen, y = 1 - freqY55)) +
#     geom_line(alpha = 0.9) +
#     facet_grid(Environment~Replicate) +
#     xlab("Generation") + ylab('Frequency of de novo mutation') +
#     theme_bw() + 
#     ylim(0,1) +
#     labs(title = paste("Frequency of variant", snpid)) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#   return(p)
# }
# 
# snptrajectory(allfreq_denovo, "BK006948.2_887373")
# allfreq_denovo %>% filter(Coordinate == "BK006934.2_275579", freqY55 <= 1 - minfreqmutation)
# deNovofixers_snpeff %>% filter(Coordinate == "BK006941.2_259037", freqY55 <= 1 - minfreqmutation)
# deNovofixers_snpeff %>% filter(Coordinate == "BK006948.2_886695", freqY55 <= 1 - minfreqmutation)
# deNovofixers_snpeff %>% filter(Coordinate == "BK006948.2_887373", freqY55 <= 1 - minfreqmutation)
# 
# ### ---
# 
# ## How many mutations are there per sample?
# 
# # How many de novo mutations (plus ancestral allele of SK1 and Y55) per environment?
# deNovofixers_snpeff %>% mutate(Condition = paste0(Environment, "_", Replicate)) %>% filter(freqY55 < selectedfreq) %>% select(Coordinate, Condition)  %>% distinct() %>% dplyr::count(Condition) #%>% arrange(gen)

# Some mutations fix to frequency of 1, in what chromosomes are those?
deNovofixers_snpeff %>% filter(freqY55 < 0.1) %>% select(gene_name, Contig) %>% distinct(.keep_all = TRUE)
# gene_name  Contig
# 1     STE12 ChrVIII
# 2      CYC8   ChrII
# 3      SUR2   ChrIV
# 4      NRG1   ChrIV
# 5      SKO1  ChrXIV
# 6      CBK1  ChrXIV
# 7      WHI2   ChrXV
# 8      ISW2   ChrXV
# 9     LDB19   ChrXV

