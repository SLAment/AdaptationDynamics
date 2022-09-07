#!/usr/bin/env Rscript

# ENAdepthCNV: Are there indications of copy number variations in the ENA gene?
#############################################################################
# Part of the Snakemake pipeline ENAdepthCNV.smk
# =======================================
# S. Lorena Ament-Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022/08/31
#############################################################################
# Load the necessary libraries
# ============================
library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff

# ============================
# Check input
# ============================
enafile <- snakemake@input$table
outputplot <- snakemake@output$plot

# ============================
# Ratio of ENA and flanks
# ============================
enadf <- read.table(enafile, header=TRUE) %>% mutate(ratio = covENA/covflanks)

RelENAp <- ggplot(enadf %>% filter(!Sample %in% c("SK1", "Y55", "S288c")), aes(ratio)) + geom_histogram() + xlim(0, 4) + 
  theme_bw() +
  theme(legend.position="none") +
  xlab(expression(paste("Relative depth of coverage of ", italic("ENA")))) +
  ylab("Frequency") +
  geom_point(data = enadf %>% filter(Sample %in% c("SK1", "Y55", "S288c")), aes(x = ratio, y = 0, colour = Sample)) +
  annotate("text", x=filter(enadf, Sample == "SK1") %>% .$ratio, y=5, label= "SK1", colour = "#00BA38") +
  annotate("text", x=filter(enadf, Sample == "Y55") %>% .$ratio, y=10, label= "Y55", colour = "#619CFF") +
  annotate("text", x=filter(enadf, Sample == "S288c") %>% .$ratio, y=5, label= "S288c", colour = "#F8766D")

ggsave(plot = RelENAp, filename = outputplot, width = 4, height = 3)
