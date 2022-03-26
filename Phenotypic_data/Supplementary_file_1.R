####Phenotypic Analyses and Plots####
#By Ciaran Gilchrist
#Updated 2022-03-22
#Analyses includes ANOVAs, linear models and linear mixed effect models
#Plots by ggplot2
#R version 4.1.0 (2021-05-18)
#Packages used - car, dplyr, emmeans, ggplot2, gridExtra and MASS

####install prerequisite packages####
install.packages(c("car", "dplyr", "emmeans", "ggplot2", "gridExtra", "MASS"))

####Import Data####
setwd("C:/Target/Directory")
#Change to directory containing Supplementary Table 1.txt
#Plots will be saved to /Plots directory and statistics summaries to /Stats directory in this directory by default

OD <- read.delim("Supplementary Table 1.txt")
str(OD)
#Population indicates environment adapted in. EvoN = NaCl 0.75M, EvoL1 = LiAc 0.01M, EvoL2 = LiAc 0.02M, EvoE = Ethanol 8%
#Environment indicates environment tested in (SC + stress); NaCl = NaCl 0.75M, LiAc0.01 = LiAc 0.01M, LiAc0.02 = LiAc 0.02M, Ethanol = Ethanol 8%, SC = no added stress.
#Replicate is which replicate population is being tested
#Generation is the number of generations populations have been evolved in their home environment
#Yield is the difference in OD after 24h Yield compared to 0h (blank corrected)
#Group is a combination of population and replicate - used as a random factor in lme tests
#Relative_Yield is the relative yield (based on OD600) vs the founder population.
OD$Population <- as.factor(OD$Population)
OD$Environment <- as.factor(OD$Environment)
OD$Replicate <- as.factor(OD$Replicate)
OD$Group <- as.factor(OD$Group)
#convert to factor to allow for use in lm etc.
str(OD)

####Mean relative yields#####
#subset by population in home environment
EvoN <- subset(OD, Population == "EvoN" & Environment == "NaCl")
EvoN <- droplevels(EvoN)
str(EvoN)


EvoL1 <- subset(OD, Population == "EvoL1" & Environment == "LiAc0.01")
EvoL1 <- droplevels(EvoL1)
str(EvoL1)

EvoL2 <- subset(OD, Population == "EvoL2" & Environment == "LiAc0.02")
EvoL2 <- droplevels(EvoL2)
str(EvoL2)


EvoE <- subset(OD, Population == "EvoE" & Environment == "Ethanol")
EvoE <- droplevels(EvoE)
str(EvoE)

#mean at each generation measured, including founder
N0 <-mean(subset(EvoN, Generation == "0")$Yield)
N1 <-mean(subset(EvoN, Generation == "100")$Yield)
N3 <-mean(subset(EvoN, Generation == "300")$Yield)
N5 <-mean(subset(EvoN, Generation == "500")$Yield)
N7 <-mean(subset(EvoN, Generation == "700")$Yield)
N1K <-mean(subset(EvoN, Generation == "1000")$Yield)
#divide increase from founder at 100 generations by increase from founder at highest yield (1000 generations)
(N1-N0)/(N1K-N0)
#value*100 is % of overall increase within first 100 generations
#75.9%

L1_0 <-mean(subset(EvoL1, Generation == "0")$Yield)
L1_1 <-mean(subset(EvoL1, Generation == "100")$Yield)
L1_3 <-mean(subset(EvoL1, Generation == "300")$Yield)
L1_5 <-mean(subset(EvoL1, Generation == "500")$Yield)
L1_7 <-mean(subset(EvoL1, Generation == "700")$Yield)
(L1_1-L1_0)/(L1_7-L1_0)
#29.7%

L2_0 <-mean(subset(EvoL2, Generation == "0")$Yield)
L2_1 <-mean(subset(EvoL2, Generation == "100")$Yield)
L2_3 <-mean(subset(EvoL2, Generation == "300")$Yield)
L2_5 <-mean(subset(EvoL2, Generation == "500")$Yield)
L2_7 <-mean(subset(EvoL2, Generation == "700")$Yield)
(L2_1-L2_0)/(L2_5-L2_0)
#100 generation has the highest mean yield. Next are 500, 700 and 300 in that order. All significantly higher than founder

E_0 <-mean(subset(EvoE, Generation == "0")$Yield)
E_1 <-mean(subset(EvoE, Generation == "100")$Yield)
E_3 <-mean(subset(EvoE, Generation == "300")$Yield)
E_5 <-mean(subset(EvoE, Generation == "500")$Yield)
E_7 <-mean(subset(EvoE, Generation == "700")$Yield)
E_1K <-mean(subset(EvoE, Generation == "1000")$Yield)
#500 generations has highest yield. 100 generations > 1000 generations. Next highest are 700, 1000 and 300 in that order. 
(E_1-E_0)/(E_5-E_0)
#89.3% (compared to max OD at 500 gens)

#####Founder Populations#####
#####Calculate Founder Relative Fitness#####
#Calculate mean in SC then compare vs. founder OD in stressful environments
OD_Founder <- subset(OD, Generation == "0")
OD_Founder <- droplevels(OD_Founder)
str(OD_Founder)
OD_Founder$Environment <- as.factor(OD_Founder$Environment)
NF <- subset(OD_Founder, Population == "EvoN")
NF <- droplevels(NF)
LEF <- subset(OD_Founder, Population == "EvoL1"| Population == "EvoL2" | Population == "EvoE")
LEF <- droplevels(LEF)

N0S <-mean(subset(OD_Founder, Environment == "SC_N")$Yield)
LE0S <- mean(subset(OD_Founder, Environment == "SC_LE")$Yield)

NF$Fitness_Cost <- NF$Yield/N0S
LEF$Fitness_Cost <- LEF$Yield/LE0S


OD_Founder <- rbind(NF, LEF)
write.table(OD_Founder, file = "Founder_Fitness.txt", sep = "\t", row.names = F)

#####Founder Population Analyses#####
#How different were the four environments in fitness cost?
#Compare vs. unstressed environment in founder populations
#NaCl & LiAc 0.02M pure massive cost, LiAc 0.01M wee cost, Ethanol 8% nae change in fitness
library(MASS)
library(car)
library(emmeans)

OD_Founder <- read.delim("Founder_Fitness.txt")
str(OD_Founder)
OD_Founder$Environment <- as.factor(OD_Founder$Environment)

NF <- subset(OD_Founder, Population == "EvoN")
NF <- droplevels(NF)

LEF <- subset(OD_Founder, Population == "EvoL1"| Population == "EvoL2" | Population == "EvoE")
LEF <- droplevels(LEF)

#EvoN founder stats
fm1 <- lm(Yield ~ Environment, data = NF)
#make model
fligner.test(Yield ~ Environment, data = NF)
#test for equal variance. p<0.05
plot(fm1)
boxcox(Yield ~ Environment, data = NF)
#not normally distributed. boxcox plot suggests transformation ^1.75
NF$YieldTransformed <- NF$Yield^1.75
fm1 <- lm(YieldTransformed ~ Environment, data = NF)
fligner.test(YieldTransformed ~ Environment, data = NF)
#>0.05
plot(fm1)
boxcox(YieldTransformed ~ Environment, data = NF)
#Improved distribution

summary(fm1)
Anova(fm1)
#Type II ANOVA for model
mmNF1 <- emmeans(fm1, ~ Environment)
tukNF <- pairs(mmNF1, simple = "Environment")
#post-hoc test for comparisons of Relative_Yield between various generations
write.table(tukNF, file = "Stats/NFounder_stats.txt", sep = "\t", row.names = F)

#EvoLE founder stats
fm1 <- lm(Yield ~ Environment, data = LEF)
#make model
fligner.test(Yield ~ Environment, data = LEF)
#test for equal variance. p<0.05
plot(fm1)
boxcox(Yield ~ Environment, data = LEF)
#not normally distributed. boxcox plot suggests transformation ^0.75
LEF$YieldTransformed <- LEF$Yield^0.75
fm1 <- lm(YieldTransformed ~ Environment, data = LEF)
fligner.test(YieldTransformed ~ Environment, data = LEF)
#<0.05
plot(fm1)
boxcox(YieldTransformed ~ Environment, data = LEF)
#Improved distribution

summary(fm1)
Anova(fm1)
#Type II ANOVA for model
mmNF1 <- emmeans(fm1, ~ Environment)
tukNF <- pairs(mmNF1, simple = "Environment")
#post-hoc test for comparisons of Relative_Yield between various generations
write.table(tukNF, file = "Stats/LEFounder_stats.txt", sep = "\t", row.names = F)

#####Founder Population Plots#####
OD_Founder <- read.delim("Founder_Fitness.txt")
str(OD_Founder)
OD_Founder$Environment <- as.factor(OD_Founder$Environment)

library(ggplot2)

p <- ggplot(OD_Founder, aes(x = Environment, y = Yield)) + 
  #Set plot parameters
  geom_violin(aes(fill = Environment)) + geom_jitter() + 
  #Plot boxplot and points
  scale_fill_manual(values=c("#AC1016", "#800080", "#228B22", "#4169E1", "#D3D3D3", "#A9A9A9")) +
  #Add colours
  labs(x="Environment", y="Fitness of founder populations") + 
  #Add axes labels and plot title
  theme_bw() + theme(panel.border = element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"), axis.ticks.length = unit(0.25, "cm"))
  #modify aesthetics (remove grids, plain background, text size etc.)

png("Plots/FounderPlots.png", res = 500, units = "in", height = 10, width = 10)
p
dev.off()

OD_Founder <- subset(OD_Founder, Environment != "SC_LE" & Environment != "SC_N")
OD_Founder <- droplevels(OD_Founder)
str(OD_Founder)

p <- ggplot(OD_Founder, aes(x = Environment, y = Fitness_Cost)) + 
  #Set plot parameters
  geom_violin(aes(fill = Environment)) + geom_jitter(size = 0.75) + 
  #Plot boxplot and points
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  #Dashed line at 1 (where fitness is equal to SC)
  scale_fill_manual(values=c("#AC1016", "#800080", "#228B22", "#4169E1")) +
  #Add colours
  labs(x="Environment", y="Fitness Cost", title = "A) Relative fitness of founder populations") + 
  #Add axes labels and plot title
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=10, face="bold"),axis.text=element_text(size=10), axis.ticks.length = unit(0.5, "cm"), 
                     axis.title=element_text(size=12,face="bold"), legend.title=element_text(size=10), legend.text=element_text(size=10))

png("Plots/FounderPlots_Relative.png", res = 500, units = "in", height = 5, width = 5)
p
dev.off()

####Adaptation in home environment####
#Testing for differences in relative yield in home environment between replicate populations at different generations
#How similar are the independent replicate populations in terms of yield?
#Answer = a wee bit

library(MASS)
library(car)
library(emmeans)

#Subset populations and drop founder population
EvoN <- subset(OD, Population == "EvoN" & Environment == "NaCl" & Generation != "0")
EvoN <- droplevels(EvoN)
EvoN$Generation <- as.factor(EvoN$Generation)
EvoN$Replicate <- as.character(EvoN$Replicate)
str(EvoN)

NF <- subset(OD, Population == "EvoN" & Environment == "NaCl" & Generation == "0")
NF$Generation <- as.factor(NF$Generation)
NF$Replicate <- as.character(NF$Replicate)
NF <- droplevels(NF)
NF$Replicate[NF$Replicate=="R0"]<- "R1"
EvoN <- rbind(NF, EvoN)
NF$Replicate[NF$Replicate=="R1"]<- "R2"
EvoN <- rbind(NF, EvoN)
NF$Replicate[NF$Replicate=="R2"]<- "R3"
EvoN <- rbind(NF, EvoN)
NF$Replicate[NF$Replicate=="R3"]<- "R4"
EvoN <- rbind(NF, EvoN)

EvoL1 <- subset(OD, Population == "EvoL1" & Environment == "LiAc0.01" & Generation != "0")
EvoL1 <- droplevels(EvoL1)
EvoL1$Generation <- as.factor(EvoL1$Generation)
EvoL1$Replicate <- as.character(EvoL1$Replicate)
str(EvoL1)

L1F <- subset(OD, Population == "EvoL1" & Environment == "LiAc0.01" & Generation == "0")
L1F$Generation <- as.factor(L1F$Generation)
L1F$Replicate <- as.character(L1F$Replicate)
L1F <- droplevels(L1F)
L1F$Replicate[L1F$Replicate=="R0"]<- "R1"
EvoL1 <- rbind(L1F, EvoL1)
L1F$Replicate[L1F$Replicate=="R1"]<- "R2"
EvoL1 <- rbind(L1F, EvoL1)
L1F$Replicate[L1F$Replicate=="R2"]<- "R3"
EvoL1 <- rbind(L1F, EvoL1)
L1F$Replicate[L1F$Replicate=="R3"]<- "R4"
EvoL1 <- rbind(L1F, EvoL1)
L1F$Replicate[L1F$Replicate=="R4"]<- "R5"
EvoL1 <- rbind(L1F, EvoL1)


EvoL2 <- subset(OD, Population == "EvoL2" & Environment == "LiAc0.02" & Generation != "0")
EvoL2 <- droplevels(EvoL2)
EvoL2$Generation <- as.factor(EvoL2$Generation)
EvoL2$Replicate <- as.character(EvoL2$Replicate)
str(EvoL2)

L2F <- subset(OD, Population == "EvoL2" & Environment == "LiAc0.02" & Generation == "0")
L2F$Generation <- as.factor(L2F$Generation)
L2F$Replicate <- as.character(L2F$Replicate)
L2F <- droplevels(L2F)
L2F$Replicate[L2F$Replicate=="R0"]<- "R1"
EvoL2 <- rbind(L2F, EvoL2)
L2F$Replicate[L2F$Replicate=="R1"]<- "R2"
EvoL2 <- rbind(L2F, EvoL2)
L2F$Replicate[L2F$Replicate=="R2"]<- "R3"
EvoL2 <- rbind(L2F, EvoL2)
L2F$Replicate[L2F$Replicate=="R3"]<- "R4"
EvoL2 <- rbind(L2F, EvoL2)
L2F$Replicate[L2F$Replicate=="R4"]<- "R5"
EvoL2 <- rbind(L2F, EvoL2)

EvoE <- subset(OD, Population == "EvoE" & Environment == "Ethanol" & Generation != "0")
EvoE <- droplevels(EvoE)
EvoE$Generation <- as.factor(EvoE$Generation)
EvoE$Replicate <- as.character(EvoE$Replicate)
str(EvoE)

EF <- subset(OD, Population == "EvoE" & Environment == "Ethanol" & Generation == "0")
EF$Generation <- as.factor(EF$Generation)
EF$Replicate <- as.character(EF$Replicate)
EF <- droplevels(EF)
EF$Replicate[EF$Replicate=="R0"]<- "R1"
EvoE <- rbind(EF, EvoE)
EF$Replicate[EF$Replicate=="R1"]<- "R2"
EvoE <- rbind(EF, EvoE)
EF$Replicate[EF$Replicate=="R2"]<- "R3"
EvoE <- rbind(EF, EvoE)
EF$Replicate[EF$Replicate=="R3"]<- "R4"
EvoE <- rbind(EF, EvoE)

#####EvoN#####
#linear model for all generations & replicates of EvoN in its home environment NaCl 0.75M
fm1 <- lm(Relative_Yield ~ Generation*Replicate, data = EvoN)
#make model
fligner.test(Relative_Yield ~ Generation, data = EvoN)
fligner.test(Relative_Yield ~ Replicate, data = EvoN)
#test for equal variance. RY~Gen - p<0.05
plot(fm1)
boxcox(Relative_Yield ~ Generation*Replicate, data = EvoN)
#boxcox suggests ^3.5 transformation
EvoN$Transformed_Relative_Yield <- EvoN$Relative_Yield^3.5
boxcox(Transformed_Relative_Yield ~ Generation*Replicate, data = EvoN)
fm1 <- lm(Transformed_Relative_Yield ~ Generation*Replicate, data = EvoN)
#make model
fligner.test(Transformed_Relative_Yield ~ Generation, data = EvoN)
fligner.test(Transformed_Relative_Yield ~ Replicate, data = EvoN)
#test for equal variance. RY~Gen - p<0.05
plot(fm1)


fm1 <- lm(Transformed_Relative_Yield ~ Generation*Replicate, data = EvoN)
fm2 <- lm(Transformed_Relative_Yield ~ Generation + Replicate, data = EvoN)
anova(fm1,fm2)
#significant difference - keep interaction
fm3 <- lm(Transformed_Relative_Yield ~ Generation:Replicate + Replicate, data = EvoN)
fm4 <- lm(Transformed_Relative_Yield ~ Generation:Replicate + Generation, data = EvoN)
anova(fm1,fm3)
anova(fm1, fm4)
fm5 <- lm(Transformed_Relative_Yield ~ Generation:Replicate, data = EvoN)
anova(fm1, fm5)
#no significant difference - keep interaction only

Anova(fm5)
#Type II ANOVA on model
summary(fm5)

N.emm <- emmeans(fm5, ~ Replicate:Generation)
Nrep <- contrast(N.emm, "pairwise", simple = "Replicate", combine = TRUE, adjust = "tukey")
#Tukey HSD post-hoc test for pairwise comparisons of replicates at all generations
write.table(Nrep, file = "Stats/EvoN_rep_stats.txt", sep = "\t", row.names = F)
Nrep <- read.delim("Stats/EvoN_rep_stats.txt")
Nrep$Significance <- ifelse(Nrep$p.value <= 0.001, "***", ifelse(Nrep$p.value <= 0.01, "**", ifelse(Nrep$p.value <= 0.05, "*", "ns")))
write.table(Nrep, file = "Stats/EvoN_rep_stats.txt", sep = "\t", row.names = F)

Nrep <- contrast(N.emm, "pairwise", simple = "Generation", combine = TRUE, adjust = "tukey")
write.table(Nrep, file = "Stats/EvoN_repgen_stats.txt", sep = "\t", row.names = F)
#Write to file to refer to later
Nrep <- read.delim("Stats/EvoN_repgen_stats.txt")
Nrep$Significance <- ifelse(Nrep$p.value <= 0.001, "***", ifelse(Nrep$p.value <= 0.01, "**", ifelse(Nrep$p.value <= 0.05, "*", "ns")))
write.table(Nrep, file = "Stats/EvoN_repgen_stats.txt", sep = "\t", row.names = F)

#####EvoL1#####
#linear model for all generations & replicates of EvoL1 in its home environment LiAc 0.01M
fm1 <- lm(Relative_Yield ~ Generation*Replicate, data = EvoL1)
#make model
fligner.test(Relative_Yield ~ Generation, data = EvoL1)
fligner.test(Relative_Yield ~ Replicate, data = EvoL1)
#test for equal variance. p<0.05 for generation
plot(fm1)
boxcox(Relative_Yield ~ Generation*Replicate, data = EvoL1)
#not normally distributed. boxcox plot suggests ^4
EvoL1$Relative_YieldTransformed <- EvoL1$Relative_Yield^4
fm1 <- lm(Relative_YieldTransformed ~ Generation*Replicate, data = EvoL1)
fligner.test(Relative_YieldTransformed ~ Generation, data = EvoL1)
fligner.test(Relative_YieldTransformed ~ Replicate, data = EvoL1)
plot(fm1)
boxcox(Relative_YieldTransformed ~ Generation*Replicate, data = EvoL1)

fm1 <- lm(Relative_YieldTransformed ~ Generation*Replicate, data = EvoL1)
fm2 <- lm(Relative_YieldTransformed ~ Generation + Replicate, data = EvoL1)
anova(fm1,fm2)
#significant difference - keep interaction
fm3 <- lm(Relative_YieldTransformed ~ Generation:Replicate + Replicate, data = EvoL1)
fm4 <- lm(Relative_YieldTransformed ~ Generation:Replicate + Generation, data = EvoL1)
anova(fm1,fm3)
anova(fm1, fm4)
fm5 <- lm(Relative_YieldTransformed ~ Generation:Replicate, data = EvoL1)
anova(fm1, fm5)
#no significant difference - keep interaction only

Anova(fm5)
#Type II ANOVA on model
summary(fm5)
L1.emm <- emmeans(fm5, ~ Replicate:Generation)
L1rep <- contrast(L1.emm, "pairwise", simple = "Replicate", combine = TRUE , adjust = "tukey")
#Tukey HSD post-hoc test for pairwise comparisons of replicates at all generations
write.table(L1rep, file = "Stats/EvoL1_rep_stats.txt", sep = "\t", row.names = F)
L1rep <- read.delim("Stats/EvoL1_rep_stats.txt")
L1rep$Significance <- ifelse(L1rep$p.value <= 0.001, "***", ifelse(L1rep$p.value <= 0.01, "**", ifelse(L1rep$p.value <= 0.05, "*", "ns")))
write.table(L1rep, file = "Stats/EvoL1_rep_stats.txt", sep = "\t", row.names = F)

L1rep <- contrast(L1.emm, "pairwise", simple = "Generation", combine = TRUE , adjust = "tukey")
write.table(L1rep, file = "Stats/EvoL1_repgen_stats.txt", sep = "\t", row.names = F)
#write table to file for later reference
L1rep <- read.delim("Stats/EvoL1_repgen_stats.txt")
L1rep$Significance <- ifelse(L1rep$p.value <= 0.001, "***", ifelse(L1rep$p.value <= 0.01, "**", ifelse(L1rep$p.value <= 0.05, "*", "ns")))
write.table(L1rep, file = "Stats/EvoL1_repgen_stats.txt", sep = "\t", row.names = F)


#####EvoL2#####
#linear model for all generations & replicates of EvoL2 in its haem environment LiAc 0.02M
fm1 <- lm(Relative_Yield ~ Generation*Replicate, data = EvoL2)
#make model
fligner.test(Relative_Yield ~ Generation, data = EvoL2)
fligner.test(Relative_Yield ~ Replicate, data = EvoL2)
#test for equal variance. p<0.05 for generation
plot(fm1)
boxcox(Relative_Yield ~ Generation*Replicate, data = EvoL2)
#not normally distributed. boxcox plot suggests ^0.8
EvoL2$Relative_YieldTransformed <- EvoL2$Relative_Yield^0.8
fm1 <- lm(Relative_YieldTransformed ~ Generation*Replicate, data = EvoL2)
fligner.test(Relative_YieldTransformed ~ Generation, data = EvoL2)
fligner.test(Relative_YieldTransformed ~ Replicate, data = EvoL2)
plot(fm1)
boxcox(Relative_YieldTransformed ~ Generation*Replicate, data = EvoL2)

fm1 <- lm(Relative_YieldTransformed ~ Generation*Replicate, data = EvoL2)
fm2 <- lm(Relative_YieldTransformed ~ Generation + Replicate, data = EvoL2)
anova(fm1,fm2)
#significant difference - keep interaction
fm3 <- lm(Relative_YieldTransformed ~ Generation:Replicate + Replicate, data = EvoL2)
fm4 <- lm(Relative_YieldTransformed ~ Generation:Replicate + Generation, data = EvoL2)
anova(fm1,fm3)
anova(fm1, fm4)
fm5 <- lm(Relative_YieldTransformed ~ Generation:Replicate, data = EvoL2)
anova(fm1, fm5)
#no significant difference - keep interaction only

Anova(fm5)
#Type II ANOVA test of final model
summary(fm5)
L2.emm <- emmeans(fm5, ~ Replicate:Generation)
L2rep <- contrast(L2.emm, "pairwise", simple = "Replicate", combine = TRUE , adjust = "tukey")
#Tukey HSD post-hoc test for pairwise comparisons of replicates at all generations
write.table(L2rep, file = "Stats/EvoL2_rep_stats.txt", sep = "\t", row.names = F)
#Write tae table for later reference
L2rep <- read.delim("Stats/EvoL2_rep_stats.txt")
L2rep$Significance <- ifelse(L2rep$p.value <= 0.001, "***", ifelse(L2rep$p.value <= 0.01, "**", ifelse(L2rep$p.value <= 0.05, "*", "ns")))
write.table(L2rep, file = "Stats/EvoL2_rep_stats.txt", sep = "\t", row.names = F)

L2rep <- contrast(L2.emm, "pairwise", simple = "Generation", combine = TRUE , adjust = "tukey")
#Tukey HSD post-hoc test for pairwise comparisons of replicates at all generations
write.table(L2rep, file = "Stats/EvoL2_repgen_stats.txt", sep = "\t", row.names = F)
L2rep <- read.delim("Stats/EvoL2_repgen_stats.txt")
L2rep$Significance <- ifelse(L2rep$p.value <= 0.001, "***", ifelse(L2rep$p.value <= 0.01, "**", ifelse(L2rep$p.value <= 0.05, "*", "ns")))
write.table(L2rep, file = "Stats/EvoL2_repgen_stats.txt", sep = "\t", row.names = F)

#####EvoE#####
#linear model for all generations & replicates of EvoE in its haem environment Ethanol 8%
fm1 <- lm(Relative_Yield ~ Generation*Replicate, data = EvoE)
#make model
fligner.test(Relative_Yield ~ Generation, data = EvoE)
fligner.test(Relative_Yield ~ Replicate, data = EvoE)
#test for equal variance. p<0.05 for generation
plot(fm1)
boxcox(Relative_Yield ~ Generation*Replicate, data = EvoE)
#not normally distributed. boxcox plot suggests ^1.75
EvoE$Relative_YieldTransformed <- EvoE$Relative_Yield^1.75
fm1 <- lm(Relative_YieldTransformed ~ Generation*Replicate, data = EvoE)
fligner.test(Relative_YieldTransformed ~ Generation, data = EvoE)
fligner.test(Relative_YieldTransformed ~ Replicate, data = EvoE)
plot(fm1)
boxcox(Relative_YieldTransformed ~ Generation*Replicate, data = EvoE)

fm1 <- lm(Relative_YieldTransformed ~ Generation*Replicate, data = EvoE)
fm2 <- lm(Relative_YieldTransformed ~ Generation + Replicate, data = EvoE)
anova(fm1,fm2)
#significant difference - keep interaction
fm3 <- lm(Relative_YieldTransformed ~ Generation:Replicate + Replicate, data = EvoE)
fm4 <- lm(Relative_YieldTransformed ~ Generation:Replicate + Generation, data = EvoE)
anova(fm1,fm3)
anova(fm1, fm4)
fm5 <- lm(Relative_YieldTransformed ~ Generation:Replicate, data = EvoE)
anova(fm1, fm5)
#no significant difference - keep interaction only

Anova(fm5)
#Type II ANOVA test on model
summary(fm5)
E.emm <- emmeans(fm5, ~ Replicate:Generation)
Erep <- contrast(E.emm, "pairwise", simple = "Replicate", combine = TRUE , adjust = "tukey")
#Tukey HSD post-hoc test for pairwise comparisons of replicates at all generations
write.table(Erep, file = "Stats/EvoE_rep_stats.txt", sep = "\t", row.names = F)
#Write table fir future reference
Erep <- read.delim("Stats/EvoE_rep_stats.txt")
Erep$Significance <- ifelse(Erep$p.value <= 0.001, "***", ifelse(Erep$p.value <= 0.01, "**", ifelse(Erep$p.value <= 0.05, "*", "ns")))
write.table(Erep, file = "Stats/EvoE_rep_stats.txt", sep = "\t", row.names = F)

Erep <- contrast(E.emm, "pairwise", simple = "Generation", combine = TRUE , adjust = "tukey")
#Tukey HSD post-hoc test for pairwise comparisons of replicates at all generations
write.table(Erep, file = "Stats/EvoE_repgen_stats.txt", sep = "\t", row.names = F)
Erep <- read.delim("Stats/EvoE_repgen_stats.txt")
Erep$Significance <- ifelse(Erep$p.value <= 0.001, "***", ifelse(Erep$p.value <= 0.01, "**", ifelse(Erep$p.value <= 0.05, "*", "ns")))
write.table(Erep, file = "Stats/EvoE_repgen_stats.txt", sep = "\t", row.names = F)


#####Replicate Population Plots#####
#Means and CI for Replicates

EvoN <- subset(OD, Population == "EvoN" & Environment == "NaCl")
EvoL1 <- subset(OD, Population == "EvoL1" & Environment == "LiAc0.01")
EvoL2 <- subset(OD, Population == "EvoL2" & Environment == "LiAc0.02")
EvoE <- subset(OD, Population == "EvoE" & Environment == "Ethanol")
library(dplyr)

EvoN_by_replicate <- EvoN %>% group_by(Group, Population, Generation, Environment, Replicate) %>% summarise(mean_relative_yield = mean(Relative_Yield),  sd = sd(Relative_Yield), n = n_distinct(Relative_Yield))
EvoN_by_replicate$ci <- 1.96*(EvoN_by_replicate$sd/sqrt(EvoN_by_replicate$n))

EvoL1_by_replicate <- EvoL1 %>% group_by(Group, Population, Generation, Environment, Replicate) %>% summarise(mean_relative_yield = mean(Relative_Yield),  sd = sd(Relative_Yield), n = n_distinct(Relative_Yield))
EvoL1_by_replicate$ci <- 1.96*(EvoL1_by_replicate$sd/sqrt(EvoL1_by_replicate$n))

EvoL2_by_replicate <- EvoL2 %>% group_by(Group, Population, Generation, Environment, Replicate) %>% summarise(mean_relative_yield = mean(Relative_Yield),  sd = sd(Relative_Yield), n = n_distinct(Relative_Yield))
EvoL2_by_replicate$ci <- 1.96*(EvoL2_by_replicate$sd/sqrt(EvoL2_by_replicate$n))

EvoE_by_replicate <- EvoE %>% group_by(Group, Population, Generation, Environment, Replicate) %>% summarise(mean_relative_yield = mean(Relative_Yield),  sd = sd(Relative_Yield), n = n_distinct(Relative_Yield))
EvoE_by_replicate$ci <- 1.96*(EvoE_by_replicate$sd/sqrt(EvoE_by_replicate$n))

Mean <- rbind(EvoN_by_replicate, EvoL1_by_replicate, EvoL2_by_replicate, EvoE_by_replicate)

#Add founder for each replicate


Mean_evo <- subset(Mean_reps, Generation != 0)
Mean_founder <- subset(Mean_reps, Generation == 0)

Mean_founder$Replicate[Mean_founder$Replicate=="R0"]<- "R1"
Mean_reps <- rbind(Mean_founder, Mean_evo)

Mean_founder$Replicate[Mean_founder$Replicate=="R1"]<- "R2"
Mean_reps <- rbind(Mean_founder, Mean_reps)

Mean_founder$Replicate[Mean_founder$Replicate=="R2"]<- "R3"
Mean_reps <- rbind(Mean_founder, Mean_reps)

Mean_founder$Replicate[Mean_founder$Replicate=="R3"]<- "R4"
Mean_reps <- rbind(Mean_founder, Mean_reps)

Mean_founder$Replicate[Mean_founder$Replicate=="R4"]<- "R5"
Mean_founder <- subset(Mean_founder, Population == "EvoL1" | Population == "EvoL2")
Mean_reps <- rbind(Mean_founder, Mean_reps)

write.table(Mean_reps, file = "Stats/replicate_means.txt", sep = "\t", row.names = F)


Mean_reps <- read.delim("Stats/replicate_means.txt")

library(ggplot2)
library(gridExtra)

p1 <- ggplot(data = subset(Mean_reps, Population == "EvoN"), aes(x = Generation, y = mean_relative_yield, colour = Replicate))  + 
  #define main plot parameters
  geom_line(size=0.75)  + ylim(0.55, 2.75) + labs(x = "Generation", y = "Relative Fitness", title = "E) NaCl adapted populations in NaCl 0.75M") + 
  #plots lines
  geom_point(data = subset(Mean_reps, Population == "EvoN"), aes(x = Generation, y = mean_relative_yield, fill = Replicate), shape = 21, colour = "black", size = 3) +
  #Add mean to figure as point
  #geom_linerange(aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, colour = Replicate), size = 0.75) + 
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, fill = Replicate)) +
  #add 95% CI shading
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  #add line at founder's relative fitness
  scale_colour_manual(values=c("#87CEFA", "#5AA3CF", "#1E90FF", "#4169E1")) + scale_fill_manual(values=c("#87CEFA", "#5AA3CF", "#1E90FF", "#4169E1")) +
  #Choose colours and fill based on Replicate. Colour = lines & error lines, fill = points
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700, 1000)) +
  #make x axis only show points measured
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=16, face="bold"),axis.text=element_text(size=16), axis.ticks.length = unit(0.5, "cm"), 
                     axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=16), legend.text=element_text(size=16))
  #Theme to remove grids, generate white background, define text size and tick size

p2 <- ggplot(data = subset(Mean_reps, Population == "EvoE"), aes(x = Generation, y = mean_relative_yield, colour = Replicate))  + 
  geom_line(size=0.75)  + ylim(0.55, 2.75) + labs(x = "Generation", y = "Relative Fitness", title = "B) EtOH adapted populations in EtOH 8%") + 
  geom_point(data = subset(Mean_reps, Population == "EvoE"), aes(x = Generation, y = mean_relative_yield, fill = Replicate), shape = 21, colour = "black", size = 3) +
  #geom_linerange(aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, colour = Replicate), size = 0.75) +
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, fill = Replicate)) +
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700, 1000)) +
  scale_colour_manual(values=c("#CD5C5C", "#FF0000", "#AC1016","#67000C")) + scale_fill_manual(values=c("#CD5C5C", "#FF0000", "#AC1016","#67000C")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=16, face="bold"),axis.text=element_text(size=16), axis.ticks.length = unit(0.5, "cm"), 
                     axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=16), legend.text=element_text(size=16))

p3 <- ggplot(data = subset(Mean_reps, Population == "EvoL1"), aes(x = Generation, y = mean_relative_yield, colour = Replicate))  + 
  geom_line(size=0.75)  + ylim(0.55, 2.75) + labs(x = "Generation", y = "Relative Fitness", title = "C) LiAc 0.01M adapted populations in LiAc 0.01M") + 
  geom_point(data = subset(Mean_reps, Population == "EvoL1"), aes(x = Generation, y = mean_relative_yield, fill = Replicate), shape = 21, colour = "black", size = 3) +
  #geom_linerange(aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, colour = Replicate), size = 0.75) +
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, fill = Replicate)) +
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700)) +
  scale_colour_manual(values=c("#D8BFD8", "#DA70D6", "#9400D3", "#800080" , "#4B0082")) + scale_fill_manual(values=c("#D8BFD8", "#DA70D6", "#9400D3", "#800080" , "#4B0082")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=16, face="bold"),axis.text=element_text(size=16), axis.ticks.length = unit(0.5, "cm"), 
                     axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=16), legend.text=element_text(size=16))

p4 <- ggplot(data = subset(Mean_reps, Population == "EvoL2"), aes(x = Generation, y = mean_relative_yield, colour = Replicate))  + 
  geom_line(size=0.75)  + ylim(0.55, 2.75) + labs(x = "Generation", y = "Relative Fitness", title = "D) LiAc 0.02M adapted populations in LiAc 0.02M") + 
  geom_point(data = subset(Mean_reps, Population == "EvoL2"), aes(x = Generation, y = mean_relative_yield, fill = Replicate), shape = 21, colour = "black", size = 3) +
  #geom_linerange(aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, colour = Replicate), size = 0.75) +
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, fill = Replicate)) +
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700)) + 
  scale_colour_manual(values=c("#98D493", "#9ACD32", "#228B22", "#6B8E23", "#00441B")) + scale_fill_manual(values=c("#98D493", "#9ACD32", "#228B22", "#6B8E23", "#00441B")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=16, face="bold"),axis.text=element_text(size=16), axis.ticks.length = unit(0.5, "cm"), 
                     axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=16), legend.text=element_text(size=16))

png("Plots/Adaptation_ByReplicate_legend.png", res = 500, units = "in", height = 10, width = 15)
grid.arrange(p2, p3, p4, p1, nrow= 2)
dev.off()

OD_Founder <- read.delim("Founder_Fitness.txt")
str(OD_Founder)
OD_Founder$Environment <- as.factor(OD_Founder$Environment)

OD_Founder <- subset(OD_Founder, Environment != "SC_LE" & Environment != "SC_N")
OD_Founder <- droplevels(OD_Founder)
str(OD_Founder)

p6 <- ggplot(OD_Founder, aes(x = Environment, y = Fitness_Cost)) + 
  #Set plot parameters
  geom_violin(aes(fill = Environment)) + geom_jitter(size = 0.9) + 
  #Plot boxplot and points
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  #Dashed line at 1 (where fitness is equal to SC)
  scale_fill_manual(values=c("#AC1016", "#800080", "#228B22", "#4169E1")) +
  #Add colours
  labs(x="Environment", y="Fitness Cost", title = "A) Relative fitness of founder populations") + 
  #Add axes labels and plot title
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=16, face="bold"),axis.text=element_text(size=16), axis.ticks.length = unit(0.5, "cm"), 
                     axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=16), legend.text=element_text(size=16))

png("Plots/FounderPlots_Relative.png", res = 500, units = "in", height = 5, width = 15)
p6
dev.off()

#Final Figure 2 plot (combined p6 with p1-4 using Inkscape for manuscript, but can use below code to generate similar plot)

png("Plots/Figure 2.png", res = 750, units = "in", height = 15, width = 15)
grid.arrange(p6, p2, p3, p4, p1, ncol= 2)
dev.off()

#####SC Plots#####
#Means and CI for Replicates in unstressed conditions
library(dplyr)

OD_SC <- subset(OD, Environment == "SC")

SC_by_reppopgen <- OD_SC %>% group_by(Population, Generation, Group,Replicate) %>% summarise(mean_relative_yield = mean(Relative_Yield),  sd = sd(Relative_Yield), n = n_distinct(Relative_Yield))
SC_by_reppopgen$ci <- 1.96*(SC_by_reppopgen$sd/sqrt(SC_by_reppopgen$n))

Mean_evo <- subset(SC_by_reppopgen, Generation != 0)
Mean_founder <- subset(SC_by_reppopgen, Generation == 0)

Mean_founder$Replicate[Mean_founder$Replicate=="R0"]<- "R1"
Mean_reps <- rbind(Mean_founder, Mean_evo)

Mean_founder$Replicate[Mean_founder$Replicate=="R1"]<- "R2"
Mean_reps <- rbind(Mean_founder, Mean_reps)

Mean_founder$Replicate[Mean_founder$Replicate=="R2"]<- "R3"
Mean_reps <- rbind(Mean_founder, Mean_reps)

Mean_founder$Replicate[Mean_founder$Replicate=="R3"]<- "R4"
Mean_reps <- rbind(Mean_founder, Mean_reps)

Mean_founder$Replicate[Mean_founder$Replicate=="R4"]<- "R5"
Mean_founder <- subset(Mean_founder, Population == "EvoL1" | Population == "EvoL2")
Mean_reps <- rbind(Mean_founder, Mean_reps)

write.table(Mean_reps, file = "Stats/SC_means.txt", sep = "\t", row.names = F)

SC_by_popgen <- OD_SC %>% group_by(Population, Generation) %>% summarise(mean_relative_yield = mean(Relative_Yield),  sd = sd(Relative_Yield), n = n_distinct(Relative_Yield))
SC_by_popgen$ci <- 1.96*(SC_by_popgen$sd/sqrt(SC_by_popgen$n))
write.table(SC_by_popgen, file = "Stats/SC_byPop_means.txt", sep = "\t", row.names = F)

#plots
Mean_reps <- read.delim("Stats/SC_means.txt")
Mean_reps$Replicate[Mean_reps$Replicate=="R1"]<- "R1_anc"
Mean_reps$Replicate[Mean_reps$Replicate=="R2"]<- "R2_anc"
Mean_reps$Replicate[Mean_reps$Replicate=="R3"]<- "R3_anc"
Mean_reps$Replicate[Mean_reps$Replicate=="R4"]<- "R4_anc"
Mean_reps$Replicate[Mean_reps$Replicate=="R5"]<- "R5_anc"
Mean_reps <- subset(Mean_reps, select = -c(Group, n))

Mean_reps_adapted <- read.delim("Stats/replicate_means.txt")
Mean_reps_adapted$Replicate[Mean_reps_adapted$Replicate=="R1"]<- "R1_adap"
Mean_reps_adapted$Replicate[Mean_reps_adapted$Replicate=="R2"]<- "R2_adap"
Mean_reps_adapted$Replicate[Mean_reps_adapted$Replicate=="R3"]<- "R3_adap"
Mean_reps_adapted$Replicate[Mean_reps_adapted$Replicate=="R4"]<- "R4_adap"
Mean_reps_adapted$Replicate[Mean_reps_adapted$Replicate=="R5"]<- "R5_adap"
Mean_reps_adapted <- subset(Mean_reps_adapted, select = -c(Group, Environment))

Mean_reps <- rbind(Mean_reps, Mean_reps_adapted)
Mean_reps$Replicate <- ordered(Mean_reps$Replicate, levels = c("R1_anc", "R2_anc", "R3_anc", "R4_anc", "R5_anc", "R1_adap", "R2_adap", "R3_adap", "R4_adap", "R5_adap"))
str(Mean_reps)

library(ggplot2)
library(gridExtra)

p1 <- ggplot(data = subset(Mean_reps, Population == "EvoN"), aes(x = Generation, y = mean_relative_yield, colour = Replicate))  + 
  #define main plot parameters
  geom_line(size=0.75) + labs(x = "Generation", y = "Relative Fitness", title = "NaCl adapted populations in ancestral media") + 
  #plots lines
  geom_point(data = subset(Mean_reps, Population == "EvoN"), aes(x = Generation, y = mean_relative_yield, fill = Replicate), shape = 21, colour = "black", size = 3) +
  #Add mean to figure as point
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, fill = Replicate)) +
  #add 95% CI shading
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  #add line at founder's relative fitness
  scale_colour_manual(values=c("#87CEFA", "#5AA3CF", "#1E90FF", "#4169E1", "#D3D3D3", "#D3D3D3","#D3D3D3","#D3D3D3")) +
  scale_fill_manual(values=c("#87CEFA", "#5AA3CF", "#1E90FF", "#4169E1", "#D3D3D3", "#D3D3D3","#D3D3D3","#D3D3D3")) +
  #Choose colours and fill based on Replicate. Colour = lines & error lines, fill = points
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700, 1000)) +
  scale_y_continuous(breaks=seq(0.5, 3 , 0.5), limit = c(0.4, 3)) +
  #make x axis only show points measured
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=10, face="bold"),axis.text=element_text(size=14), axis.ticks.length = unit(0.25, "cm"), 
                     axis.title=element_text(size=14,face="bold"))
#Theme to remove grids, generate white background, define text size and tick size

p2 <- ggplot(data = subset(Mean_reps, Population == "EvoE"), aes(x = Generation, y = mean_relative_yield, colour = Replicate))  + 
  geom_line(size=0.75) + labs(x = "Generation", y = "Relative Fitness", title = "EtOH adapted populations in ancestral media") + 
  geom_point(data = subset(Mean_reps, Population == "EvoE"), aes(x = Generation, y = mean_relative_yield, fill = Replicate), shape = 21, colour = "black", size = 3) +
  #geom_linerange(aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, colour = Replicate), size = 0.75) +
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, fill = Replicate)) +
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700, 1000)) +
  scale_y_continuous(breaks=seq(0.5, 3 , 0.5), limit = c(0.4, 3)) +
  scale_colour_manual(values=c("#CD5C5C", "#FF0000",  "#AC1016", "#67000C", "#D3D3D3", "#D3D3D3","#D3D3D3","#D3D3D3")) +
  scale_fill_manual(values=c("#CD5C5C", "#FF0000",  "#AC1016", "#67000C", "#D3D3D3", "#D3D3D3","#D3D3D3","#D3D3D3")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=10,face="bold"), axis.text=element_text(size=14), axis.ticks.length = unit(0.25, "cm"), 
                     axis.title=element_text(size=14,face="bold"))

p3 <- ggplot(data = subset(Mean_reps, Population == "EvoL1"), aes(x = Generation, y = mean_relative_yield, colour = Replicate))  + 
  geom_line(size=0.75) + labs(x = "Generation", y = "Relative Fitness", title = "LiAc 0.01M adapted populations in ancestral media") + 
  geom_point(data = subset(Mean_reps, Population == "EvoL1"), aes(x = Generation, y = mean_relative_yield, fill = Replicate), shape = 21, colour = "black", size = 3) +
  #geom_linerange(aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, colour = Replicate), size = 0.75) +
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, fill = Replicate)) +
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700, 1000)) +
  scale_y_continuous(breaks=seq(0.5, 3 , 0.5), limit = c(0.4, 3)) +
  scale_colour_manual(values=c("#D8BFD8", "#DA70D6", "#9400D3", "#800080", "#4B0082", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3")) + 
  scale_fill_manual(values=c("#D8BFD8", "#DA70D6", "#9400D3", "#800080", "#4B0082", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=10,face="bold"), axis.text=element_text(size=14), axis.ticks.length = unit(0.25, "cm"),
                     axis.title=element_text(size=14,face="bold"))

p4 <- ggplot(data = subset(Mean_reps, Population == "EvoL2"), aes(x = Generation, y = mean_relative_yield, colour = Replicate))  + 
  geom_line(size=0.75) + labs(x = "Generation", y = "Relative Fitness", title = "LiAc 0.02M adapted populations in ancestral media") + 
  geom_point(data = subset(Mean_reps, Population == "EvoL2"), aes(x = Generation, y = mean_relative_yield, fill = Replicate), shape = 21, colour = "black", size = 3) +
  #geom_linerange(aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, colour = Replicate), size = 0.75) +
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci, fill = Replicate)) +
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700, 1000)) +
  scale_y_continuous(breaks=seq(0.5, 3 , 0.5), limit = c(0.4, 3)) +
  scale_colour_manual(values=c("#98D493", "#9ACD32", "#228B22", "#6B8E23", "#00441B", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3")) + 
  scale_fill_manual(values=c("#98D493", "#9ACD32", "#228B22", "#6B8E23", "#00441B", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=10,face="bold"),  axis.text=element_text(size=14), axis.ticks.length = unit(0.25, "cm"),
                     axis.title=element_text(size=14,face="bold"))

png("Plots/SC_ByReplicate_withAdapted.png", res = 500, units = "in", height = 10, width = 10)
grid.arrange(p1, p2, p3, p4, nrow= 2)
dev.off()


library(dplyr)

p5 <- ggplot(data = Mean_pops, aes(x = Generation, y = mean_relative_yield, colour = Population, fill = Population))  + 
  #define main plot parameters
  geom_line(size=0.75)  + ylim(0.9, 1.6) + labs(x = "Generation", y = "Relative Fitness", title = "All populations in SC") + 
  #plots lines
  geom_point(shape = 21, colour = "black", size = 3) +
  #Add mean to figure as point
  geom_ribbon(alpha=0.1, colour = NA, aes(ymin = mean_relative_yield-ci, ymax = mean_relative_yield+ci)) +
  #add 95% CI shading
  geom_hline(yintercept=1, linetype="dashed", size = 0.5) +
  #add line at founder's relative fitness
  scale_colour_manual(values=c("#AC1016", "#800080", "#228B22", "#4169E1")) + scale_fill_manual(values=c("#AC1016", "#800080", "#228B22", "#4169E1")) +
  #Choose colours and fill based on Replicate. Colour = lines & error lines, fill = points
  scale_x_continuous(breaks=c(0, 100, 300, 500, 700, 1000)) +
  #make x axis only show points measured
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=12, face="bold"),axis.text=element_text(size=14), axis.ticks.length = unit(0.25, "cm"), 
                     axis.title=element_text(size=14,face="bold"))
#Theme to remove grids, generate white background, define text size and tick size

png("Plots/SC_Combined_legend.png", res = 500, units = "in", height = 5, width = 5)
p5
dev.off()