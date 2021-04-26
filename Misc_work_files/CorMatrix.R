library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(rsq)
library(emmeans)
library(MCMCglmm)
library(dhglm)
library(corrplot)
#### Cleaning Data ####


# put into rda file once sure is completely clean 

raw_wing_table <- read_csv("NEW_CD_DGRP_Subset_Data_2019_V2.csv")

str(raw_wing_table)

wing_table_clean <- raw_wing_table %>% 
  mutate(wing_size_mm = sqrt(TotalArea.px * 0.00215^2),
         TotalArea.px = NULL,
         Temperature = NULL, 
         Allele_1 = factor(Allele_1, levels = c("OREw", "bx[1]", "bx[2]", 
                                                "bx[3]","sd[29.1]", "sd[1]", 
                                                "sd[E3]", "sd[ETX4]", "sd[58d]")), 
         WT_Background = factor(WT_Background))

levels(wing_table_clean$WT_Background)
levels(wing_table_clean$Allele_1)

source("scripts/ID_LeveneStat_V1_2016.R")


wing_table_lev_raw <- wing_table_clean %>%
  mutate ( lev_stat = LeveneDeviates(y = wing_size_mm, 
                                     group = Allele_1:WT_Background, med = TRUE))

str(wing_table_clean)

saveRDS(wing_table_lev_raw, file ="wing_table_lev_raw.rds")

bxdat <- wing_table_lev_raw %>% 
  filter(Allele_1 %in% c("OREw", "bx[1]", "bx[2]", "bx[3]")) %>%
  droplevels()

str(bxdat)

levels(bxdat$Allele_1)

sddat<- wing_table_lev_raw %>% 
  filter(Allele_1 %in% c("OREw", "sd[1]", "sd[29.1]", "sd[58d]", "sd[E3]", "sd[ETX4]")) %>%
  droplevels() 

str(sddat)

levels(sddat$Allele_1)

prior <- list( R = list(V=diag(9)/9, nu=0.004),  
               G = list(G1=list(V=diag(9)/9, nu=0.004)))

MCMCglmmmTotal1 <- MCMCglmm(fixed = lev_stat ~ 1 + Allele_1,
                            random =~ us(Allele_1):WT_Background,
                            rcov = ~idh(Allele_1):units,
                            data = wing_table_lev_raw, prior=prior,
                            nitt = 20000, burnin = 5000, thin = 10)
summary(MCMCglmmmTotal1)
summary(MCMCglmmmTotal1$Sol) # fixed effects
summary(MCMCglmmmTotal1$VCV)
a <- summary(MCMCglmmmTotal1$VCV)$statistics[,"Mean"]; a
length(a)# extracting posterior means from the variances and covariances
a <- a[1:81];a 
G_mat1 <- matrix(a, nrow = 9, ncol = 9, byrow = T); G_mat1

G_cor1 <- cov2cor(G_mat1); G_cor1 # genetic variance covariance matrix for the (co)variation among DGRP across mutant alleles

#visualization
colnames(G_cor1) <- c("Wild Type", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")
rownames(G_cor1) <- c("Wild Type", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")

col1 <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(G_cor1, type = "lower", method = "color", col=col1(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)

