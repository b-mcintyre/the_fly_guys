library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(car)
library(rsq)
library(emmeans)
library(corrplot)
library(corpcor)
library(jtools)
library(glmmTMB)

counts <- wing_table_clean %>% group_by(Allele_1, WT_Background) %>% summarise(n = n())


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


sddat <- sddat %>% mutate(Replicate = factor(Replicate))


bxdat <- bxdat %>% mutate(Replicate = factor(Replicate))


m1 <- lmer(wing_size_mm ~ 1 + (1|WT_Background), data = wing_table_clean)
summary(m1)
#not a huge amount of differece in background means

summ(m1)

m2 <-lmer(wing_size_mm ~ 1 + (1|Allele_1), data = wing_table_clean)
summary(m2)
#decent amount of variation in allele means
summ(m2)

m3 <- lmer(wing_size_mm ~ 1 + Allele_1 + (1|WT_Background), data = wing_table_clean)
summary(m3)

summ(m3)

plot(m3)

m4 <- lmer(wing_size_mm ~ 1 + Allele_1 + (1|WT_Background) + (1|Replicate), data = wing_table_clean)

summary(m4)
#all very close to the wild type background

c1 <- cov2cor(VarCorr(m3)[[1]])

corrplot(c1, type = "lower", method = "color", col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)

model1 <- lm(wing_size_mm ~ 1 + Allele_1 + WT_Background + Replicate + Allele_1:WT_Background, 
             data = wing_table_clean)

summary(model1)
summ(model1)
#lota sort through..
Anova(model1)


model2 <- lmer(wing_size_mm ~ 1 + Allele_1 + (1|WT_Background:Allele_1) + (1|Replicate), data = wing_table_clean)

summary(model2)
#need to think about this model and figure out exactly what interactions in glms mean.



model3 <- lmer(wing_size_mm ~ 1 +Allele_1 + (0 + Allele_1||WT_Background) + (1|Replicate), data = wing_table_clean)

summary(model3)

model4 <- lmer(wing_size_mm ~ 1 + Allele_1 + (0 + Allele_1||WT_Background) + Replicate, data = wing_table_clean)
#both still singular fits...


sdModel <- lmer(win_size_mm ~ 1 + Allele_1 + (0 + Allele_1|W_Background) + Replicate, data= )
line_means_sd_var_cv <- wing_table_lev_raw %>% 
  group_by(Allele_1, WT_Background) %>% 
  summarise(length_means = mean(wing_size_mm), 
            length_sd = sd(wing_size_mm),
            length_cv = (sd(wing_size_mm)/mean(wing_size_mm)),
            lev_stat = mean(lev_stat),
            Individuals = n()
  )

# could try centering?

grand_mean <- mean(wing_table_clean$wing_size_mm)

cline_means_sd_var_cv <- line_means_sd_var_cv %>% mutate(clength_means = length_means - grand_mean)
cline_means_sd_var_cv

cwing_cleaned <- wing_table_clean %>% mutate(cwing_size_mm = wing_size_mm - grand_mean)

cwing_cleaned

model1 <- lmer(cwing_size_mm ~1 + Allele_1 + (0 + Allele_1|WT_Background) + (1|Replicate), data = cwing_cleaned)
#still singular with centered values 




#### WHAT I THINK I ACTUALLY NEED! ####
#random effects anova in R

rand.mutanteffects.anova = lmer(wing_size_mm ~ 1 + (1|Allele_1), data = wing_table_clean)

summary(rand.mutanteffects.anova)

ggplot(wing_table_clean, aes(x=Allele_1, y=wing_size_mm)) + geom_violin
(position = "jitter") + geom_hline(yintercept = 1.01658)

install.packages("devtools")
library(devtools)
devtools::install_github("dustinfife/fifer")
devtools::install_github("dustinfife/flexplot")

# install.packages("devtools")
# install the stable version
devtools::install_github("dustinfife/flexplot")
# install the development version
devtools::install_github("dustinfife/flexplot", ref="development")
library(flexplot)

flexplot(wing_size_mm~Allele_1, data=wing_table_clean, jitter = T, spread="stdev" ) + geom_hline(yintercept=1.01658)
data("relationship_satisfaction")
relationship_satisfaction

icc(rand.mutanteffects.anova)
# huge amount of the variability is due to the clusters around the mutants...
# which duh

#do same for the wing area stuff
rand.DGRPeffects.anova <- lmer(wing_size_mm ~ 1 + (1|WT_Background), data = wing_table_clean)
summary(rand.DGRPeffects.anova)
icc(rand.DGRPeffects.anova)
#3% variability due to the clustering around DGRPs so it is present...

rand.replicateeffects.anova <- lmer(wing_size_mm ~ 1 + (1|Replicate), data = wing_table_clean)
summary(rand.replicateeffects.anova)
icc(rand.replicateeffcts.anova)
# says large amount of the varibity is due to clusering which again duh

#now considerations on both... need to think about this

rand.intercept <- lmer(wing_size_mm ~ 1 + Allele_1 + (1|WT_Background), data = wing_table_clean)

summary(rand.intercept)


rand.slope <- lmer(wing_size_mm ~ 1 + Allele_1 + (-1 + Allele_1|WT_Background), data = wing_table_clean)

isSingular(rand.slope)

summary(rand.slope)


wtmodel <-lmer(wing_size_mm ~ 1 + WT_Background + (0 + WT_Background|Allele_1) + (1|Replicate), data = wing_table_clean)
