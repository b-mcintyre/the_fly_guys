library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(rsq)
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

#Create levene statistic
lev_stat <- with(wing_table_clean, 
                 LeveneDeviates(y = wing_size_mm, group = Allele_1:WT_Background, med = TRUE))

wing_table_lev_raw <- wing_table_clean %>% mutate ( lev_stat = lev_stat)

str(wing_table_clean)

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


#getting line means/sd/var and group sd/bx means/sd/var
#maybe plot the variance over the means? 

line_means_sd_var_cv <- wing_table_lev_raw %>% 
  group_by(Allele_1, WT_Background) %>% 
  summarise(length_means = mean(wing_size_mm), 
            length_sd = sd(wing_size_mm),
            length_cv = (sd(wing_size_mm)/mean(wing_size_mm)),
            lev_stat = median(lev_stat),
            Individuals = n()
            )

sd_means_sd_var_cv <- sddat %>% 
  group_by(Allele_1, WT_Background) %>%
  summarise(length_means = mean(wing_size_mm), 
            length_sd = sd(wing_size_mm),
            length_cv = (sd(wing_size_mm)/mean(wing_size_mm)),
            lev_stat = median(lev_stat),
            Individuals = n()
  )

bx_means_sd_var_cv <- bxdat %>%
  group_by(Allele_1, WT_Background) %>%
  summarise(length_means = mean(wing_size_mm), 
            length_sd = sd(wing_size_mm),
            length_cv = (sd(wing_size_mm)/mean(wing_size_mm)),
            lev_stat = median(lev_stat),
            Individuals = n()
  )


ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_sd, color=Allele_1, size=Individuals)) +geom_point()
  # kind of looks like a quadratic relationship? 

ggplot(sd_means_sd_var_cv, aes(x=length_means, y=length_sd, color=Allele_1, size=Individuals)) + geom_point()
# one maybe two outliers in the sd29.1 data set and wild type data set? 
# maybe not log transformation is not applicable? looks like there may be a pattern there however,
# would another transformation work better? quadratic? 

ggplot(bx_means_sd_var_cv, aes(x=length_means, y=length_sd, colour=Allele_1, size=Individuals)) + geom_point()
# this looks to have no observable pattern...

ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_sd, color=WT_Background)) + 
  geom_point() + 
  facet_wrap(~Allele_1)
# maybe not log transform? 

#display sample size for the dots 
#bx log trasform would be better --> using levenes in raw form would be fine 


ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_sd)) + 
  geom_point(aes(color=Allele_1, size=Individuals)) +
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) + 
  geom_smooth(method = lm, formula = y ~ poly(x, 3), color="red")

ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_cv)) +
  geom_point(aes(color=Allele_1, size=Individuals)) + 
  geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
  geom_smooth(method = lm, formula = y ~ poly(x, 3), color = "red")

#RXN norm plots for means

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=length_means)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

ggplot(sd_means_sd_var_cv, aes(x=Allele_1, y=length_means)) +
  geom_line(aes(color=WT_Background, group=WT_Background)) +
  geom_point(aes(color=WT_Background))


ggplot(bx_means_sd_var_cv, aes(x=Allele_1, y=length_means)) +
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

#RXN norm plots for CV

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=length_cv)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

ggplot(sd_means_sd_var_cv, aes(x=Allele_1, y=length_cv)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

ggplot(bx_means_sd_var_cv, aes(x=Allele_1, y=length_cv)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

## RXN norm plots for median levene's statistic

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=lev_stat)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))


ggplot(sd_means_sd_var_cv, aes(x=Allele_1, y=lev_stat)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

ggplot(bx_means_sd_var_cv, aes(x=Allele_1, y=lev_stat)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

#### testing ####

#Two step method, using the polynomial 2 model fit to the previous data using sd

lmmodel1 <- lm(length_sd ~ poly(length_means, 2), data = line_means_sd_var_cv)

lmmodel2 <- lm(length_sd ~ poly(length_means, 3), data = line_means_sd_var_cv)

anova(lmmodel1, lmmodel2)

#this glm model won't fit at the second step/level --> why?
#use log link because an gamma stuggests that no parameters can be negatie
#identity says that yes some parameters can be negative
#therefore use a log link as it also says no parameters can be negative 

glmmodel1 <- glm(lm(length_sd ~ poly(length_means, 2),
                    family = Gamma(link = log),
                    start = coef(lmmodel1),
                    data = sd_means_sd_var_cv))
#model also shows singularity, uses sd as the comparison statistic.

summary(glmmodel1)

Anova(glmmodel1)

#### rough among line SD/CV, compare from model later ####

rough_among_line_sd1 <- line_means_sd_var_cv %>% group_by(Allele_1) %>%
  summarize(line_sd = sd(length_means))

rough_among_genotype_sd1 <- line_means_sd_var_cv %>% group_by(Allele_1) %>%
  summarise(genotype_sd = mean(length_sd))


rough_among_line_mean1 <- line_means_sd_var_cv %>% group_by(Allele_1) %>% 
  summarise(line_mean = mean(length_means))

rough_among_line_cv1 <- line_means_sd_var_cv %>% group_by(Allele_1) %>% 
  summarise(line_cv = sd(length_means)/length_means)



#### multilevel model ####
factor(wing_table_clean$Replicate)

#JD code suggestion results in a warning message of checking convergence 
JDall_glm_wing_size_lev <- lmer(lev_stat ~  1 + Allele_1 + Replicate + (0 + Allele_1 | WT_Background),
                       data = wing_table_lev_raw)

#Original coding:
all_glm_wing_size_lev <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
                              + (1|Replicate),
                              data = wing_table_clean)

# this is where JD said replicate should be treated as a fixed effect because there is only two levels(3ish)
# which he believes is causing the singularity 
#BB suggests using (1|background/allele) making the assumption that each of the alleles covary the same
#this seems like a very strong assumption.
#BB wants to know at what level we are looking into/what the goal of this whole analysis is

all_glm_wing_size_size <- lmer(wing_size_mm ~ 1 + Allele_1 + (0 + Allele_1 | WT_Background) 
                               + (1 | Replicate),
                               data = wing_table_clean)

#over paramaterized for both models? 

summary(all_glm_wing_size_lev)

summary(all_glm_wing_size_size)

Anova(all_glm_wing_size_lev)

Anova(all_glm_wing_size_size)

#R2 values
rsq(all_glm_wing_size_lev)
#model doesn't explain a lot of the variability (less than half? a bunch of noise?)

rsq(all_glm_wing_size_size)

#model explains a lot of the size effects