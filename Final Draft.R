library(plyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(corrplot)
library(MASS)
library(blme)
library(rstanarm)
library(broom.mixed)
library(dotwhisker)
library(Matrix)
library(cowplot)
library(emmeans)
library(effects)
remotes::install_github("glmmTMB/glmmTMB/glmmTMB#670")
library(glmmTMB)
library(lattice)
library(coda)









#### Cleaning Data ####

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

wd <- readRDS("wing_table_lev_raw.rds")

#getting line means/sd/var and group sd/bx means/sd/var
#maybe plot the variance over the means? 


line_means_sd_var_cv <- wing_table_lev_raw %>% 
  group_by(Allele_1, WT_Background) %>% 
  summarise(length_means = mean(wing_size_mm), 
            length_sd = sd(wing_size_mm),
            length_cv = (sd(wing_size_mm)/mean(wing_size_mm)),
            lev_stat = mean(lev_stat),
            Individuals = n()
            )

#### visualizing data ####

ggplot(line_means_sd_var_cv, 
       aes(x=length_means, y=length_cv, color=Allele_1, size=Individuals)) +
  geom_point() + 
  labs(x = "Line means (mm)", 
       y = "Line Standard Deviation", 
       size = "Fly Sample Size",
       color = "Mutant Allele") +
theme(axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(size= 12))


ggplot(line_means_sd_var_cv, 
       aes(x=length_means, y=length_sd, color=Allele_1, size=Individuals)) + 
  geom_point() +
  labs(x = "Line Means (mm)", 
       y = "Line Coefficient of Variation", 
       size = "Fly Sample Size", 
       color = "Mutant Allele") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size= 12))


ggplot(line_means_sd_var_cv, 
       aes(x=length_means, y=length_cv, color=WT_Background)) + 
  geom_point() + 
  facet_wrap(~Allele_1) +
  labs(x = "Line Means (mm)", 
       y = "Line coefficent of Variation", 
       color = "DGRP Background",
       size = "Fly Sample Size") +
  theme(axis.text.x = element_text(angle = 90, size = 10), 
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size= 12))

  
# seeing by each mutant allele 

## RXN norm plot for raw mean levene's statistic

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=lev_stat)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background)) + 
  labs(y="Levene's Statistic", x= "Mutant Allele") + 
  theme(legend.position = "")

#boxplots for sd1, SdE3, and sd58d or wild type, moderate and severe mutational effects 

boxdat <- wd %>% filter(Allele_1 == c("sd[1]", "sd[E3]", "sd[58d]"))

ggplot(boxdat, aes(x=WT_Background, y=wing_size_mm)) + 
  facet_wrap(~Allele_1) + 
  geom_boxplot() + 
  theme(axis.text.x=element_text(angle=90, size = 9)) + 
  labs(x= "DGRP Line", y="Wing size mm")

# see more variation in the moderate allele than other alleles 

#### multilevel modeling ####

#JD code suggestion results in a warning message of checking convergence 
JDall_glm_wing_size_lev <- lmer(lev_stat ~  1 + Allele_1 + Replicate + (0 + Allele_1 | WT_Background),
                       data = wing_table_lev_raw)

allFit(JDall_glm_wing_size_lev)
#majority still singular 

#try nesting within the mutants 
JD2all_glm_wing_size_lev <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
                              + (1|Allele_1/Replicate),
                              data = wing_table_lev_raw)
#still singular

#try nesting within backgrounds
JD2all_glm_wing_size_lev1 <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
                                + (1|WT_Background/Replicate),
                                 data = wing_table_lev_raw)
#still singular 


#Original coding:
all_glm_wing_size_lev <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
                              + (1|Replicate),
                              data = wing_table_lev_raw)


# this is where JD said replicate should be treated as a fixed effect because there is only two levels(3ish)
# which he believes is causing the singularity/problems 

#### BB's rank reduced model with 5 components ####


m5 <- glmmTMB(lev_stat ~ Allele_1 + rr(0 + Allele_1 | WT_Background,5) +
                Replicate,
              data=wd,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))

summary(m5)
summary(wd)

plot(emtrends(m5))


# correlation matrix
v5 <- cov2cor(VarCorr(m5)$cond[[1]])

rr_mat <- matrix(v5, nrow = 9, ncol = 9, byrow = T); rr_mat

rr_cor <- cov2cor(rr_mat); rr_cor

#visualization of covariance matrix
colnames(rr_cor) <- c("Wild Type", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")
rownames(rr_cor)<- c("Wild Type", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))


corrplot(rr_cor, type = "lower", method = "color", col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)
#### plotting variance ####

#Anova
Anova(m5, type = "III")

## emmeans plot
plot(emmeans(m5, "Allele_1", "Replicate"),
     ylab = "Mutant Allele",
     xlab = "Within line variability",
     comparison = T,
     horizontal = F)

#Rxn norm of the mean plot

summary(m5)

ranef(m5)
coef(m5)


estimates <- coef(m5)$cond
estimates2 <- data.frame(estimates[1])

estimates_df <-  data.frame(DGRP =rownames(estimates2),
                          WT =  (estimates2[,2] + estimates2[,1]),
                          bx1 = (estimates2[,2] + estimates2[,3]),
                          bx2 = (estimates2[,2] + estimates2[,4]),
                          bx3 = (estimates2[,2] + estimates2[,5]),
                          sd29.1 = (estimates2[,2]  + estimates2[,6]), 
                          sd1 = (estimates2[,2] + estimates2[,7]), 
                          sdE3 = (estimates2[,2]  + estimates2[,8]),
                          sdETX4 = (estimates2[,2]  + estimates2[,9]),
                          sd58d = (estimates2[,2]  + estimates2[,10]))


dat_for_rxnnorm <- estimates_df %>% 
  pivot_longer(c(WT, bx1, bx2, bx3, sd29.1, sd1, sdE3, sdETX4, sd58d), 
               names_to = "Genotype", values_to = "lev_stat")


dat_for_rxnnorm <-dat_for_rxnnorm %>% 
  mutate(Genotype = factor(Genotype, 
                           levels = c("WT", "bx1", "bx2","bx3",
                                      "sd29.1", "sd1", "sdE3", 
                                      "sdETX4", "sd58d")))
levels(dat_for_rxnnorm$Genotype)

ggplot(dat_for_rxnnorm, aes(x=Genotype, y=lev_stat)) + 
  geom_line(aes(color=DGRP, group=DGRP)) + 
  geom_point(aes(color=DGRP)) + 
  labs(y="Levene's Statistic", x= "Mutant Allele") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# these are huge and they are all the same for sd1 then go missing for ETX4/sd58d? 
# is this due to the use of rank reduction? 
# are the massive changes due to an outlier? 
#curious...


#### Rank Reduced sd/bx might get rid of later #### 

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


m3 <- glmmTMB(lev_stat ~ Allele_1 + rr(0 + Allele_1 | WT_Background,6) +
                Replicate,
              data=sddat,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))
summary(m3)

coef(m5)$cond
ranef(m5)$cond
#still same issue is it something to do with my dataset?

plot(emmeans(m3, "Allele_1"),
     ylab = "Mutant Allele",
     xlab = "Within line variability",
     comparisons = T)



sdestimates <- coef(m3, complete = F)$cond
sdestimates2 <- data.frame(sdestimates[1],row.names = )

sdestimates_df <-  data.frame(DGRP =rownames(estimates2),
                            WT =  (estimates2[,2] + estimates2[,1]),
                            sd29.1 = (estimates2[,2]  + estimates2[,3]), 
                            sd1 = (estimates2[,2] + estimates2[,4]), 
                            sdE3 = (estimates2[,2]  + estimates2[,5]),
                            sdETX4 = (estimates2[,2]  + estimates2[,6]),
                            sd58d = (estimates2[,2]  + estimates2[,7]))


sddat_for_rxnnorm <- sdestimates_df %>% 
  pivot_longer(c(WT, sd29.1, sd1, sdE3, sdETX4, sd58d), 
               names_to = "Genotype", values_to = "lev_stat")

sddat_for_rxnnorm <- sddat_for_rxnnorm %>% 
  mutate(Genotype = factor(Genotype, 
                           levels = c("WT","sd29.1", "sd1", "sdE3", 
                                      "sdETX4", "sd58d")))
levels(sddat_for_rxnnorm$Genotype)
  
ggplot(sddat_for_rxnnorm, aes(x=Genotype, y=lev_stat)) + 
  geom_line(aes(color=DGRP, group=DGRP)) + 
  geom_point(aes(color=DGRP)) + 
  labs(y="Variability(Levene's Statistic)", x= "Mutant Allele") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# how come there is negative numbers.... levene's stat shouldn't alow for this, telling
# directionality of variability?

m4 <- glmmTMB(lev_stat ~ Allele_1 + rr(0 + Allele_1 | WT_Background,4) +
                Replicate,
              data=bxdat,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))


summary(m4)

plot(emmeans(m4, "Allele_1"),
     ylab = "Mutant Allele",
     xlab = "Within line variability",
     comparisons = T)


