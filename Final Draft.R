library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(corrplot)
library(MASS)
library(blme)
library(MCMCglmm)
library(rstanarm)
library(broom.mixed)
library(dotwhisker)
library(Matrix)
library(cowplot)
remotes::install_github("glmmTMB/glmmTMB/glmmTMB#670")
library(glmmTMB)
library(lattice)


pkgs_CRAN <- c("lme4","MCMCglmm","blme",
               "pbkrtest","coda","aods3","bbmle","ggplot2",
               "reshape2","plyr","numDeriv","Hmisc",
               "plotMCMC","gridExtra","R2admb",
               "broom.mixed","dotwhisker")
install.packages(pkgs_CRAN)
rr <- "http://www.math.mcmaster.ca/bolker/R"
install.packages("glmmADMB",type="source",repos=rr)
library("devtools")


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

ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_cv, color=Allele_1, size=Individuals)) +geom_point()
# kind of looks like a quadratic relationship, where alleles with moderate phenotypic effect have the most within line
# variation 

ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_sd, color=Allele_1, size=Individuals)) + geom_point()
#See that alleles with moderate phenotypic effct on mean wing length have the largest amoug of variation

ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_sd, color=WT_Background)) + 
  geom_point() + 
  facet_wrap(~Allele_1)
# seeing by each mutant allele 

## RXN norm plot for raw mean levene's statistic

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=lev_stat)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background)) + 
  labs(y="Levene's Statistic", x= "Mutant Allele") + 
  theme(legend.position = "")

#boxplots for ORE, SdE3, and sd58d or wild type, moderate and severe mutational effects 

boxdat <- wd %>% filter(Allele_1 == c("OREw", "sd[E3]", "sd[58d]"))

ggplot(boxdat, aes(x=WT_Background, y=wing_size_mm)) + facet_wrap(~Allele_1) + 
  geom_boxplot() + theme(axis.text.x=element_text(angle=90)) + labs(x= "DGRP Line", y="Wing size mm")

# see more variation in wing size in the moderate allele than other alleles 

#### multilevel modeling ####
factor(wing_table_clean$Replicate)

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

#try nesting within both
JD2all_glm_wing_size_lev3 <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
                                 + (1|Allele_1/Replicate) + (1|WT_Background/Replicate),
                                 data = wing_table_lev_raw)

#still singular

#Original coding:
all_glm_wing_size_lev <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
                              + (1|Replicate),
                              data = wing_table_lev_raw)

#BB suggestion:

BBall_glm_wing_size_lev <- lmer(lev_stat ~ Allele_1 + (1|WT_Background/Allele_1)
                                + (1|Replicate),
                                data= wing_table_lev_raw)
#still singular 

# this is where JD said replicate should be treated as a fixed effect because there is only two levels(3ish)
# which he believes is causing the singularity/problems 

#### BB's rank reduced model with  components ####


m5 <- glmmTMB(lev_stat ~ Allele_1 + rr(0 + Allele_1 | WT_Background,5) +
                Replicate,
              data=wd,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))


# correlation matrix
v5 <- cov2cor(VarCorr(m5)$cond[[1]])

rr_mat <- matrix(v5, nrow = 9, ncol = 9, byrow = T); rr_mat

rr_cor <- cov2cor(rr_mat); rr_cor

#visualization
colnames(rr_cor) <- c("Wild Type", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")
rownames(rr_cor)<- c("Wild Type", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))


corrplot(rr_cor, type = "lower", method = "color", col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)
#### plotting variance ####

Anova(m5)

summary(m5)

fixed_effects <- fixef(m5)
fixed_effects

fixed_effects1 <-fixed_effects[1]
fixed_effects2 <- unlist(fixed_effects1)
fixed_effects3 <- fixed_effects2[c("cond.(Intercept)","cond.Allele_1bx[1]","cond.Allele_1bx[2]","cond.Allele_1bx[3]","cond.Allele_1sd[29.1]","cond.Allele_1sd[1]","cond.Allele_1sd[E3]","cond.Allele_1sd[ETX4]", "cond.Allele_1sd[58d]")]

names(fixed_effects3) <- c("sd[+]","bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")

fixed_effects3

dotchart(fixed_effects3, xlab = "Levene's statistic", ylab = "Mutant Allele")

rand_eff <-ranef(m5)
rand_eff1<- rand_eff[1]
rand_eff2<- unlist(rand_eff1)
rand_eff3 <- rand_eff2[-(121:180)]
rand_eff3                       

dotchart(rand_eff3, groups = 1:20)
## super uninformative....

?dotchart()

vc <- VarCorr(m5)
print(vc, comp = c("Std.Dev."))



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



m4 <- glmmTMB(lev_stat ~ Allele_1 + rr(0 + Allele_1 | WT_Background,4) +
                Replicate,
              data=bxdat,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))

summary(m4)
