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
while (!require("glmmTMB")) {
    remotes::install_github("glmmTMB/glmmTMB")
}
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
  labs(y="Within-line variation\n(Levene's statistic)", x= "Mutant Allele") + 
  theme(legend.position = "")

#boxplots for sd1, SdE3, and sd58d for weak, moderate and severe mutational effects 

boxdat <- wd %>% filter(Allele_1 == c("sd[1]", "sd[E3]", "sd[58d]"))

ggplot(boxdat, aes(x=WT_Background, y=wing_size_mm)) + 
  facet_wrap(~Allele_1) + 
  geom_boxplot() + 
  theme(axis.text.x=element_text(angle=90, size = 9)) + 
  labs(x= "DGRP line", y="Wing size (mm)")

# see more variation in the moderate allele than other alleles 

#### multilevel modeling ####

#Original coding:
all_glm_wing_size_lev <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
                              + (1|Replicate),
                              data = wing_table_lev_raw)

#JD code suggestion results in a warning message of checking convergence 
JDall_glm_wing_size_lev <- lmer(lev_stat ~  1 + Allele_1 + Replicate + (0 + Allele_1 | WT_Background),
                       data = wing_table_lev_raw)

allFit(JDall_glm_wing_size_lev)
#majority still singular 


# this is where JD said replicate should be treated as a fixed effect because there is only two levels(3ish)
# which he believes is causing the singularity/problems 

#### BMB's rank reduced model with 5 components ####


## first try: full model
m1 <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
           + (1|Replicate),
           data = wd)

VarCorr(m1)

## this shows that among-replicate variance is zero: as JD says, three
## levels is too few for fitting RE. Replace random effect of replicate
## with fixed effect:

m2 <- update(m1, . ~ . - (1|Replicate) + Replicate)

## but, still singular.

## compound symmetry
m3 <- update(m2, . ~ . - (0 + Allele_1 | WT_Background)
             + (1|WT_Background/Allele_1))

## lev_stat ~ Allele_1 + Replicate + (1 | WT_Background/Allele_1)
## = ... + (1|WT_Background) + (1|WT_Background:Allele_1)

## expected value of allele A within WT B =
##   eps_{WT,B} + eps_{a:WT,A}
## eps_WT ~ N(0,sigma^2_WT); eps_allele:WT ~ N(0, sigma^2_{a:WT}
## for the covariance of i and j where i and j have the same WT & A
## (variance) E[(eps_{WT,B} + eps_{a:WT,A})^2] = sigma^2_WT + sigma^2_{a:WT}
## 

VarCorr(m3)
## profile confidence intervals are fairly narrow
confint(m3, parm="theta_", oldNames=FALSE)

library(bbmle)
AICtab(m1,m2,m3, logLik=TRUE)
## hmm ... this says that the more complex models are *way* better

## help("image-methods", package="Matrix")
## -> 
ifun <- function(v,blank_diag=TRUE,...) {
  nm <- rownames(v)
  nm <- gsub("Allele_","",nm)
  vv <- Matrix(v)
  if (blank_diag) diag(vv) <- NA_real_
  Matrix::image(Matrix(vv),
                useAbs=FALSE,
                scales=list(x=list(at=seq(nrow(v)),labels=nm, rot=90),
                            y=list(at=seq(ncol(v)),labels=nm)),
                xlab="",ylab="",sub="",...)
}

## correlation matrix for full model
v2 <- cov2cor(VarCorr(m2)[[1]])

## construct compound symmetric matrix
v <- c(VarCorr(m3)[["Allele_1:WT_Background"]],
       VarCorr(m3)[["WT_Background"]])
v3 <- matrix(v[2]/(v[1]+v[2]),nrow=nrow(v2),ncol=ncol(v2),
             dimnames=dimnames(v2))
diag(v3) <- 1

plot_grid(ifun(v2,main="unstructured"),ifun(v3, main="compound symm"))

## principal components 
rr <- rePCA(m2)
plot(rr[[1]]$sdev,type="b",pch=16)

m4 <- glmmTMB(formula(m2), data=wd)
v4 <- cov2cor(VarCorr(m4)$cond[[1]])

attr(v4,"blockCode") <- NULL ## strip for comparison
plot_grid(ifun(v2,main="lme4"),ifun(v4, main="glmmTMB"))
all.equal(v2,v4) ## 3% difference ...
## slightly different ... but probably close enough to trust
eigen(v2)$values
eigen(v4)$values

## rank-reduced model: 5 components
m5 <- glmmTMB(lev_stat ~ Allele_1 + rr(0 + Allele_1 | WT_Background,5) +
                Replicate,
              data=wd,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))
v5 <- cov2cor(VarCorr(m5)$cond[[1]])
plot_grid(ifun(v2,main="full"),ifun(v5, main="rr5"))
zapsmall(eigen(v5)$values)
AICtab(m1,m2,m3,m4,m5,
       mnames=c("full+rand rep","full","compsymm","TMB full+rep","rr5"),
       logLik=TRUE)

saveRDS(m5, file = "5Rank_Reduced_Model.rds")

#### Anova ####
Anova(m5, type = "III")


#### correlation matrix ####
v5 <- cov2cor(VarCorr(m5)$cond[[1]])

rr_mat <- matrix(v5, nrow = 9, ncol = 9, byrow = T); rr_mat

rr_cor <- cov2cor(rr_mat); rr_cor

#visualization of covariance matrix
colnames(rr_cor) <- rownames(rr_cor) <-
    c("OREw", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))


corrplot(rr_cor, type = "lower", method = "color", col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)

#### emmeans plot ####
plot(emmeans(m5, "Allele_1", "Replicate"),
     ylab = "Mutant Allele",
     xlab = "Within-line variability\n(Levene statistic)",
     comparison = TRUE,
     horizontal = TRUE)

#### Rxn norm of the mean plot using estimates ####

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
# is this due to the use of rank reduction I believe.

#### BMB method utilizing the model to predict the outcomes to create a reaction norm graph ####

RNdat <- with(wd,
              expand.grid(Allele_1=levels(Allele_1), WT_Background=levels(WT_Background)))

## this takes a little while
pp <- predict(m5, newdata=data.frame(RNdat, Replicate="R1"), se.fit=TRUE)

RNdat <- with(pp,
              data.frame(RNdat,
                         lev_stat=fit,
                         lev_stat_lwr=fit-2*se.fit,
                         lev_stat_upr=fit+2*se.fit))




library(colorspace)
theme_set(theme_bw())
gg1 <- (ggplot(RNdat, aes(x=Allele_1, y=lev_stat,color=WT_Background))
        + geom_line(aes(group=WT_Background))
        + geom_point()
        + geom_ribbon(aes(ymin=lev_stat_lwr, ymax=lev_stat_upr, fill=WT_Background,
                          group=WT_Background),
                      colour=NA, alpha=0.05)
        + labs(y="Within-line variability\n(Levene statistic)", x= "Mutant Allele")
        + theme(legend.position = "")
)

print(gg1 
      + scale_colour_discrete_qualitative()
      + scale_fill_discrete_qualitative()
)

library(hues)
print(gg1
      + scale_colour_iwanthue()
)

## nice idea but needs work.
library(directlabels)
direct.label(gg1, "last.bumpup")
