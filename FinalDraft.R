library(plyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(corrplot)
library(MASS)
## library(blme)
library(rstanarm)
library(broom.mixed)
library(dotwhisker)
library(Matrix)
library(cowplot)
library(emmeans)
library(effects)
## don't auto-install
library(glmmTMB)
if (packageVersion("glmmTMB") < "1.0.2.9000") {
  stop("Need devel version of glmmTMB for rr() model: use 'remotes::install_github(\"glmmTMB/glmmTMB/glmmTMB\")'")
}
library(lattice)
library(coda)

## JD: library() should be at top
library(hues)
library(bbmle)
library(colorspace)

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


## RXN norm plot for raw mean levene's statistic

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=lev_stat)) +
  geom_line(aes(color=WT_Background, group=WT_Background)) +
  geom_point(aes(color=WT_Background)) +
  labs(y="Within-line variation\n(Levene's statistic)", x= "Mutant Allele") +
  theme(legend.position = "")

## RXN orm plot for raw mean wing size

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=length_means)) +
  geom_line(aes(color=WT_Background,group=WT_Background)) +
  geom_point(aes(color=WT_Background)) +
    labs(y="Line means\n (wing_size_mm)", x= "Mutant Allele") +
  theme(legend.position = "")

#boxplots for sd1, SdE3, and sd58d for weak, moderate and severe mutational effects

boxdat <- wd %>% filter(Allele_1 == c("sd[1]", "sd[E3]", "sd[58d]"))

ggplot(boxdat, aes(x=WT_Background, y=wing_size_mm)) +
  facet_wrap(~Allele_1) +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle=90, size = 9)) +
  labs(x= "DGRP line", y="Wing size (mm)")

# see more variation in the moderate allele than other alleles

##Checking to see if log transformation is appropriate for levene's statistic values or not by using Box-Cox analysis


str(wing_table_lev_raw)

X<-wing_table_lev_raw$lev_stat
X
X <- X + 1    ##offset for Box-Cox (log)

bc <- boxcox ( X ~ 1,
               lambda = seq ( -20, 10, 0.01 ),
               plotit = TRUE
)

bc$x
cbind ( bc$x, bc$y )
bc$x [ which.max ( bc$y ) ]                                       # -- optimum lambda

LAllowance <- qchisq ( 0.95, 1 ) / 2
CI.Limits  <- bc$x [ bc$y > max ( bc$y ) - LAllowance ]           # -- approx 95% confidence limits

min ( CI.Limits )                                                 # -- approx lower 95% confidence limit
max ( CI.Limits )                                                 # -- approx upper 95% confidence limit

##Checking to see if log transformation is appropriate for wing size values or not, using Box-Cox analysis
X2<-wing_table_lev_raw$wing_size_mm
## X2 ## BMB: don't print ...
X2 <- X2 + 1    ##offset for Box-Cox (log)
bc2 <- boxcox ( X2 ~ 1,
                lambda = seq ( -20, 10, 0.01 ),
                plotit = TRUE
)

bc2$x
cbind ( bc2$x, bc2$y )
bc2$x [ which.max ( bc2$y ) ]                                       # -- optimum lambda

LAllowance <- qchisq ( 0.95, 1 ) / 2
CI.Limits  <- bc2$x [ bc2$y > max ( bc2$y ) - LAllowance ]           # -- approx 95% confidence limits

min ( CI.Limits )                                                 # -- approx lower 95% confidence limit
max ( CI.Limits )                                                 # -- approx upper 95% confidence limit


#### multilevel modeling ####

#Original coding:
all_glm_wing_size_lev <- lmer(lev_stat ~ Allele_1 + (0 + Allele_1 | WT_Background)
                              + (1|Replicate),
                              data = wing_table_lev_raw)

#JD code suggestion results in a warning message of checking convergence
JDall_glm_wing_size_lev <- lmer(lev_stat ~  1 + Allele_1 + Replicate + (0 + Allele_1 | WT_Background),
                       data = wing_table_lev_raw)

## BMB: might want to skip/comment out this step, as it takes a long time and isn't actually used in any
##  of the subsequent steps. When people like me or JD try to rerun the code, it's tedious. (And you *should*
##  be running your code start-to-finish periodically to make sure everything works consistently ...)
if (FALSE) {
  allFit(JDall_glm_wing_size_lev)
}
#majority still singular
## BMB: singularity is **not** a problem you should expect to fit with a different optimizer.
## Failure to convergence is a property of/problem with _numerical optimization/model fitting_
## Singularity is a property of the _model & data_

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

print(gg1
      + scale_colour_iwanthue()
)

#### mean wing size ####

m6 <- lmer(wing_size_mm ~ Allele_1 + (0 + Allele_1 | WT_Background)
           + (1|Replicate),
           data = wd)

VarCorr(m6)

## this shows that among-replicate variance is zero: as JD says, three
## levels is too few for fitting RE. Replace random effect of replicate
## with fixed effect:

m7 <- update(m6, . ~ . - (1|Replicate) + Replicate)

## but, still singular.

## compound symmetry
m8 <- update(m7, . ~ . - (0 + Allele_1 | WT_Background)
             + (1|WT_Background/Allele_1))

VarCorr(m7)
## profile confidence intervals are fairly narrow
## BMB: consider skipping or caching? (slow)
if (FALSE) confint(m7, parm="theta_", oldNames=FALSE)

AICtab(m6,m7,m8, logLik=TRUE)
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
v6 <- cov2cor(VarCorr(m7)[[1]])

## construct compound symmetric matrix
v7 <- c(VarCorr(m8)[["Allele_1:WT_Background"]],
       VarCorr(m8)[["WT_Background"]])
v8 <- matrix(v7[2]/(v7[1]+v7[2]),nrow=nrow(v6),ncol=ncol(v6),
             dimnames=dimnames(v6))
diag(v8) <- 1

plot_grid(ifun(v7,main="unstructured"),ifun(v8, main="compound symm"))

## principal components
rr1 <- rePCA(m7)
plot(rr1[[1]]$sdev,type="b",pch=16)

m9 <- glmmTMB(formula(m7), data=wd)
v9 <- cov2cor(VarCorr(m9)$cond[[1]])

attr(v9,"blockCode") <- NULL ## strip for comparison
plot_grid(ifun(v6,main="lme4"),ifun(v9, main="glmmTMB"))
all.equal(v6,v9) ## 3% difference ...
## slightly different ... but probably close enough to trust
eigen(v6)$values

## rank-reduced model: 3 components
m10 <- glmmTMB(wing_size_mm ~ Allele_1 + rr(0 + Allele_1 | WT_Background,3) +
                Replicate,
              data=wd,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))
v10 <- cov2cor(VarCorr(m10)$cond[[1]])
plot_grid(ifun(v2,main="full"),ifun(v5, main="rr5"))

zapsmall(eigen(v5)$values)
AICtab(m6,m7,m8,m9,m10,
       mnames=c("full+rand rep","full","compsymm","TMB full+rep","rr5"),
       logLik=TRUE)


#added correlation plot for rr3 model
rr_mat1 <- matrix(v10, nrow = 9, ncol = 9, byrow = T); rr_mat1

rr_cor1 <- cov2cor(rr_mat1); rr_cor1

#visualization
colnames(rr_cor1) <- rownames(rr_cor1) <-
  c("OREw", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")

corrplot(rr_cor1, type = "lower", method = "color", col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)

# Anova
Anova(m10, type = "III")

#emmeans
plot(emmeans(m10, "Allele_1", "Replicate"),
     ylab = "Mutant Allele",
     xlab = "Wing size (mm)",
     comparison = TRUE,
     horizontal = TRUE)

#RXNM
RNdat1 <- with(wd,
               expand.grid(Allele_1=levels(Allele_1), WT_Background=levels(WT_Background)))
pp1 <- predict(m10, newdata=data.frame(RNdat1, Replicate="R1"), se.fit=TRUE)

RNdat1 <- with(pp1,
               data.frame(RNdat1,
                          wing_size_mm=fit,
                          wing_size_mm_lwr=fit-2*se.fit,
                          wing_size_mm_upr=fit+2*se.fit))
theme_set(theme_bw())
gg2 <- (ggplot(RNdat1, aes(x=Allele_1, y=wing_size_mm,color=WT_Background))
        + geom_line(aes(group=WT_Background))
        + geom_point()
        + geom_ribbon(aes(ymin=wing_size_mm_lwr, ymax=wing_size_mm_upr, fill=WT_Background,
                          group=WT_Background),
                      colour=NA, alpha=0.05)
        + labs(y="Wing size (mm)", x= "Mutant Allele")
        + theme(legend.position = "")
)

print(gg2
      + scale_colour_discrete_qualitative()
      + scale_fill_discrete_qualitative()
)

# Deviations code
wRNdat <- with(wd,
               expand.grid(Allele_1=levels(Allele_1), WT_Background=levels(WT_Background)))

## PREDICTED wing size for each combination
wRNdat <- with(pp,
               data.frame(wRNdat,
                          wing_size_mm=fit))

wideRNdat <- wRNdat %>% pivot_wider(names_from = "Allele_1", "WT_Background", values_from = "wing_size_mm")

#deviations by mutant allele
bx1 <- wideRNdat[,3] -  wideRNdat[,2]
bx2 <-  wideRNdat[,4] -  wideRNdat[,2]
bx3 <-  wideRNdat[,5] -  wideRNdat[,2]
sd291 <-  wideRNdat[,6] -  wideRNdat[,2]
sd1 <-  wideRNdat[,7] -  wideRNdat[,2]
sde3 <-  wideRNdat[,8] -  wideRNdat[,2]
sdetx4 <-  wideRNdat[,9] -  wideRNdat[,2]
sd58d <-  wideRNdat[,10] -  wideRNdat[,2]
## BMB: could do this with a mutate(across(...))
wing_size_deviations <- cbind(wideRNdat[,1, drop = FALSE],bx1,bx2,bx3,sd291,sd1,sde3,sdetx4,sd58d)

wing_size_deviations <- wing_size_deviations %>%
  pivot_longer(c("bx[1]", "bx[2]", "bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]"),
               names_to = "Genotype", values_to = "Deviations")

wing_size_deviations$WT_Background <- as.factor(wing_size_deviations$WT_Background)

modRNdat <- RNdat %>% filter(Allele_1 != "OREw")

## attach lev_stat
wing_size_deviations <- cbind(wing_size_deviations,modRNdat[,3, drop = FALSE])

gg3 <- ggplot(wing_size_deviations,
       aes(y = lev_stat, x = Deviations, colour = Genotype)) +
  geom_point( size = 2, alpha = 0.5) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Deviations From Wild type, mm",
       y = "Within-Line Variability\n(Levene Stat)")

gg3

cor.test(wing_size_deviations$Deviations,
         wing_size_deviations$lev_stat,
         method = "pearson",
         conf.level = 0.95)

#Mixed modeling

deviationlmer <- lmer(lev_stat ~  1 + Deviations +
                        Genotype  + (1|WT_Background),
                      data = wing_size_deviations)

allFit(deviationlmer)

coef(deviationlmer)

wing_size_deviations <- wing_size_deviations %>%
  mutate( scaledDevs = scale(wing_size_deviations$Deviations))


deviationlmer1 <- lmer(lev_stat ~  1  + scaledDevs +
                         Genotype  + (1|WT_Background),
                       data = wing_size_deviations)

allFit(deviationlmer1)

fixef(deviationlmer)
fixef(deviationlmer1)
#scaling doesn't seem to make much of a difference
ranef(deviationlmer)
ranef(deviationlmer1)

llikAIC(deviationlmer, deviationlmer1)

deviationlmer2 <- lmer(lev_stat ~ 1 + Deviations + (1|WT_Background),data = wing_size_deviations)

allFit(deviationlmer2)


deviationlm <- lm(lev_stat ~ 1 + Deviations + Genotype + WT_Background,
                       data = wing_size_deviations)

summary(deviationlm)

with(wd, table(Allele_1,WT_Background, Replicate))

## BMB: an alternate calculation, based on raw data.  Trying to figure out what the difference is ... ???
wd2 <- (wd
  %>% filter(Replicate=="R1")
  %>% dplyr::select(Allele_1,WT_Background,lev_stat,wing_size_mm, Replicate)
  %>% bind_cols(wing_size_pred=predict(m10, newdata=.))
)

  %>% mutate(across(c(wing_size_mm, wing_size_pred), ~. - wing_size_mm[Allele_1=="OREw"])))
  %>% filter(Allele_1 != "OREw")
)

ggplot(wd2, aes(wing_size_pred, lev_stat, colour=Allele_1)) + geom_point(alpha=0.5)

## BMB: giving up now, I can't figure out what's going on ...
