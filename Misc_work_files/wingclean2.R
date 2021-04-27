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
summary(wing_table_clean)
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

sd_means_sd_var_cv <- sddat %>% 
  group_by(Allele_1, WT_Background) %>%
  summarise(length_means = mean(wing_size_mm), 
            length_sd = sd(wing_size_mm),
            length_cv = (sd(wing_size_mm)/mean(wing_size_mm)),
            lev_stat = mean(lev_stat),
            Individuals = n()
  )

bx_means_sd_var_cv <- bxdat %>%
  group_by(Allele_1, WT_Background) %>%
  summarise(length_means = mean(wing_size_mm), 
            length_sd = sd(wing_size_mm),
            length_cv = (sd(wing_size_mm)/mean(wing_size_mm)),
            lev_stat = mean(lev_stat),
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

#RXN norm plots for line means

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=length_means)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background)) + 
  labs(y = "Wing Size mm", x = "Genetic Enviornment") + 
  theme(legend.position = "none")

ggplot(sd_means_sd_var_cv, aes(x=Allele_1, y=length_means)) +
  geom_line(aes(color=WT_Background, group=WT_Background)) +
  geom_point(aes(color=WT_Background))


ggplot(bx_means_sd_var_cv, aes(x=Allele_1, y=length_means)) +
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

#RXN norm plots for raw CV

ggplot(line_means_sd_var_cv, aes(x=Allele_1, y=length_cv)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background)) +
  labs(y = "Line CV", x = "Genetic Enviornment") + 
  theme(legend.position = "none")

ggplot(sd_means_sd_var_cv, aes(x=Allele_1, y=length_cv)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

ggplot(bx_means_sd_var_cv, aes(x=Allele_1, y=length_cv)) + 
  geom_line(aes(color=WT_Background, group=WT_Background)) + 
  geom_point(aes(color=WT_Background))

## RXN norm plots for raw mean levene's statistic

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
                    data = line_means_sd_var_cv))
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
#BB suggests using (1|background/allele) making the assumption that each of the alleles covary the same
#this seems like a very strong assumption and not what we are looking for.
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



#### Add BB's compound symetric model here ####


### DHGLM

model_mu <- DHGLMMODELING(Model="mean",
                          LinPred = Area_mmsq ~ 1 + mutant + replicate + (1|DGRP))

model_phi <- DHGLMMODELING(Model="dispersion", Link = "log",
                           LinPred = Area_mmsq ~ 1 + mutant + (1|DGRP))

dh_mod_fit <- dhglmfit(RespDist = "gaussian",
                       MeanModel = model_mu,
                       DispersionModel = model_phi,
                       DataMain = Allele_1)

... ##need to figure out how to make this work 

##### MCMCGLMM for Total Dataset ####
prior <- list( R = list(V=diag(9)/9, nu=0.004),  
               G = list(G1=list(V=diag(9)/9, nu=0.004)))

library(MCMCglmm)

MCMCglmmmTotal <- MCMCglmm(fixed =  ~ 1 + Allele_1,
                           random =~ idh(Allele_1):WT_Background,
                           rcov = ~idh(Allele_1):units,
                           data = wing_table_clean, prior=prior,
                           nitt = 20000, burnin = 5000, thin = 10)
summary(MCMCglmmmTotal)
summary(MCMCglmmmTotal$Sol) # fixed effects
summary(MCMCglmmmTotal$VCV)
s <- summary(MCMCglmmmTotal$VCV)$statistics[,"Mean"]; s
length(s)# extracting posterior means from the variances and covariances
s <- s[1:81];s 
G_mat <- matrix(s, nrow = 9, ncol = 9, byrow = T); G_mat

G_cor <- cov2cor(G_mat); G_cor # genetic variance covariance matrix for the (co)variation among DGRP across mutant alleles

#visualization
colnames(G_cor) <- c("Wild Type", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")
rownames(G_cor) <- c("Wild Type", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(G_cor, type = "lower", method = "color", col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)



###MCMCGLMM for SD Dataset

MCMCGlmmSD <- MCMCglmm(fixed = wing_size_mm ~ 1 + Allele_1,
                      random =~ us(Allele_1):WT_Background,
                      rcov = ~idh(Allele_1):units,
                      data = sddat, 
                      nitt = 20000, burnin = 5000, thin = 10)
summary(MCMCGlmmSD)
summary(MCMCGlmmSD$Sol) # fixed effects
summary(MCMCGlmmSD$VCV)
sSD <- summary(MCMCGlmmSD$VCV)$statistics[,"Mean"]; sSD
length(sSD)# extracting posterior means from the variances and covariances
sSD <- sSD[1:36];sSD
G_matSD <- matrix(sSD, nrow = 6, ncol = 6, byrow = T); G_matSD

G_corSD <- cov2cor(G_matSD); G_corSD # genetic variance covariance matrix for the (co)variation among DGRP across mutant alleles



###MCMCGLMM for BX Dataset

MCMCGlmmBX <- MCMCglmm(fixed = wing_size_mm ~ 1 + Allele_1,
                       random =~ us(Allele_1):WT_Background,
                       rcov = ~idh(Allele_1):units,
                       data = bxdat, 
                       nitt = 20000, burnin = 5000, thin = 10)
summary(MCMCGlmmBX)
summary(MCMCGlmmBX$Sol) # fixed effects
summary(MCMCGlmmBX$VCV) #random effects
sBX <- summary(MCMCGlmmBX$VCV)$statistics[,"Mean"]; sBX
length(sBX)# extracting posterior means from the variances and covariances
sBX<- sBX[1:16];sBX
G_matBX <- matrix(sBX, nrow = 4, ncol = 4, byrow = T); G_matBX

G_corBX <- cov2cor(G_matBX); G_corBX # genetic variance covariance matrix for the (co)variation among DGRP across mutant alleles



#### Simplified Model Start ####
###Changing focus from this point. Decided to simplify what we were looking at. Just going to do a GLM, looking at the fixed effects
###Going to just focus on the between mutants variation. Can do one of two things. Either we do GLM using with log transformation and 
###and that's it, or we do GLM with and without log transformation and compare the results. Might want to do AICs as well.

m1 <- lmer(wing_size_mm ~ Allele_1 
           + (1|Replicate) 
           + (1|WT_Background),
           data = wing_table_clean)


summary(m1)

VarCorr(m1)
covcor <- (VarCorr(m1[[1]]))

G_matsimple <- matrix(m5, nrow = 9, ncol = 9, byrow = T); G_matsimple

G_corsimple <- cov2cor(G_matsimple); G_corsimple


corrplot(mat, type = "lower", method = "color", col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)



install.packages("MASS")
library(MASS)
wing_table_clean$log.y <- log10( wing_table_clean$wing_size_mm)  
X<-wing_table_clean$wing_size_mm
X
X <- X + 1  
#
bc <- boxcox ( X ~ 1,
               lambda = seq ( -2, 10, 0.01 ),
               plotit = TRUE
)

bc$x
#
cbind ( bc$x, bc$y )
bc$x [ which.max ( bc$y ) ]                                       # -- optimum lambda

LAllowance <- qchisq ( 0.95, 1 ) / 2
CI.Limits  <- bc$x [ bc$y > max ( bc$y ) - LAllowance ]           # -- approx 95% confidence limits

min ( CI.Limits )                                                 # -- approx lower 95% confidence limit
max ( CI.Limits )                                                 # -- approx upper 95% confidence limit

boxcox ( CRP ~ 1,
         lambda = seq ( -2, 2, 0.01 ),
         plotit = TRUE
)


if ( Save.Graphics == TRUE  ) dev.off ( )
if ( Save.Graphics == FALSE ) Figure_03_Box_Cox_CRP_Data <- recordPlot ( )
#
if ( Save.Graphics == TRUE ) postscript ( file = '03_Box_Cox_CRP_Data.ps' )

par ( mfrow = c ( 1, 1 ) )
par ( las = 1 )
#  log-transform the observations

attach (wing_table_clean)

wing_table_clean
CRP <- CRP + 1  

Table.1.Data

round ( tapply ( wing_size_mm, Allele_1:WT_Background, mean ), 0 )                             #  actual: y ... means and SDs
round ( tapply ( wing_size_mm, Allele_1:WT_Background, sd   ), 0 )

round ( tapply ( log.y, Allele_1:WT_Background, mean ), 3 )                         #  transformed: log y ... means and SDs
round ( tapply ( log.y, Allele_1:WT_Background, sd   ), 3 )
#
# ------------------------------------------------------------------------------         # -- Table 1: last line


# -- The Log Transformation: An Overview :: Figure 1 ---------------------------         # -- Figure 1: first line
#
if ( Save.Graphics == TRUE ) postscript ( file = '01_Plankton_Data.ps' )

par ( mfrow = c ( 2, 2 ) )
par ( las   = 1          )

plot  ( Type, y,                                                  #  actual: y
        frame.plot =  FALSE,
        main = 'Table 1 data: y',
        xlab = ' ',
        ylab = ' ',
        pch  = 19,
        ylim =  c ( 0, 50000 ),
)

plot  ( Type, log.y,                                              #  transformed: log y
        frame.plot =  FALSE,
        main = 'Table 1 data: log10 ( y )',
        xlab = ' ',
        ylab = ' ',
        pch  = 19,
        ylim =  c ( 2.5, 5.0 ),
)

if ( Save.Graphics == TRUE  ) dev.off ( )
if ( Save.Graphics == FALSE ) Figure_01_Plankton_Data <- recordPlot ( )
#
# ------------------------------------------------------------------------------         # -- Figure 1: last line


# -- The Log Transformation: An Overview :: Figure 3 ---------------------------         # -- Figure 3: first line
#

# -- Actual CRP data and Box-Cox point and interval estimates of lambda --------
#
CRP <- c (  0,    3.9,  5.64, 8.22, 0,    5.62,  3.92,  6.81, 30.61,  0,
            73.2,  0,   46.7,  0,    0,    26.41, 22.82,  0,     0,     3.49,
            0,    0,    4.81, 9.57, 5.36,  0,     5.66,  0,    59.76, 12.38,
            15.74, 0,    0,    0,    0,     9.37, 20.78,  7.1,   7.89,  5.53
)
#
CRP <- CRP + 1                                                    #  offset for Box-Cox (log)


# -- Box-Cox estimates of lambda -----------------------------------------------
#
bc <- boxcox ( CRP ~ 1,
               lambda = seq ( -2, 2, 0.01 ),
               plotit = FALSE
)
#
cbind ( bc$x, bc$y )
bc$x [ which.max ( bc$y ) ]                                       # -- optimum lambda

LAllowance <- qchisq ( 0.95, 1 ) / 2
CI.Limits  <- bc$x [ bc$y > max ( bc$y ) - LAllowance ]           # -- approx 95% confidence limits

min ( CI.Limits )                                                 # -- approx lower 95% confidence limit
max ( CI.Limits )                                                 # -- approx upper 95% confidence limit


if ( Save.Graphics == TRUE ) postscript ( file = '03_Box_Cox_CRP_Data.ps' )

par ( mfrow = c ( 1, 1 ) )
par ( las = 1 )

# -- Plot Box-Cox transformation -----------------------------------------------
#
boxcox ( CRP ~ 1,
         lambda = seq ( -2, 2, 0.01 ),
         plotit = TRUE
)

if ( Save.Graphics == TRUE  ) dev.off ( )
if ( Save.Graphics == FALSE ) Figure_03_Box_Cox_CRP_Data <- recordPlot ( )
#
# ------------------------------------------------------------------------------         # -- Figure 3: last line


# -- The Log Transformation: A Possible Counterargument? :: Figure 5 -----------         # -- Figure 5: first line
#
if ( Save.Graphics == TRUE ) postscript ( file = '05_Feng_Distributions.ps' )

par ( mfrow = c ( 2, 2 ) )
par ( las   = 1          )

nObs                        <- 10000                              #  sample size from Feng 2013 response

Uniform.Data.Feng           <- round ( runif ( nObs, min = 0, max = 1 ), 4 )

Feng.Data.Figure.1.Response <- 100 * ( exp ( Uniform.Data.Feng ) - 1 ) + 1
round ( skewness ( Feng.Data.Figure.1.Response ), 3 )

Log.Feng.Data.Figure.1.Response <- log ( Feng.Data.Figure.1.Response )
round ( skewness ( Log.Feng.Data.Figure.1.Response ), 3 )

hist   ( Feng.Data.Figure.1.Response,                             #  actual: y
         main   = paste ( 'Feng Fig 1 (p 3773): n = ', nObs ),
         freq   =  FALSE,
         axes   =  FALSE,
         xlim   = c ( 0, 200 ),
         xlab   = 'y',
         ylab   = '',
         nclass = 20
)
axis ( 1 )

hist   ( Log.Feng.Data.Figure.1.Response,                         #  transformed: ln y
         main   = paste ( 'Feng Fig 1 (p 3773): n = ', nObs ),
         freq   =  FALSE,
         axes   =  FALSE,
         xlim   = c ( 0, 6 ),
         xlab   = 'ln y',
         ylab   = '',
         nclass = 20
)
axis ( 1 )

if ( Save.Graphics == TRUE  ) dev.off ( )
if ( Save.Graphics == FALSE ) Figure_05_Feng_Distributions <- recordPlot ( )
#
# ------------------------------------------------------------------------------         # -- Figure 5: last line


# -- The Log Transformation: A Possible Counterargument? :: Figure 6 -----------         # -- Figure 6: first line
#

if ( Save.Graphics == TRUE ) postscript ( file = '06_Feng_Bootstrap_Distributions.ps' )

par ( mfrow = c ( 2, 2 ) )
par ( las   = 1          )

nBoot            <- 10000
boot.statistic   <- function ( x, i ) mean ( x[i] )
boot.t.statistic <- function ( x, i ) t.test ( x[i], mu = mean ( f ) )$statistic

f                <- Feng.Data.Figure.1.Response

# -- Bootstrap the sample mean -------------------------------------------------
#
Bootstrap_Data_ybar <- boot ( f, boot.statistic, nBoot )
#
hist ( Bootstrap_Data_ybar$t,
       main   = paste ( 'Boot Feng Fig 1 data ybar :: n =', nObs ),
       ylab   = '',
       xlab   = 'Sample mean',
       nclass = 16,
       axes   = FALSE
)
#
axis ( 1 )

qqnorm ( Bootstrap_Data_ybar$t,
         main  = paste ( 'Boot Feng Fig 1 data ybar :: n =', nObs ),
         frame = FALSE
)
qqline ( Bootstrap_Data_ybar$t )

# -- Bootstrap t ---------------------------------------------------------------
#
Bootstrap_Data_t <- boot ( f, boot.t.statistic, nBoot )
#
hist ( Bootstrap_Data_t$t,
       main   = paste ( 'Boot Feng Fig 1 data t :: n =', nObs ),
       ylab   = '',
       xlab   = 't',
       nclass = 16,
       axes   = FALSE
)
#
axis ( 1 )

qqnorm ( Bootstrap_Data_t$t,
         main   = paste ( 'Boot Feng Fig 1 data t :: n =', nObs ),
         frame = FALSE
)
qqline ( Bootstrap_Data_t$t )

if ( Save.Graphics == TRUE  ) dev.off ( )
if ( Save.Graphics == FALSE ) Figure_06_Feng_Bootstrap_Distributions <- recordPlot ( )
#
# ------------------------------------------------------------------------------         # -- Figure 6: last line


# -- Get % below and above t percentiles ---------------------------------------
#
Alpha   <- 0.05                                                   #  critical significance level
t.Coeff <- qt ( Alpha / 2, nObs - 1, lower.tail = FALSE ), 4 )    #  critical value of t
t.Coeff

mean ( Bootstrap_Data_t$t < -t.Coeff ) * 100                      #  % bootstrap distribution of t < -t.Coeff
mean ( Bootstrap_Data_t$t >  t.Coeff ) * 100                      #  % bootstrap distribution of t >  t.Coeff


# -- The Log Transformation: A Possible Counterargument? :: Figure 7 -----------         # -- Figure 7: first line
#

# -- Box-Cox estimates of lambda -----------------------------------------------
#
bc <- boxcox ( Feng.Data.Figure.1.Response ~ 1,
               lambda = seq ( -2, 2, 0.01 ),
               plotit = FALSE
)

cbind ( bc$x, bc$y )
bc$x [ which.max ( bc$y ) ]                                       # -- optimum lambda

LAllowance <- qchisq ( 0.95, 1 ) / 2
CI.Limits  <- bc$x [ bc$y > max ( bc$y ) - LAllowance ]           # -- approx 95% confidence limits

min ( CI.Limits )                                                 # -- approx lower 95% confidence limit
max ( CI.Limits )                                                 # -- approx upper 95% confidence limit


if ( Save.Graphics == TRUE ) postscript ( file = '07_Box_Cox_Feng_Data.ps' )

par ( mfrow = c ( 1, 1 ) )
par ( las   = 1          )

# -- Plot Box-Cox transformation -----------------------------------------------
#
boxcox ( Feng.Data.Figure.1.Response ~ 1,
         lambda = seq ( -2, 2, 0.01 ),
         plotit = TRUE
)

if ( Save.Graphics == TRUE  ) dev.off ( )
if ( Save.Graphics == FALSE ) Figure_07_Box_Cox_Feng_Data <- recordPlot ( )
#
# ------------------------------------------------------------------------------         # -- Figure 7: last line


# -- Practical Considerations :: Figure 8 --------------------------------------         # -- Figure 8: first line
#

if ( Save.Graphics == TRUE ) postscript ( file = '08_Plankton_Residuals.ps' )

par ( mfrow = c ( 2, 2 ) )
par ( las   = 1          )

lm.y <- lm ( y ~ factor ( Type ) )                                # -- statistical model: actual plankton data

plot ( Type, resid ( lm.y ),                                      # -- residual plot: actual plankton data
       main  = 'Actual: residuals vs Plankton Type',
       ylim  = c ( -15000, +15000 ),
       pch   = 19,
       col   = 'black',
       frame = FALSE
)
abline ( h = 0, col = 'gray75' )


lm.log.y <- lm ( log.y ~ factor ( Type ) )                        # -- statistical model: transformed plankton data

plot ( Type, resid ( lm.log.y ),                                  # -- residual plot: transformed plankton data
       main  = 'Log: residuals vs Plankton Type',
       ylim  = c ( -0.30, +0.30 ),
       pch   = 19,
       col   = 'black',
       frame = FALSE
)
abline ( h = 0, col = 'gray75' )

if ( Save.Graphics == TRUE  ) dev.off ( )
if ( Save.Graphics == FALSE ) Figure_08_Plankton_Residuals <- recordPlot ( )
#
# ------------------------------------------------------------------------------         # -- Figure 8: last line


# -- Replay each data graphic --------------------------------------------------
#
if ( Save.Graphics == FALSE ) replayPlot ( Figure_01_Plankton_Data                )
if ( Save.Graphics == FALSE ) replayPlot ( Figure_03_Box_Cox_CRP_Data             )
if ( Save.Graphics == FALSE ) replayPlot ( Figure_05_Feng_Distributions           )
if ( Save.Graphics == FALSE ) replayPlot ( Figure_06_Feng_Bootstrap_Distributions )
if ( Save.Graphics == FALSE ) replayPlot ( Figure_07_Box_Cox_Feng_Data            )
if ( Save.Graphics == FALSE ) replayPlot ( Figure_08_Plankton_Residuals           )