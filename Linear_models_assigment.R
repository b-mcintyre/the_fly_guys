#linear models assignment

#### Prelims ####
library("lmPerm")
library("coin")
library("gtools")
library("tidyverse")
library("car")
library("emmeans")

# Please us the data set found in my repo in the main page titled: 
# NEW_CD_DGRP_Subset_Data_2019_V2

wing_table <- read.csv("NEW_CD_DGRP_Subset_Data_2019_V2.csv")

wing_table$Allele_1 <- factor(wing_table$Allele_1)

wing_table$WT_Background <- factor(wing_table$WT_Background)

# OREw as the first variable as it is wild type is a good biological starting point
# reordering the background is not applicable, no reason to have one as a zero comparison 
# point over the other 
wing_table$Allele_1 <- relevel(wing_table$Allele_1, "OREw")
levels(wing_table$Allele_1)


px.mmsqr_conversion <- 0.00005375

wing_table_mmsqr <- wing_table %>%
  mutate(
    TA_mmsqr = TotalArea.px * px.mmsqr_conversion,
    TotalArea.px = NULL
  )




#### Diagnostic plots for Allele_1 and WT_Background ####
Allele_means <- with(wing_table_mmsqr, tapply(TA_mmsqr,Allele_1,mean))
print(Allele_means)

# in model have not accounted for the wild type background variation (residuals should show
#this as unexplained variation)
lmAllele_1 <- lm(TA_mmsqr ~ Allele_1, 
         data = wing_table_mmsqr)

plot(lmAllele_1,id.n = 10)
# Residuals vs fitted see some heteroscedacity, a lot of variation, have not accounted for
#the wild type background yet, and the alleles are in an allelic series that should have
#large variation between groups 
# QQ see less conformation to normal distribution
# Scale-location see huge bump indicating some heteroscedasticity, biologically this is good
#as we hope to see changes in variance in the center, though again have not accounted for 
#wild type background or background dependence so would expect a decreased when accounted
#for wild type background, though model treats the backgrounds as each allele having the
#same effect which I believe is not true and possibly another source for heteroscedasticity 
#either way variation is not constant across the data set 
# Residuals vs leverage looks like all points are within 0.5

summary(lmAllele_1)

Allele_1_em <- emmeans(lmAllele_1,~Allele_1)
print(Allele_1_em)
plot(Allele_1_em, comparisons = TRUE)
#Not much overlap which is expected from the allelic series.

anova(lmAllele_1)
# a lot of variation between groups 

lmWT_Background <-lm(TA_mmsqr ~ WT_Background,
        data = wing_table_mmsqr)

plot(lmWT_Background, id.n=10)
# Residuals vs fitted much larger residuals because have not taken into account the variation
#due to allele_1
# QQ plot has very low values compared to expectations because have not taken into account
#a major source of variation 
# Scale-Location see a trend to decrease variation with increasing size 
# Residuals vs leverage again looks like all points are within cook's distance 


summary(lmWT_Background)

lmWT_Background_em <- emmeans(lmWT_Background,~WT_Background)
lmWT_Background_em
plot(lmWT_Background_em, comparisons = TRUE)
#see overlap of various wild type lines which is to be expected.

anova(lmWT_Background)
#see large amount of residuals unaccounted for.

lmboth <- lm(TA_mmsqr~Allele_1 + WT_Background,
             data= wing_table_mmsqr)

summary(lmboth)

length(coef(lmboth))
#should have 28 (9 Alleles, 20 DGRP lines)

Anova(lmboth)
#observe much more variation between allele groups than wild type backgrounds 

plot(lmboth)
# Residuals vs fitted looks like much less variation in residuals and very close to linear
# QQ closer to normal, more variation is being accounted for in the model
# Scale-location still see bend in the center
# All points fall within cooks distance 0.5

#### Interaction model #### 

interact <-lm(formula=TA_mmsqr ~ WT_Background * Allele_1,
        data=wing_table_mmsqr)

summary(interact)
#What I'm most interested in because I believe the Alleles and Backgrounds are interacting

plot(interact)
# Residual vs Fitted still variation but is decreasing, almost completely linear 
# QQ still not completely normal distribution 
# Scale location, observe hump in the center, still have not taken into acount background
#genetic effects
# Residuals vs Leverage, points within cooks distance 0.5

interact_em <- emmeans(interact,~Allele_1|WT_Background)
print(interact_em)
# presented as each background.


plot(interact_em,~Allele_1|WT_Background, comparisons = TRUE) + 
  facet_wrap(~WT_Background) +
  labs(x= "wing area, mm^2", y = "Mutant Allele")
# emmeans of the wild type and mutant alleles, shows 95% CI of each of the means
# see some lines have more variability that others but this is hard to interpret with
# this many lines and mutant alleles. 

length(coef(interact))
# 180 different genotypes. 

Anova(interact)
# majority of variation found between allele groups, still some residuals possibly from
# genetic background effects, in our model treated all backgrounds as having the same effect
# on mutant alleles, do not think this is the case and can possibly account for at least
# a portion of the residuals. 
