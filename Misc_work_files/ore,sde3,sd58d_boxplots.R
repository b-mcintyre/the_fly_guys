library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(rsq)
library(emmeans)
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


#boxplots for ORE, SdE3, and sd58d 

boxdat <- wing_table_clean %>% filter(Allele_1 == c("OREw", "sd[E3]", "sd[58d]"))

ggplot(boxdat, aes(x=WT_Background, y=wing_size_mm)) + facet_wrap(~Allele_1) + 
  geom_boxplot() + theme(axis.text.x=element_text(angle=90)) + labs(x= "DGRP Line", y="Wing size mm")
