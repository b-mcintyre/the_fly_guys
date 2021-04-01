library(tidyverse)
library(ggplot2)

#### Cleaning Data ####

# put into rda file once sure is completely clean 

raw_wing_table <- read_csv("NEW_CD_DGRP_Subset_Data_2019_V2.csv")

str(raw_wing_table)

wing_table_clean <- raw_wing_table %>% 
  mutate(wing_length_mm = sqrt(TotalArea.px * 0.00215^2),
         TotalArea.px = NULL,
         Temperature = NULL, 
         Allele_1 = factor(Allele_1, levels = c("OREw", "sd[1]", "sd[29.1]", 
                                                "sd[ETX4]", "sd[E3]", "sd[58d]", 
                                                "bx[1]", "bx[2]", "bx[3]")), 
         WT_Background = factor(WT_Background))

levels(wing_table_clean$WT_Background)
levels(wing_table_clean$Allele_1)

str(wing_table_clean)

bxdat <- wing_table_clean %>% 
  filter(Allele_1 %in% c("OREw", "bx[1]", "bx[2]", "bx[3]")) %>%
  droplevels()

str(bxdat)

levels(bxdat$Allele_1)

sddat<- wing_table_clean %>% 
  filter(Allele_1 %in% c("OREw", "sd[1]", "sd[29.1]", "sd[58d]", "sd[E3]", "sd[ETX4]")) %>%
  droplevels() 

str(sddat)

levels(sddat$Allele_1)



#getting line means/sd/var and group sd/bx means/sd/var
#maybe plot the variance over the means? 

line_means_sd_var_cv <- wing_table_clean %>% 
  group_by(Allele_1, WT_Background) %>% 
  summarise(length_means = mean(wing_length_mm), 
            length_sd = sd(wing_length_mm), 
            length_var = var(wing_length_mm),
            length_cv = (sd(wing_length_mm)/mean(wing_length_mm))*100)

sd_means_sd_var_cv <- sddat %>% 
  group_by(Allele_1, WT_Background) %>%
  summarise(length_means = mean(wing_length_mm), 
            length_sd = sd(wing_length_mm), 
            length_var = var(wing_length_mm),
            length_cv = (sd(wing_length_mm)/mean(wing_length_mm))*100)

bx_means_sd_var_cv <- bxdat %>%
  group_by(Allele_1, WT_Background) %>%
  summarise(length_means = mean(wing_length_mm), 
            length_sd = sd(wing_length_mm), 
            length_var = var(wing_length_mm),
            length_cv = (sd(wing_length_mm)/mean(wing_length_mm))*100)

ggplot(sd_means_sd_var, aes(y=length_var, x=length_means, color=Allele_1)) + geom_point()
# interesting that the E3 and Etx4 seem to have the most variability works well with our prediction 

ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_sd, color=Allele_1, size=length_var)) +geom_point()
# kind of looks like a quadratic relationship? 

ggplot(sd_means_sd_var, aes(x=length_means, y=length_sd, color=Allele_1, size=length_var)) + geom_point()
  # one maybe two outliers in the sd29.1 data set and wild type data set? 
# maybe not log transformation is not applicable? looks like there may be a pattern there however,
# would another transformation work better? quadratic? 

ggplot(bx_means_sd_var_cv, aes(x=length_means, y=length_sd, colour=Allele_1, size=length_var)) + geom_point()
# this looks to have no observable pattern...

ggplot(line_means_sd_var_cv, aes(x=length_means, y=length_sd, color=WT_Background)) + 
  geom_point() + 
  facet_wrap(~Allele_1)
# maybe not log transform? 