#TITLE: Code for data checking, cleaning, coding an averages table based on genotypes, and plotting.  
#code tested on R-4.0.3

## JD: It's good to have good file names, but let us know the filename for your weekly submission
## I removed the comma from the name of your other R script (while I was confused). Commas in filenames seem non-standard and thus potentially bad.

#load tidyverse ## JD: changed from install
library(tidyverse)

#ensure tidyr and ggplot2 are in the library
## JD: No, this is just a waste of time.
## Also, you don't need to comment everything you do ☺
## library(tidyr)
## library(ggplot2)

# create a wing_table
Wing_Table <- read_csv("NEW_CD_DGRP_Subset_Data_2019_V2.csv")

###############################################################################
# Beginning of data checking and cleaning


summary(Wing_Table)
options(dplyr.summarise.inform=FALSE) 

#make list of what occurs with what and how many times
## JD: Not clear the purpose here: you are not saving these, and only looking at the top of the output
print(Wing_Table %>%
        group_by(Allele_1, WT_Background)%>%
        summarise(count = n())
      )

print(Wing_Table %>%
        group_by(Allele_1, WT_Background, Replicate) %>%
      summarise(count = n())
      )


## JD: Use unite with care: it's often better practice to keep the logical columns and just group by your "key" columns.
## It's OK to use for convenience if you're just checking as below
# Create genotype and add ID numbers to give relational info to raw data
Wing_Table_Geno <- Wing_Table %>% unite(genotype, Allele_1,WT_Background) %>%
  mutate(obs=seq(n()))

#Add ID numbers to each row
#Wing_Table_Geno$ID <-1:nrow(Wing_Table_Geno)
  
summary(Wing_Table_Geno)

## Nice; you can also use stopifnot here to flag in case an error is introduced
# Ensure at least 2 replicates per genotype
print(Wing_Table_Geno %>%
        group_by(Replicate, genotype)%>%
        summarise(count= n()) %>%
        filter(n()< 2)
      )


## Check if number of replicates is the same across each genotype
print(Wing_Table_Geno %>%  
        group_by(Replicate) %>%
        summarise(count = n()))

## what is R2a? Go back to raw data input and find out
## JD: ☺

R2a <- Wing_Table_Geno %>%
  group_by(Replicate,genotype) %>%
  summarise(count = n())

R2a <- dplyr::filter(R2a, Replicate == "R2a") #count of number R2a individuals per genotype 


R2a_pos <- Wing_Table_Geno %>%
  group_by(Replicate, genotype, obs) %>%
  summarise()

R2a_pos <- dplyr::filter(R2a_pos, Replicate == "R2a") # positional info on R2a individuals


## what is R2a? Go back to raw data input and find out
## R2a is a small replicate block added

# see number of replicate individuals 
print(Wing_Table_Geno %>%
        group_by(genotype,Individuals) %>%
      summarise(count = n ())) %>%
      filter(count>2)
# appears only one genotype is in replicate blocks 1,2, and 2a 


# Find out if there is any duplicates of the wing area in px which is unlikely
# due to the method of calculation so more likely an input error.
# not applicable to SQ because they were measured semi-quantitatively so I would 
# have many duplicates
dup_dat <- Wing_Table_Geno %>% 
        group_by(TotalArea.px,) %>%
        filter(n()>1)

dup_dat <- dplyr::select(dup_dat, genotype, TotalArea.px, obs)
# double check wings to see if they are duplicate or if very unlikely 

# write a csv in order to print and reference for manual check 
## JD: provide a file extension (see change below)
## JD: ./ doesn't help here – it just means write it where you were already planning to
write.csv(dup_dat, "dup_dat.csv")

# convert total area in pixels to total area in mm^2
## JD: Nice use of a name to annotate what you're doing
px.mmsqr_conversion <- 0.00005375
Wing_Table_Geno_mmsqr <- Wing_Table_Geno %>%
  mutate(
    TA_mmsqr = TotalArea.px * px.mmsqr_conversion,
    TotalArea.px = NULL
  )

###############################################################################
# now trying to get average table showing both the SQ and the total area


Wing_Table_Geno_av <- Wing_Table_Geno_mmsqr %>%
  group_by(genotype)%>%
  dplyr::summarise(across(c(SQ_Measure, TA_mmsqr), mean)) %>%
  dplyr::rename(SQ_avg = SQ_Measure, TA_avg_mmsqr = TA_mmsqr) %>%
  print(Wing_Table_Geno_av)


summary(Wing_Table_Geno_av)


###############################################################################
###############################################################################
# plots+graphs code

# see how total area follows SQ measures 
TA_avgVsSQ_avgp <- plot (TA_avg_mmsqr~SQ_avg, data=Wing_Table_Geno_av)
# good that it follows a downwards trend makes sense

# see distribution of wing size over genotype looks like
geno_wsb <- ggplot(data = Wing_Table_Geno_av,aes(x=genotype, y=TA_avg_mmsqr)) + 
  geom_bar(stat = "identity")

print(geno_wsb)

# Same for SQ
geno_sqb <- ggplot(data = Wing_Table_Geno_av, aes(x=genotype, y=SQ_avg)) +
  geom_bar (stat = "identity")

print(geno_sqb)

# interesting that distribution is different between the two measurements 

## Grade 2.2/3

###############################################################################
###############################################################################
#TITLE: Code for data checking, cleaning, coding an averages table based on genotypes, and plotting.  
#code tested on R-4.0.3

## JD: Um, what's going on here; did you accidentally just save the code twice? Please figure out and clean up.

#download tidyverse
library(tidyverse)

#ensure tidyr and ggplot2 are in the library
library(tidyr)
library(ggplot2)

# create a wing_table
Wing_Table <- read_csv("NEW_CD_DGRP_Subset_Data_2019_V2.csv")

###############################################################################
# Beginning of data checking and cleaning


summary(Wing_Table)
options(dplyr.summarise.inform=FALSE) 

#make list of what occurs with what and how many times
print(Wing_Table %>%
        group_by(Allele_1, WT_Background)%>%
        summarise(count = n())
      )

print(Wing_Table %>%
        group_by(Allele_1, WT_Background, Replicate) %>%
      summarise(count = n())
      )


# Create genotype and add ID numbers to give relational info to raw data
Wing_Table_Geno <- Wing_Table %>% unite(genotype, Allele_1,WT_Background) %>%
  mutate(obs=seq(n()))

#Add ID numbers to each row
#Wing_Table_Geno$ID <-1:nrow(Wing_Table_Geno)
  
summary(Wing_Table_Geno)


# Ensure at least 2 replicates per genotype
print(Wing_Table_Geno %>%
        group_by(Replicate, genotype)%>%
        summarise(count= n()) %>%
        filter(n()< 2)
      )


## Check if number of replicates is the same across each genotype
print(Wing_Table_Geno %>%  
        group_by(Replicate) %>%
        summarise(count = n()))

## what is R2a? Go back to raw data input and find out

R2a <- Wing_Table_Geno %>%
  group_by(Replicate,genotype) %>%
  summarise(count = n())

R2a <- dplyr::filter(R2a, Replicate == "R2a") #count of number R2a individuals per genotype 


R2a_pos <- Wing_Table_Geno %>%
  group_by(Replicate, genotype, obs) %>%
  summarise()

R2a_pos <- dplyr::filter(R2a_pos, Replicate == "R2a") # positional info on R2a individuals


## what is R2a? Go back to raw data input and find out
## R2a is a small replicate block added

# see number of replicate individuals 
print(Wing_Table_Geno %>%
        group_by(genotype,Individuals) %>%
      summarise(count = n ())) %>%
      filter(count>2)
# appears only one genotype is in replicate blocks 1,2, and 2a 


# Find out if there is any duplicates of the wing area in px which is unlikely
# due to the method of calculation so more likely an input error.
# not applicable to SQ because they were measured semi-quantitatively so I would 
# have many duplicates
dup_dat <- Wing_Table_Geno %>% 
        group_by(TotalArea.px,) %>%
        filter(n()>1)

dup_dat <- dplyr::select(dup_dat, genotype, TotalArea.px, obs)
# double check wings to see if they are duplicate or if very unlikely 

# write a csv in order to print and reference for manual check 
write.csv(dup_dat, "./dup_dat")

# convert total area in pixels to total area in mm^2
px.mmsqr_conversion <- 0.00005375
Wing_Table_Geno_mmsqr <- Wing_Table_Geno %>%
  mutate(
    TA_mmsqr = TotalArea.px * px.mmsqr_conversion,
    TotalArea.px = NULL
  )

###############################################################################
# now trying to get average table showing both the SQ and the total area


Wing_Table_Geno_av <- Wing_Table_Geno_mmsqr %>%
  group_by(genotype)%>%
  dplyr::summarise(across(c(SQ_Measure, TA_mmsqr), mean)) %>%
  dplyr::rename(SQ_avg = SQ_Measure, TA_avg_mmsqr = TA_mmsqr) %>%
  print(Wing_Table_Geno_av)


summary(Wing_Table_Geno_av)


###############################################################################
###############################################################################
# plots+graphs code

# see how total area follows SQ measures 
TA_avgVsSQ_avgp <- plot (TA_avg_mmsqr~SQ_avg, data=Wing_Table_Geno_av)
# good that it follows a downwards trend makes sense

# see distribution of wing size over genotype looks like
geno_wsb <- ggplot(data = Wing_Table_Geno_av,aes(x=genotype, y=TA_avg_mmsqr)) + 
  geom_bar(stat = "identity")

print(geno_wsb)

# Same for SQ
geno_sqb <- ggplot(data = Wing_Table_Geno_av, aes(x=genotype, y=SQ_avg)) +
  geom_bar (stat = "identity")

print(geno_sqb)
# interesting that distribution is different between the two measurements 
###############################################################################
###############################################################################
