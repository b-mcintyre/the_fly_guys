#Levene's statistic functions
#### Note a more common approach is to use Levene's test as you can deviates among individuals (see Dworkin 2005)
### Below I have written a general purpose function to generate Levene's deviates
LeveneDeviates <- function(y, group, med=TRUE, log_trans=TRUE) {
    
    #log transform data?
    if (log_trans==TRUE)
        y = log(y)
    
    # allows for use of mean or median as measure of central tendency
    if (med==TRUE)
        meds <- tapply(y, group, median, na.rm = TRUE)
    else 
        meds <- tapply(y, group, mean, na.rm = TRUE) 
    
    # calculates deviates for each observation from a measure of central tendency for a "group"
    abs(y - meds[group])}
    
    
# download the dll.csv data and do the following

#dll.data <- read.csv("~/Downloads/dll.csv")  # set this for yourself.  
    
    
# lev_stat <- with(dll.data, 
#     LeveneDeviates(y = SCT, group = genotype:line, med=TRUE, log_trans=FALSE))

# ls.lm <- lm(lev_stat ~ genotype*temp*line, data=dll.data)
# anova(ls.lm)  # of course this ANOVA is not particularly valid given the lack of normality, use bootstrapping or permutation.

