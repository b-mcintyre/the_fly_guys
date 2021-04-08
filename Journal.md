
JD: Why is journal empty?

began with calculating the means, standard deviations, coefficients of variations, indviduals counts, and levene's statistic for grouped Allele_1 and wild type background

created two sepearate data frames with the two differet allelic series of bx and sd

plotted the sd vs the means for each of the grouped lines

can maybe fit a quadratic model to the model? - as talked to by BB

Goal --> find the different between individuals within group/line variabilites   (sd) (commonly called within individual variation)
SEE if there is genetic variation(or change in sd) in how much variation(sd) there is in among those individuals within group/lines 
Does variability or (more standard deviations on standard deviations) increase on mutants with moderate phenotypic effects? 


i.e. 1) measure change in mean between environments(mutations): go from OREw --> sd58d how the line means change
<<<<<<< HEAD
2) measure change in standard deviation between environments: go from OREw --> sd58d how does the line standard deviations change

looking at variation between individuals within lines/genotypes
looking at variation between lines 

utilize median form of levene's statistic to give a measurement of variability among individuals who are genetically similar, median form of levene's statistic does not rely on the normality of data 
Test for equality of levene's statistic using an F-test
can utilize means and standard deviation/variance to see how among/between line means

for a given quantiative trait experimentally measured across conditions/treatments the mean and more generally the distribution can vary among data sets(groups)
- need appropriate measures of the measurements between mean and variance
- method of doing this = regress variance over the mean of the dataset then use the residual deviations from this relationship 

to compare robustness of different traits, measured with different phenotypic scales to the same pertubation(our case mutant enviornent/allele) need a dimensionless metric
CV often but this only applies to normal distributions with non-negative values 
Comparing CV on different traits remains problematic when relationship cannot be assumed to be proportional 


Fixed effects = mutants (x9)
random effects = DGRP lines (x20) and replicate blocks (x3)

=======
2) measure change in standard deviation between envionrments: go from OREw --> sd58d how does the line standard deviations change

looking at variation between individual within lines/genotypes
looking at variation between lines 


7 April 21

- create reaction norm plots via utilization of spaghetti plots
x-axis = mutant allele
color = DGRP line
y-axis = line means 

This will allow me to see the change in means indicating that the lines means are changing due to genetic background effects in conjunction with mutational effects

- create more reaction norm plots
x-axis = mutant allele
color = DGRP line
y-axis = our measure of variation (CV or levene's statistic)

- can see that while means of lines can change so can the amount of variation around that mean depending on the genetic background and mutation 


do for all lines then separate plots via allelic series bx and sd

Have multi-level dataset --> can do 2 step, use a multilevel model (like Double hiearchy?)
see if can do a two step model

1) fit a regression to the within group between individuals
2) fit a regression to the among groups

- Will this be the same predictors for both levels of the model = no at the second level the predictors are the models' that were created for each of the predictors
- What are my predictors for this model = DGRP line, mutant allele, and Replicate blocks

- question then becomes what models do I fit at the 1st levels to bring up to the second level (how do I fit these models how does that make sense??)


>>>>>>> 11afad4bbc6a005fb41b1b49d67afe76485f6bda

