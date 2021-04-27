# the_fly_guys

Investigating the Genetic Background Effects on variability in *D. melanogaster*

Data set located in the .csv file `NEW_CD_DGRO_Subset_Data_2019_V2`
Code located in .R file titled "Final Draft"
Our write up on the work we completed is titled "Bhargava_and_McIntyre_Final_Project_Write_Up.pdf".

The data is a subset of data from a previous experiment conducted in Dr. Ian Dworkin's lab by MSc Caitlyn Daley using *D. melanogaster*. The data contains crosses of 8 mutant alleles in two allelic series in the Oregon-R genetic background with 20 wild type backgrounds from the Drosophila Genetic Research Panel (DGRP). The alleles range from very weak to very severe phenotypic effects on wing phenotype. The beadex mutant allelic series from  moderately weak to moderately severe is bx[1], bx[2], bx[3]. The scalloped mutant allelic series from weak to severe is sd[1], sd[29.1], sd[ETX4], sd[E3], sd[58d]. Both the bx and the sd mutations are on the x chromosome, as such only males were examined in the F1 generation from the cross. Severity was measured quantitatively through a Fiji image analysis macro which is the data we will be using
for our analysis. Each mutant allele exists in an oregon genetic background. 

The biological question we attempted to answer is how the variation around the trait mean of wing size varies between DGRP lines crossed with different mutant alleles. We hypothesize that mutant alleles with moderate phenotypic effects will display the greatest variability in within line wing size variability for dgrp lines crossed with mutant alleles. 

The goal of the statistical test will be to measure changes in variability both within DGRP backgrounds
for various mutations from weak to severe. We utilized the median form levene's statistic as a measure of within line variation and mean wing size deviations from wild type to assess the magnitude of phenotypic effect. We utilized mixed models to create predicted means for levene's statistic and wing size. We then attempted to model the effects of deviations from wild type wing size on within-line variability to understand how the mangnitude of phenotypic effect affects within-line variation. 


All crosses were reared at 24C in the Percival incubator (RH ~ 60%) on a 12:12 hour day:night cycle (unless otherwise noted).

CD_Wing_Macro:

```
// Macro to compute wing area for Drosophila images. macro "WingSizeMacro2 [q]" { run("Scale...", "x=0.25 y=0.25 width=1020 height=768 interpolation=Bilinear average create"); run("Find Edges"); run("Enhance Contrast...", "saturated=0.4"); run("Sharpen"); run("Sharpen"); run("Make Binary"); run("Close-"); run("Dilate"); run("Fill Holes"); run("Despeckle"); run("Remove Outliers...", "radius=50 threshold=50 which=Dark"); run("Analyze Particles...", "size=50-Infinity pixel display summarize"); // run("Close"); //
```

#Next steps --> decide on if we want to use a DHGLM vs a two step method
We currently plan on testing using a DHGLM model, with the comparison statistic for variability as CVp
The DHGLM works well with non-normal unbalanced data sets. 
OR
We currently plan on looking into the two-step method to see if that method is applicable to our data set, 
where we estimate within-individual variability and fit it as a response variable for a follow up analysis.

If we have time it would be interesting to utilize both models and examine the differences between the two 

 
Figure out how to incorporate the replicate blocks into the model.

Tentative variables:

Random-effects: Replicate blocks, WT DGRP strain
Fixed-effects: mutant allele
Response variable: wing size in mm

Looking into variation at 2 levels
1) Investigating variation among individuals who are genetically identical, variation around the mean
2) Investigating variation among means among each of our groups 


**BMB**: 

I think I would suggest starting with the two-stage model.  We can ask Ian about recommendations for readings.

- replicate blocks are certainly a random effect.
- I think you should include a *severity score* for the WT DGRP strain and mutant allele as a fixed effect (you could model these as quadratic), but *also* include the actual identity/label of WT and mutant allele as random effects; this allows for an overall trend but doesn't constrain the WT (for example) to fit a linear or quadratic curve *exactly* (there is a tricky idea about treating these levels as *ordered*, but it would be more complicated: we can talk about it after you have some basic stuff done)
- sex differences would be fixed, although you have to think about what things you want to allow to interact (i.e. do you expect sex differences to vary across lines etc.?)


