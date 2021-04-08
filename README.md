# the_fly_guys

Investigating the Genetic Background Effects on variability in *D. melanogaster*

Data set located in the .csv file `NEW_CD_DGRO_Subset_Data_2019_V2`

The data is a subset of data from a previous experiment conducted in Dr. Ian Dworkin's lab by MSc Caitlyn Daley using *D. melanogaster*. The data contains crosses of 8 mutant alleles in two allelic series in the Oregon-R genetic background with 20 wild type backgrounds from the Drosophila Genetic Research Panel (DGRP). The alleles range from very weak to very severe phenotypic effects on wing phenotype. The beadex mutant allelic series from weak to moderate is bx[1], bx[2], bx[3]. The scalloped mutant allelic series from weak to severe is sd[1], sd[29.1], sd[ETX4], sd[E3], sd[58d]. Both he bx and the sd mutations are on the x chromosome, as such only males were examined in the F1 generation from the cross. Severity was measured quantitatively as total area in pixels through a Fiji image analysis macro and semi-quantitatively as morphological changes on a scale of 1-21. A score of 1 was given to wings appearing morphologically wild type and 21 given to wings with more severe morphological changes. 

**JD:** total area in pixels seems like it needs to be carefully standardized

The Biological questions we are trying to answer include:
1) Do mutant alleles with weak or severe phenotypic effects display decreased sensitivity to genetic background effects
(indicated by reduced variation in wing total area measurements between and among strains) ?
2) Do mutant alleles with moderate phenotypic effects display increased sensitivity to genetic background effects
(indicated by increased phenotypic variation of wing total area) ?

#Can distill the above questions into one. "Do mutant alleles with moderate phenotypic effects display increased sensitivity to genetic background effects between and within genetic wildtype backgrounds?"
#(will quantify this by looking at wing total area measurements)

**JD:** This could depend a lot on how you choose to define variation: variance? CV? something else?

Oregon-R (ORE) was the background used for beadex and scalloped mutants.

The goal of the statistical test will be to measure changes in variability both within and between genetic backgrounds 
for various mutations from severe to weak. In order to test this change in variability we proposed using either 
Levene’s statistical test or a Brown-Forsythe test within the context of an ANOVA. Both tests can be found in the `cars` 
package in R. The reason for this is because we are interested in looking at the relationship between relative variation
of mutant phenotypes and if it fits the biological model predictions including alleles of moderate phenotypic effects 
having increased phenotypic variation of wing size and SQ measure relative to the weak and severe mutant alleles. For 
these comparisons we are interested in the absolute deviation from the mean or median as a measure of variation.

**BMB**: one way or the other you will also want to quantify the *magnitude* of the effect. Classic tests like Levene's/Brown-Forsythe may be challenging to implement in a complex context (are there covariates?) 


Cleasby, Ian R., Shinichi Nakagawa, and Holger Schielzeth. “Quantifying the Predictability of Behaviour: Statistical Approaches for the Study of between-Individual Variation in the within-Individual Variance.” Methods in Ecology and Evolution 6, no. 1 (2015): 27–37. https://doi.org/10.1111/2041-210X.12281.

#Shoud switch to a HGLM/DHGLM using the CVp as our statistic?   

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


