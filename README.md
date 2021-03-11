# the_fly_guys
Investigating the Genetic Background Effects on variability in D. melanogaster 
Data set located in the .csv file tited: "NEW_CD_DGRO_Subset_Data_2019_V2"

The data is a subset of data from a previous experiment conducted in Dr. Ian Dworkin's lab by MSc Caitlyn Daley utilizing D. melanogaster. The data contains crosses of 8 mutant alleles in two allelic series in the Oregon-R genetic background with 20 wild type backgrounds from the Drosophila Genetic Research Panel (DGRP). The alleles range from very weak to very severe phenotypic effects on wing phenotype. Severity was measured quantitatively as total area in pixels and semi-quantitatively as morphological changes on a scale of 1-21. A score of 1 was given to wings appearing morphologically wild type and 21 given to wings with more severe morphological changes.

The Biological questions we are trying to answer include:
1) Do mutant alleles with weak or severe phenotypic effects display decreased sensitivity to genetic background effects
indicated by reduced variation in wing total area measurements between and among strains.
2)Do mutant alleles with moderate phenotypic effects display increased sensitivity to genetic background effects
indicated by increased phenotypic variation of wing total area.


beadex mutant allelic series (weak to moderate): bx[1], bx[2], bx[3]

scalloped mutant allelic series (weak to severe): sd[1], sd[29.1], sd[ETX4], sd[E3], sd[58d]

Oregon-R (ORE) was the background used for beadex and scalloped mutants.

The goal of the statistical test will be to measure changes in variability both within and between genetic backgrounds 
for various mutations from severe to weak. In order to test this change in variability we proposed utilizing either 
Levene’s statistical test or a Brown-Forsythe test within the context of an ANOVA. Both tests can be found in the cars 
package in R. The reason for this is because we are interested in looking at the relationship between relative variation
of mutant phenotypes and if it fits the biological model predictions including alleles of moderate phenotypic effects 
having increased phenotypic variation of wing size and SQ measure relative to the weak and severe mutant alleles. For 
these comparisons we are interested in the absolute deviation from the mean or median as a measure of variation.

All crosses were reared at 24C in the Percival incubator (RH ~ 60%) on a 12:12 hour day:night cycle (unless otherwise noted).

CD_Wing_Macro:

// Macro to compute wing area for Drosophila images. macro "WingSizeMacro2 [q]" { run("Scale...", "x=0.25 y=0.25 width=1020 height=768 interpolation=Bilinear average create"); run("Find Edges"); run("Enhance Contrast...", "saturated=0.4"); run("Sharpen"); run("Sharpen"); run("Make Binary"); run("Close-"); run("Dilate"); run("Fill Holes"); run("Despeckle"); run("Remove Outliers...", "radius=50 threshold=50 which=Dark"); run("Analyze Particles...", "size=50-Infinity pixel display summarize"); // run("Close"); //