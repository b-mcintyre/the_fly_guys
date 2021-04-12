library(tidyverse)
library(lme4)
library(blme)
library(MCMCglmm)
library(rstanarm)
library(broom.mixed)
library(dotwhisker)
## 
library(Matrix) ## image plots
library(cowplot) ## plot_grid

wd <- readRDS("wing_table_lev_raw.rds")

## first try: full model
m1 <- lmer(wing_size_mm ~ Allele_1 + (0 + Allele_1 | WT_Background)
           + (1|Replicate),
           data = wd)

VarCorr(m1)

## this shows that among-replicate variance is zero: as JD says, three
## levels is too few for fitting RE. Replace random effect of replicate
## with fixed effect:

m2 <- update(m1, . ~ . - (1|Replicate) + Replicate)

## but, still singular.

## compound symmetry
m3 <- update(m2, . ~ . - (0 + Allele_1 | WT_Background)
             + (1|WT_Background/Allele_1))

## lev_stat ~ Allele_1 + Replicate + (1 | WT_Background/Allele_1)
## = ... + (1|WT_Background) + (1|WT_Background:Allele_1)

## expected value of allele A within WT B =
##   eps_{WT,B} + eps_{a:WT,A}
## eps_WT ~ N(0,sigma^2_WT); eps_allele:WT ~ N(0, sigma^2_{a:WT}
## for the covariance of i and j where i and j have the same WT & A
## (variance) E[(eps_{WT,B} + eps_{a:WT,A})^2] = sigma^2_WT + sigma^2_{a:WT}
## 

VarCorr(m3)
## profile confidence intervals are fairly narrow
confint(m3, parm="theta_", oldNames=FALSE)

library(bbmle)
AICtab(m1,m2,m3, logLik=TRUE)
## hmm ... this says that the more complex models are *way* better

## help("image-methods", package="Matrix")
## -> 
ifun <- function(v,blank_diag=TRUE,...) {
    nm <- rownames(v)
    nm <- gsub("Allele_","",nm)
    vv <- Matrix(v)
    if (blank_diag) diag(vv) <- NA_real_
    Matrix::image(Matrix(vv),
                  useAbs=FALSE,
                  scales=list(x=list(at=seq(nrow(v)),labels=nm, rot=90),
                              y=list(at=seq(ncol(v)),labels=nm)),
                  xlab="",ylab="",sub="",...)
}

## correlation matrix for full model
v2 <- cov2cor(VarCorr(m2)[[1]])

## construct compound symmetric matrix
v <- c(VarCorr(m3)[["Allele_1:WT_Background"]],
       VarCorr(m3)[["WT_Background"]])
v3 <- matrix(v[2]/(v[1]+v[2]),nrow=nrow(v2),ncol=ncol(v2),
             dimnames=dimnames(v2))
diag(v3) <- 1

plot_grid(ifun(v2,main="unstructured"),ifun(v3, main="compound symm"))

## principal components 
rr <- rePCA(m2)
plot(rr[[1]]$sdev,type="b",pch=16)

## install from pull request (need compilation tools installed)
remotes::install_github("glmmTMB/glmmTMB/glmmTMB#670")
library(glmmTMB)

m4 <- glmmTMB(formula(m2), data=wd)
v4 <- cov2cor(VarCorr(m4)$cond[[1]])

attr(v4,"blockCode") <- NULL ## strip for comparison
plot_grid(ifun(v2,main="lme4"),ifun(v4, main="glmmTMB"))
all.equal(v2,v4) ## 3% difference ...
## slightly different ... but probably close enough to trust
eigen(v2)$values
eigen(v4)$values

## rank-reduced model: 5 components
m5 <- glmmTMB(wing_size_mm ~ Allele_1 + rr(0 + Allele_1 | WT_Background,5) +
                  Replicate,
              data=wd,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))
v5 <- cov2cor(VarCorr(m5)$cond[[1]])
plot_grid(ifun(v2,main="full"),ifun(v5, main="rr5"))
zapsmall(eigen(v5)$values)
AICtab(m1,m2,m3,m4,m5,
       mnames=c("full+rand rep","full","compsymm","TMB full+rep","rr5"),
       logLik=TRUE)

## something a little funny going on here ...
