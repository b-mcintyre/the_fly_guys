#AHHH FUCK

m6 <- glmmTMB(wing_size_mm ~ Allele_1 + rr(0 + Allele_1 | WT_Background,5) +
                Replicate,
              data=wd,
              control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))

saveRDS(m6, file = "5Rank_Reduced_Model.rds")

Anova(m6, type = "III")

v6 <- cov2cor(VarCorr(m6)$cond[[1]])

rr_mat1 <- matrix(v6, nrow = 9, ncol = 9, byrow = T); rr_mat1

rr_cor1 <- cov2cor(rr_mat1); rr_cor1

colnames(rr_cor1) <- rownames(rr_cor1) <-
  c("OREw", "bx[1]","bx[2]","bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]")

col1 <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))


corrplot(rr_cor1, type = "lower", method = "color", col=col(200),
         addCoef.col = "black",
         tl.col = "black", tl.srt = 45)

#### emmeans plot ####
plot(emmeans(m6, "Allele_1", "Replicate"),
     ylab = "Mutant Allele",
     xlab = "Wing size (mm)",
     comparison = TRUE,
     horizontal = TRUE)

RNdat1 <- with(wd,
              expand.grid(Allele_1=levels(Allele_1), WT_Background=levels(WT_Background)))

## this takes a little while
pp1 <- predict(m6, newdata=data.frame(RNdat1, Replicate="R1"), se.fit=TRUE)

RNdat1 <- with(pp1,
              data.frame(RNdat1,
                         wing_size_mm=fit,
                         wing_size_mm_lwr=fit-2*se.fit,
                         wing_size_mm_upr=fit+2*se.fit))

wRNdat <- with(wd,
               expand.grid(Allele_1=levels(Allele_1), WT_Background=levels(WT_Background)))

wRNdat <- with(pp,
               data.frame(wRNdat,
                          wing_size_mm=fit))

wRNdat

wideRNdat <- wRNdat %>% pivot_wider(names_from = "Allele_1", "WT_Background", values_from = "wing_size_mm")


wideRNdat


#wideRNdat <- wideRNdat[,3] -  wideRNdat[,2]

x <- wideRNdat[,2, drop=FALSE]
y <- cbind(x,x,x,x,x,x,x,x,x,x)


wideRNdat1 <- wideRNdat[,3:10] - y[,3:10]
wideRNdat1


wing_size_devaations <- data.frame(DGRP = wRNdat,wideRNdat)
wing_size_devaations

head(wing_size_deviations)

wing_size_deviations <- wing_size_deviations %>% 
  pivot_longer(c("bx[1]", "bx[2]", "bx[3]", "sd[29.1]", "sd[1]", "sd[E3]", "sd[ETX4]", "sd[58d]"),
               names_to = "Genotype", values_to = "Deviations")

library(colorspace)
theme_set(theme_bw())
gg2 <- (ggplot(RNdat1, aes(x=Allele_1, y=wing_size_mm,color=WT_Background))
        + geom_line(aes(group=WT_Background))
        + geom_point()
        + geom_ribbon(aes(ymin=wing_size_mm_lwr, ymax=wing_size_mm_upr, fill=WT_Background,
                          group=WT_Background),
                      colour=NA, alpha=0.05)
        + labs(y="Wing size (mm)", x= "Mutant Allele")
        + theme(legend.position = "")
)

print(gg2 
      + scale_colour_discrete_qualitative()
      + scale_fill_discrete_qualitative()
)



#### ah fuck ah fuck ah fuck ####


RNdatboth <- merge(RNdat, RNdat1)


ggplot(RNdatboth,
       aes(y = lev_stat, x = wing_size_mm, colour = Allele_1)) + geom_point()

ggplot(RNdatboth, 
       aes(y = lev_stat, x = wing_size_mm, colour = Allele_1)) +
  geom_point( size = 4, alpha = 0.6) +
  guides(shape = FALSE) +
  theme(legend.position = "top", 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  scale_colour_discrete(name = "Genotype",
                        breaks = c("WT", "sd291","sd1", "sdE3", "sdETX4", "sd58d"),
                        labels =
                          c(expression(paste(italic(sd)^"+")), 
                            expression(paste(italic(sd)^"29.1")),
                            expression(paste(italic(sd)^"1")),
                            expression(paste(italic(sd)^"E3")),
                            expression(paste(italic(sd)^"ETX4")),
                            expression(paste(italic(sd)^"58d")))) +
  labs(x = expression(paste(plain("wing area, mm")^2)), y = "Variability") +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.title.x = element_text(size = 16))

#not much change from original.
