################################################################################
## Calcuate BLUPS for GWAS in GAPIT
## Authors: Sarah Turner-Hissong & Makenzie Mabry
## Date: 25 May 2020
################################################################################
## 1. Set Working Directory
setwd("~/Box Sync/BoleraceaLeafScans/GWAS/")
library(ggplot2)
library(ggpubr)
library(forcats)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)

## 2. Read in dataset
qualdat = read.csv("/Users/mem2c2/OneDrive - University of Missouri/Computer/Projects/BoleraceaMorph/leaf_phenos_GWAS.csv", header=T)


## 3. Examine distribution of data, change for element of choice
hist1 <- hist(qualdat$circularity)
hist2 <- hist(qualdat$area)
hist3 <- hist(qualdat$convex_area)
hist4 <- hist(qualdat$solidity)
hist5 <- hist(qualdat$width)
hist6 <- hist(qualdat$length)
hist7 <- hist(qualdat$aspect_ratio)

## 4. check if data needs to be transformed
library(dplyr)
#install.packages("ggpubr")
library(ggpubr)
library(MASS)
library(lme4)

sample = as.factor(qualdat$sample)
rep = as.factor(qualdat$rep)
leafNum = as.factor(qualdat$leaf_num)
plantout = as.factor(qualdat$plantout)
morphotype = as.factor(qualdat$morphotype)

apply(qualdat[ ,6:16],2, shapiro.test)

#from the shapiro output, if pvalue is > 0.05 than it is normally distrubutied

## 6. Tranform data using box cox transfromation 
#install.packages("remotes")
#remotes::install_github("jonathon-love/car3")
library(car3)

#run linear mixed model
Circularity_model = lmer(qualdat$circularity ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
Area_model = lmer(qualdat$area ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
Convex_model = lmer(qualdat$convex_area ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
Solidity_model = lmer(qualdat$solidity ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
Width_model = lmer(qualdat$width ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
Length_model = lmer(qualdat$length ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
AspectRatio_model = lmer(qualdat$aspect_ratio ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))

PC1_model = lmer(qualdat$PC1 ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
PC2_model = lmer(qualdat$PC2 ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
PC3_model = lmer(qualdat$PC3 ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))
PC4_model = lmer(qualdat$PC4 ~ (1|sample) + (1|rep) + (1|leafNum) + (1|plantout))


## 7. Calcualting heritiability

# Extract variance components from models using the boxcox transformation
summCircularity <- summary(Circularity_model)
Circularity_sampleVar <- summCircularity$varcor$sample[1]
Circularity_repVar <- summCircularity$varcor$rep[1]
Circularity_leafNumVar <- summCircularity$varcor$leafNum[1]
Circularity_plantoutVar <- summCircularity$varcor$plantout[1]
Circularity_Residual <- attr(VarCorr(Circularity_model), "sc")^2

summArea <- summary(Area_model)
Area_sampleVar <- summArea$varcor$sample[1]
Area_repVar <- summArea$varcor$rep[1]
Area_leafNumVar <- summArea$varcor$leafNum[1]
Area_plantoutVar <- summArea$varcor$plantout[1]
Area_Residual <- attr(VarCorr(Area_model), "sc")^2

summConvex <- summary(Convex_model)
Convex_sampleVar <- summConvex$varcor$sample[1]
Convex_repVar <- summConvex$varcor$rep[1]
Convex_leafNumVar <- summConvex$varcor$leafNum[1]
Convex_plantoutVar <- summConvex$varcor$plantout[1]
Convex_Residual <- attr(VarCorr(Convex_model), "sc")^2

summSolidity <- summary(Solidity_model)
Solidity_sampleVar <- summSolidity$varcor$sample[1]
Solidity_repVar <- summSolidity$varcor$rep[1]
Solidity_leafNumVar <- summSolidity$varcor$leafNum[1]
Solidity_plantoutVar <- summSolidity$varcor$plantout[1]
Solidity_Residual <- attr(VarCorr(Solidity_model), "sc")^2

summWidth <- summary(Width_model)
Width_sampleVar <- summWidth$varcor$sample[1]
Width_repVar <- summWidth$varcor$rep[1]
Width_leafNumVar <- summWidth$varcor$leafNum[1]
Width_plantoutVar <- summWidth$varcor$plantout[1]
Width_Residual <- attr(VarCorr(Width_model), "sc")^2

summLength <- summary(Length_model)
Length_sampleVar <- summLength$varcor$sample[1]
Length_repVar <- summLength$varcor$rep[1]
Length_leafNumVar <- summLength$varcor$leafNum[1]
Length_plantoutVar <- summLength$varcor$plantout[1]
Length_Residual <- attr(VarCorr(Length_model), "sc")^2

summAspectRatio <- summary(AspectRatio_model)
AspectRatio_sampleVar <- summAspectRatio$varcor$sample[1]
AspectRatio_repVar <- summAspectRatio$varcor$rep[1]
AspectRatio_leafNumVar <- summAspectRatio$varcor$leafNum[1]
AspectRatio_plantoutVar <- summAspectRatio$varcor$plantout[1]
AspectRatio_Residual <- attr(VarCorr(AspectRatio_model), "sc")^2

#estimating heritability
#var(sample)/(var(sample)+(var(rep)/4)+(var(plantout)/2)+(var(residual)/8))

Circularity_heritability <- Circularity_sampleVar/(Circularity_sampleVar + (Circularity_repVar/4)+ (Circularity_leafNumVar/3) + (Circularity_plantoutVar/2) + (Circularity_Residual/24))
Area_heritability <- Area_sampleVar/(Area_sampleVar + (Area_repVar/4)+ (Area_leafNumVar/3) + (Area_plantoutVar/2) + (Area_Residual/24))
Convex_heritability <- Convex_sampleVar/(Convex_sampleVar + (Convex_repVar/4)+ (Convex_leafNumVar/3) + (Convex_plantoutVar/2) + (Convex_Residual/24))
Solidity_heritability <- Solidity_sampleVar/(Solidity_sampleVar + (Solidity_repVar/4)+ (Solidity_leafNumVar/3) + (Solidity_plantoutVar/2) + (Solidity_Residual/24))
Width_heritability <- Width_sampleVar/(Width_sampleVar + (Width_repVar/4)+ (Width_leafNumVar/3) + (Width_plantoutVar/2) + (Width_Residual/24))
Length_heritability <- Length_sampleVar/(Length_sampleVar + (Length_repVar/4)+ (Length_leafNumVar/3) + (Length_plantoutVar/2) + (Length_Residual/24))
AspectRatio_heritability <- AspectRatio_sampleVar/(AspectRatio_sampleVar + (AspectRatio_repVar/4)+ (AspectRatio_leafNumVar/3) + (AspectRatio_plantoutVar/2) + (AspectRatio_Residual/24))

Heritability <- data.frame(Circularity_heritability, 
                           Area_heritability,
                           Convex_heritability,
                           Solidity_heritability,
                           Width_heritability,
                           Length_heritability,
                           AspectRatio_heritability)

Measurement <- c("Circularity", "Area", "Convex", "Solidity", "Width","Length", "AspectRatio")
heritability <- c(Circularity_heritability, Area_heritability, Convex_heritability, Solidity_heritability, Width_heritability, Length_heritability, AspectRatio_heritability)

table1 <- data.frame(Measurement,heritability)
p <- ggplot(data = table1, aes(x=reorder(Measurement,heritability), y=heritability))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%0.2f", round(heritability, digits = 2))))

## 8. estimate BLUPS
Circularity_blup = ranef(Circularity_model)
Area_blup = ranef(Area_model)
Convex_blup = ranef(Convex_model)
Solidity_blup = ranef(Solidity_model)
Width_blup = ranef(Width_model)
Length_blup = ranef(Length_model)
AspectRatio_blup = ranef(AspectRatio_model)

PC1_blup = ranef(PC1_model)
PC2_blup = ranef(PC2_model)
PC3_blup = ranef(PC3_model)
PC4_blup = ranef(PC4_model)


# extract blup for line
Circularity_sampleblup = Circularity_blup$sample
Area_sampleblup = Area_blup$sample
Convex_sampleblup = Convex_blup$sample
Solidity_sampleblup = Solidity_blup$sample
Width_sampleblup = Width_blup$sample
Length_sampleblup = Length_blup$sample
AspectRatio_sampleblup = AspectRatio_blup$sample

PC1_sampleblup = PC1_blup$sample
PC2_sampleblup = PC2_blup$sample
PC3_sampleblup = PC3_blup$sample
PC4_sampleblup = PC4_blup$sample

#make datafram for all elements
Morpho_calcuated_BLUPS <- data.frame(Circularity_sampleblup = Circularity_sampleblup,
                                     Area_sampleblup = Area_sampleblup,
                                     Convex_sampleblup = Convex_sampleblup,
                                     Solidity_sampleblup = Solidity_sampleblup,
                                     Width_sampleblup = Width_sampleblup,
                                     Length_sampleblup = Length_sampleblup,
                                     AspectRatio_sampleblup = AspectRatio_sampleblup,
                                     PC1_sampleblup = PC1_sampleblup, 
                                     PC2_sampleblup = PC2_sampleblup, 
                                     PC3_sampleblup = PC3_sampleblup, 
                                     PC4_sampleblup = PC4_sampleblup)

colnames(Morpho_calcuated_BLUPS) <- c("Circularity_BLUP", "Area_BLUP", "Convex_BLUP", "Solidity_BLUP", "Width_BLUP", "Length_BLUP","AspectRatio_BLUP", "PC1_BLUP", "PC2_BLUP", "PC3_BLUP", "PC4_BLUP")


##use the following to write one file with all BLUPS for each element
write.csv(Morpho_calcuated_BLUPS, "Morpho_calcuated_BLUPS_no188.csv")

