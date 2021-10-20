---
title: "Chile Arsenic Study, epigenetic clocks"
output:
  html_document:
    toc: true
    toc_float: true
---

```{r, message=FALSE, warning=FALSE, echo = F}
# Required Packages

library(tidyverse)
library(minfi)
library(ggplot2)
library(stringr)
library(sqldf)
library(psych)
library(kableExtra)
library(ENmix)
library(gridExtra)
```

## PBMCs

### Normalization using funnorm
```{r, eval=FALSE, echo = TRUE}
# load RGset
load("pmbcs_RGset.rdata")

# Funnorm
funnorm = preprocessFunnorm(RGset)
betas_fun = minfi::getBeta(funnorm)

pheno = pData(funnorm)
pheno$sample_name = rownames(pheno)
# adding Sentrix ID for subsequent batch correction in downstream analysis
sentrix_id = as.factor(
  stringr::str_extract_all(pheno$sample_name, "^([0-9]{12})",
                           simplify = TRUE))
pheno = data.frame(pheno, sentrix_id)
# remove features without any variation
pheno = pheno[, apply(pheno, 2, function(x) length(unique(x)) > 1)]

# one sample with sex mismatch
table(pheno$sex, pheno$predictedSex)
     # F  M
  # F 19  0
  # M  1 20  
  
# save the processed betas and pheno data
save(betas_fun, pheno, file = "PBMC_betas_funnorm_pheno.RData")
```

### Process data for dnamage
```{r, eval=F, echo = TRUE}
# load data
# load('PBMC_betas_funnorm_pheno.RData')

# filtering beta values using CpGs from https://dnamage.genetics.ucla.edu/new
datMiniAnnotation=read.csv("datMiniAnnotation3.csv")
dim(datMiniAnnotation)
# 30084     7

match1=match(datMiniAnnotation[,1], rownames(betas_fun))

# filter probes
betas_fun_reduced = betas_fun[match1,]
dim(betas_fun_reduced)
# 30084    40

betas_fun_reduced = data.frame(ProbeID = datMiniAnnotation$Name, betas_fun_reduced)
dim(betas_fun_reduced)
# 30084    41

betas_fun_reduced[,1]=as.character(betas_fun_reduced[,1])

betas_fun_reduced = data.frame(betas_fun_reduced)

for (i in 2:dim(betas_fun_reduced)[[2]] ){betas_fun_reduced[,i]=
as.numeric(as.character(gsub(x= betas_fun_reduced[,i],pattern="\"",replacement=""))) }

colNames = colnames(betas_fun_reduced)[2:ncol(betas_fun_reduced)]
colNames = substr(colNames, 2, nchar(colNames))
colnames(betas_fun_reduced)[2:ncol(betas_fun_reduced)] = colNames

# number of missing probes
sum(is.na(betas_fun_reduced[,2]))
# 2558

# save filtered probes
write.table(betas_fun_reduced,"PBMC_Betas_funnorm_reduced_for_clock.csv", row.names=F, sep="," )

# pheno/annotation file
pheno$Female[pheno$sex == 'F'] = 1
pheno$Female[pheno$sex == 'M'] = 0

ageAnno = data.frame(ID = rownames(pheno), Age = pheno$age, Female = pheno$Female, Tissue = 'Blood PBMC')

all(ageAnno$ID == colnames(betas_fun_reduced)[-1])
# TRUE
identical(ageAnno$ID, colnames(betas_fun_reduced)[-1])
# TRUE

# save annotation file
write.table(ageAnno,"PBMCs_anno_for_clock.csv", row.names=F, sep="," )
```

### DNA methylation age calculator

https://dnamage.genetics.ucla.edu/submit

upload and select "normalize data" and "advanced analysis"

Six clocks: (1) DNAmAge = Horvath clock; (2) DNAmAgeHannum = Hannum clock; (3) DNAmPhenoAge = PhenoAge Clock; (4) DNAmAgeSkinBloodClock = Skin and blood clock; (5) DNAmGrimAge = GrimAge; DNAmTL = (6) DNAm telomere length

Acceleration measurements:(1) DNAmAge = AgeAccelerationResidual; (2) DNAmAgeHannum = AgeAccelerationResidualHannum; (3) DNAmPhenoAge = AgeAccelPheno; (4) DNAmAgeSkinBloodClock = DNAmAgeSkinBloodClockAdjAge; (5) DNAmGrimAge = AgeAccelGrim; (6) DNAmTL = DNAmTLAdjAge

IEAA = intrinsic epigenetic age (derived from Horvath's acceleration; adjusting for cell-types); EAA = extrinsic epigenetic age (derived from Hannum's clock) 

```{r, message=FALSE, warning=FALSE, eval=F}
# load result
PBMC_age = read.csv('PBMC_Betas_funnorm_reduced_for_clock.output.csv')
PBMC_age$SampleID = substr(PBMC_age$SampleID, 2, nchar(PBMC_age$SampleID))

all(PBMC_age$SampleID == pheno$sample_name)
# TRUE
identical(PBMC_age$SampleID, pheno$sample_name)
# TRUE

# merge result with pheno data
pheno = cbind(pheno, PBMC_age[,-1])

# predicted tissue
table(pheno$predictedTissue)
# Blood PBMC 
        # 40 
        
# drop duplicated column
pheno = pheno[,-94]
```

### Save pheno data with clocks
```{r, message=FALSE, warning=FALSE, eval=F}
write.csv(pheno, file = 'PBMC_fun_pheno_withAges.csv')

# drop sample with sex mismatch 
pheno = pheno[pheno$sex == pheno$predictedSex,]
dim(pheno)
# 39 154

table(pheno$sex, pheno$predictedSex)
     # F  M
  # F 19  0
  # M  0 20
  
# rename pheno 
PBMC_pheno = pheno
```

### Scatter plots of Horvath clocks and age
```{r, message=FALSE, warning=FALSE, eval=F, echo = F}
# function for calculating MAE
mae2 = function(df, clockVar, ageVar){
	mae = median(abs(df[[ageVar]] - df[[clockVar]])); return(mae)
}

# DNAmAge
cor.test(PBMC_pheno$age, PBMC_pheno$DNAmAge)
# r = 0.6172555, p = 2.845e-05
mae2(PBMC_pheno, 'DNAmAge', 'age')
# 0.8125248
scatter_DNAmAge_PBMC = ggplot(PBMC_pheno, aes(age, DNAmAge)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'Horvath DNAm age') + coord_cartesian(ylim = c(35,60), clip = 'off') + theme(axis.title.x = element_blank()) + labs(title = expression(paste('r = 0.62, ', italic('p'), ' < 0.001, ', 'MAE = 0.81'))) + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmAgeHannum
cor.test(PBMC_pheno$age, PBMC_pheno$DNAmAgeHannum)
# r = 0.7099599, p = 4.159e-07
mae(PBMC_pheno, 'DNAmAgeHannum', 'age')
# 9.387462
scatter_DNAmAgeHannum_PBMC = ggplot(PBMC_pheno, aes(age, DNAmAgeHannum)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'Hannum DNAm age') + coord_cartesian(ylim = c(23,56), clip = 'off') + theme(axis.title.x = element_blank()) + labs(title = expression(paste('r = 0.71, ', italic('p'), ' < 0.001, ', 'MAE = 9.39'))) + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmPhenoAge
cor.test(PBMC_pheno$age, PBMC_pheno$DNAmPhenoAge)
# r = 0.6219137, p = 2.376e-05
mae(PBMC_pheno, 'DNAmPhenoAge', 'age')
# 14.92139
scatter_DNAmPhenoAge_PBMC = ggplot(PBMC_pheno, aes(age, DNAmPhenoAge)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'DNAm PhenoAge') + coord_cartesian(ylim = c(15,50), clip = 'off') + theme(axis.title.x = element_blank()) + labs(title = expression(paste('r = 0.62, ', italic('p'), ' < 0.001, ', 'MAE = 14.92'))) + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmAgeSkinBloodClock
cor.test(PBMC_pheno$age, PBMC_pheno$DNAmAgeSkinBloodClock)
# r = 0.8571135, p = 3.302e-12
mae(PBMC_pheno, 'DNAmAgeSkinBloodClock', 'age')
# 1.106939
scatter_DNAmAgeSkinBloodClock_PBMC = ggplot(PBMC_pheno, aes(age, DNAmAgeSkinBloodClock)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'Skin and blood DNAm age') + coord_cartesian(ylim = c(35,64), clip = 'off') + labs(title = expression(paste('r = 0.86, ', italic('p'), ' < 0.001, ', 'MAE = 1.11')), x = 'Chronological age') + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmGrimAge
cor.test(PBMC_pheno$age, PBMC_pheno$DNAmGrimAge)
# r = 0.751369, p = 3.557e-08
mae(PBMC_pheno, 'DNAmGrimAge', 'age')
# 2.14003
scatter_DNAmGrimAge_PBMC = ggplot(PBMC_pheno, aes(age, DNAmGrimAge)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'DNAm GrimAge') + coord_cartesian(ylim = c(35,61), clip = 'off') + labs(title = expression(paste('r = 0.75, ', italic('p'), ' < 0.001, ', 'MAE = 2.14')), x = 'Chronological age') + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmTL
cor.test(PBMC_pheno$age, PBMC_pheno$DNAmTL)
# r = -0.6868437, p = 1.377e-06
scatter_DNAmTL_PBMC = ggplot(PBMC_pheno, aes(age, DNAmTL)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'DNAm telomere length') + coord_cartesian(ylim = c(6.2,7.5), clip = 'off') + labs(title = expression(paste('r = -0.69, ', italic('p'), ' < 0.001')), x = 'Chronological age') + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(scatter_DNAmAge_PBMC, scatter_DNAmAgeHannum_PBMC, scatter_DNAmPhenoAge_PBMC, scatter_DNAmAgeSkinBloodClock_PBMC, scatter_DNAmGrimAge_PBMC, scatter_DNAmTL_PBMC, ncol = 3)
quartz.save(file = 'PBMC_clocks_scatterPlots.png', type = 'png', dpi = 300)
```

```{r, echo = F, out.width = '100%', echo = F}
knitr::include_graphics("PBMC_clocks_scatterPlots.png")
```


### Linear models of associations between age acceleration and exposure

Adjusted for predicted sex and smoking

```{r, message=FALSE, warning=FALSE, eval=F, echo = F}
# function to run lm and create dataframe of output
lmAge = function(df, expose, ageVars, covars = NULL){
	dfAge = data.frame(matrix(ncol = 8, nrow = length(ageVars)))
	colnames(dfAge) = c('B_uadj', 'lower_unadj', 'upper_unadj', 'p_unadj', 'B_adj', 'lower_adj', 'upper_adj', 'p_adj')
	rownames(dfAge) = ageVars
	for (i in 1:length(ageVars)){
		mod_unadj = lm(as.formula(paste(ageVars[i], ' ~ ', expose, sep = '')), data = df)
		dfAge[i,1] = mod_unadj$coefficients[2]
		dfAge[i,2] = confint(mod_unadj)[2,1]
		dfAge[i,3] = confint(mod_unadj)[2,2]
		dfAge[i,4] = summary(mod_unadj)[4][[1]][2,4]
		mod_adj = lm(as.formula(paste(ageVars[i], ' ~ ', expose, ' + ', paste(covars, collapse = ' + '), sep = '')), data = df)
		dfAge[i,5] = mod_adj$coefficients[2]
		dfAge[i,6] = confint(mod_adj)[2,1]
		dfAge[i,7] = confint(mod_adj)[2,2]
		dfAge[i,8] = summary(mod_adj)[4][[1]][2,4]
	}
	return(dfAge)
}

PBMC_pheno$exposed = factor(PBMC_pheno$exposed, levels = c('no', 'yes'))
PBMC_pheno$predictedSex = factor(PBMC_pheno$predictedSex, levels = c('F', 'M'))
PBMC_pheno$smoking = factor(PBMC_pheno$smoking, levels = c('no', 'yes'))

PBMC_dfAge = lmAge(PBMC_pheno, 'exposed', c('AgeAccelerationResidual', 'AgeAccelerationResidualHannum', 'AgeAccelPheno', 'DNAmAgeSkinBloodClockAdjAge', 'AgeAccelGrim', 'DNAmTLAdjAge'), c('predictedSex', 'smoking'))
```

```{r, echo = F}
PBMC_dfAge %>% kable() %>% kable_styling(font_size = 14) 
```

```{r, message=FALSE, warning=FALSE, eval=F, echo = F}
label = c('Horvath', 'Hannum', 'PhenoAge', 'Skin and Blood', 'GrimAge')
mean = PBMC_dfAge$B_adj[-6]
lower = PBMC_dfAge$lower_adj[-6]
upper = PBMC_dfAge$upper_adj[-6]

df <- data.frame(label, mean, lower, upper)

# reverses the factor level ordering for labels after coord_flip()
df$label = factor(df$label, levels=rev(df$label))

fp_PBMCs = ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_point(size=2) + 
  geom_errorbar(size = 0.5, width = 0.3) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x='Epigenetic Age Acceleration (EAA)')+
  labs(y="Mean age acceleration (years), exposed vs. unexposed") +
  ggtitle("PBMCs")+
  theme_classic()        + # use a white background
  theme(plot.title = element_text(size=8.6, hjust = 0.5),panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.title.y=element_text(size=8.5))
quartz.save(file = 'PBMC_ageAccelAssoc_forestPlots.png', type = 'png', dpi = 300)
```

```{r, echo = F, out.width = '75%', echo = F}
knitr::include_graphics("PBMC_ageAccelAssoc_forestPlots.png")
```


## Buccal cells

### Normalization using funnorm
```{r, eval=FALSE, echo = TRUE}
# remove PBMC objects
rm(ageAnno, betas_fun, betas_fun_reduced, funnorm, pheno, RGset)

# load RGset
load("buccal_RGset.rdata")

# B14's detection p-values are significantly larger
RGset = RGset[,-26]
pheno = pData(RGset)
dim(pheno)
# 39 34

# Funnorm
funnorm = preprocessFunnorm(RGset)
betas_fun = minfi::getBeta(funnorm)

pheno = pData(funnorm)
pheno$sample_name = rownames(pheno)
# adding Sentrix ID for subsequent batch correction in downstream analysis
sentrix_id = as.factor(
  stringr::str_extract_all(pheno$sample_name, "^([0-9]{12})",
                           simplify = TRUE))
pheno = data.frame(pheno, sentrix_id)
# remove features without any variation
pheno = pheno[, apply(pheno, 2, function(x) length(unique(x)) > 1)]

# one sample with sex mismatch
table(pheno$sex, pheno$predictedSex)
     # F  M
  # F 19  0
  # M  1 20  d

# save the processed betas and pheno data
save(betas_fun, pheno, file = "buccalCells_betas_funnorm_pheno.RData")
```

### Process data for dnamage
```{r, eval=F, echo = TRUE}
# load data
# load('/Users/annebozack/Box/Chile_arsenic_epigenetic_aging/buccalCells_betas_funnorm_pheno.RData')

match1=match(datMiniAnnotation[,1], rownames(betas_fun))

# filter probes
betas_fun_reduced = betas_fun[match1,]
dim(betas_fun_reduced)
# 30084    39

betas_fun_reduced = data.frame(ProbeID = datMiniAnnotation$Name, betas_fun_reduced)
dim(betas_fun_reduced)
# 30084    41

betas_fun_reduced[,1]=as.character(betas_fun_reduced[,1])

betas_fun_reduced = data.frame(betas_fun_reduced)

for (i in 2:dim(betas_fun_reduced)[[2]] ){betas_fun_reduced[,i]=
as.numeric(as.character(gsub(x= betas_fun_reduced[,i],pattern="\"",replacement=""))) }

colNames = colnames(betas_fun_reduced)[2:ncol(betas_fun_reduced)]
colNames = substr(colNames, 2, nchar(colNames))
colnames(betas_fun_reduced)[2:ncol(betas_fun_reduced)] = colNames

# number of missing probes
sum(is.na(betas_fun_reduced[,2]))
# 2558

# save filtered probes
write.table(betas_fun_reduced,"buccalCells_Betas_funnorm_reduced_for_clock.csv", row.names=F, sep="," )

# pheno/annotation file
pheno$Female[pheno$sex == 'F'] = 1
pheno$Female[pheno$sex == 'M'] = 0

ageAnno = data.frame(ID = rownames(pheno), Age = pheno$age, Female = pheno$Female, Tissue = 'Blood PBMC')

all(ageAnno$ID == colnames(betas_fun_reduced)[-1])
# TRUE
identical(ageAnno$ID, colnames(betas_fun_reduced)[-1])
# TRUE

# save annotation file
write.table(ageAnno,"buccalCells_anno_for_clock.csv", row.names=F, sep="," )
```

### DNA methylation age calculator

```{r, message=FALSE, warning=FALSE, eval=F}
# load result
buccalCell_age = read.csv('buccalCells_Betas_funnorm_reduced_for_clock.output.csv')
buccalCell_age$SampleID = substr(buccalCell_age$SampleID, 2, nchar(buccalCell_age$SampleID))

all(buccalCell_age$SampleID == pheno$sample_name)
# TRUE
identical(buccalCell_age$SampleID, pheno$sample_name)
# TRUE

# merge result with pheno data
pheno = cbind(pheno, buccalCell_age[,-1])

# predicted tissue
table(pheno$predictedTissue)
# Buccal Saliva 
    # 22     17 
        
# drop duplicated column
pheno = pheno[,-92]
```

### Save pheno data with clocks
```{r, message=FALSE, warning=FALSE, eval=F}
write.csv(pheno, file = 'buccalCell_fun_pheno_withAges.csv')

# drop sample with sex mismatch 
pheno = pheno[pheno$sex == pheno$predictedSex,]
dim(pheno)
# 38 142

table(pheno$sex, pheno$predictedSex)
     # F  M
  # F 19  0
  # M  0 19
  
# rename pheno 
buccalCell_pheno = pheno
```

### Scatter plots of Horvath clocks and age
```{r, message=FALSE, warning=FALSE, eval=F, echo = F}
# DNAmAge
cor.test(buccalCell_pheno$age, buccalCell_pheno$DNAmAge)
# r = 0.4612412, p = 0.003562
mae(buccalCell_pheno, 'DNAmAge', 'age')
# 8.055889
scatter_DNAmAge_buccalCell = ggplot(buccalCell_pheno, aes(age, DNAmAge)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'Horvath DNAm age') + coord_cartesian(ylim = c(21,52), clip = 'off') + theme(axis.title.x = element_blank()) + labs(title = expression(paste('r = 0.46, ', italic('p'), ' = 0.004, ', 'MAE = 8.06'))) + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmAgeHannum
cor.test(buccalCell_pheno$age, buccalCell_pheno$DNAmAgeHannum)
# r = 0.4971326, p = 0.001497
mae(buccalCell_pheno, 'DNAmAgeHannum', 'age')
# 0.8818825
scatter_DNAmAgeHannum_buccalCell = ggplot(buccalCell_pheno, aes(age, DNAmAgeHannum)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'Hannum DNAm age') + coord_cartesian(ylim = c(36,62), clip = 'off') + theme(axis.title.x = element_blank()) + labs(title = expression(paste('r = 0.50, ', italic('p'), ' = 0.001, ', 'MAE = 0.88'))) + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmPhenoAge
cor.test(buccalCell_pheno$age, buccalCell_pheno$DNAmPhenoAge)
# r = 0.3858514, p = 0.01673
mae(buccalCell_pheno, 'DNAmPhenoAge', 'age')
# 2.052611
scatter_DNAmPhenoAge_buccalCell = ggplot(buccalCell_pheno, aes(age, DNAmPhenoAge)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'DNAm PhenoAge') + coord_cartesian(ylim = c(19,74), clip = 'off') + theme(axis.title.x = element_blank()) + labs(title = expression(paste('r = 0.39, ', italic('p'), ' = 0.017, ', 'MAE = 2.05'))) + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmAgeSkinBloodClock
cor.test(buccalCell_pheno$age, buccalCell_pheno$DNAmAgeSkinBloodClock)
# r = 0.7206731, p = 3.354e-07
mae(buccalCell_pheno, 'DNAmAgeSkinBloodClock', 'age')
# 1.687192
scatter_DNAmAgeSkinBloodClock_buccalCell = ggplot(buccalCell_pheno, aes(age, DNAmAgeSkinBloodClock)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'Skin and blood DNAm age') + coord_cartesian(ylim = c(34,63), clip = 'off') + labs(title = expression(paste('r = 0.72, ', italic('p'), ' < 0.001, ', 'MAE = 1.69')), x = 'Chronological age') + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmGrimAge
cor.test(buccalCell_pheno$age, buccalCell_pheno$DNAmGrimAge)
# r = 0.6625884, p = 5.847e-06
mae(buccalCell_pheno, 'DNAmGrimAge', 'age')
# 11.11323
scatter_DNAmGrimAge_buccalCell = ggplot(buccalCell_pheno, aes(age, DNAmGrimAge)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'DNAm GrimAge') + coord_cartesian(ylim = c(47,73), clip = 'off') + labs(title = expression(paste('r = 0.66, ', italic('p'), ' < 0.001, ', 'MAE = 11.11')), x = 'Chronological age') + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

# DNAmTL
cor.test(buccalCell_pheno$age, buccalCell_pheno$DNAmTL)
# r = -0.3875084, p = 0.01623
scatter_DNAmTL_buccalCell = ggplot(buccalCell_pheno, aes(age, DNAmTL)) + geom_smooth(method = 'lm', alpha = 0.2, size = 0.5, color = 'black') + geom_point(size = 0.75) + theme_classic() + labs(x = 'Chronological age (years)', y = 'DNAm telomere length') + coord_cartesian(ylim = c(5.25, 7.5), clip = 'off') + labs(title = expression(paste('r = -0.39, ', italic('p'), ' = 0.016')), x = 'Chronological age') + theme(plot.title = element_text(size=8.5)) + theme(text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(scatter_DNAmAge_buccalCell, scatter_DNAmAgeHannum_buccalCell, scatter_DNAmPhenoAge_buccalCell, scatter_DNAmAgeSkinBloodClock_buccalCell, scatter_DNAmGrimAge_buccalCell, scatter_DNAmTL_buccalCell, ncol = 3)
quartz.save(file = 'buccalCell_clocks_scatterPlots.png', type = 'png', dpi = 300)
```

```{r, echo = F, out.width = '100%', echo = F}
knitr::include_graphics("PBMC_clocks_scatterPlots.png")
```


### Linear models of associations between age acceleration and exposure

Adjusted for predicted sex and smoking

```{r, message=FALSE, warning=FALSE, eval=F, echo = F}
buccalCell_pheno$exposed = factor(buccalCell_pheno$exposed, levels = c('no', 'yes'))
buccalCell_pheno$predictedSex = factor(buccalCell_pheno$predictedSex, levels = c('F', 'M'))
buccalCell_pheno$smoking = factor(buccalCell_pheno$smoking, levels = c('no', 'yes'))

buccalCell_dfAge = lmAge(buccalCell_pheno, 'exposed', c('AgeAccelerationResidual', 'AgeAccelerationResidualHannum', 'AgeAccelPheno', 'DNAmAgeSkinBloodClockAdjAge', 'AgeAccelGrim', 'DNAmTLAdjAge'), c('predictedSex', 'smoking'))
```

```{r, echo = F}
buccalCell_dfAge %>% kable() %>% kable_styling(font_size = 14) 
```

```{r, message=FALSE, warning=FALSE, eval=F, echo = F}
label = c('Horvath', 'Hannum', 'PhenoAge', 'Skin and Blood', 'GrimAge')
mean = buccalCell_dfAge$B_adj[-6]
lower = buccalCell_dfAge$lower_adj[-6]
upper = buccalCell_dfAge$upper_adj[-6]

df <- data.frame(label, mean, lower, upper)

# reverses the factor level ordering for labels after coord_flip()
df$label = factor(df$label, levels=rev(df$label))

fp_buccalCell = ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_point(size=2) + 
  geom_errorbar(size = 0.5, width = 0.3) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x='Epigenetic Age Acceleration (EAA)')+
  labs(y="Mean age acceleration (years), exposed vs. unexposed") +
  ggtitle("Buccal Cells")+
  theme_classic()        + # use a white background
  theme(plot.title = element_text(size=8.6, hjust = 0.5),panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.title.y=element_text(size=8.5))
quartz.save(file = 'buccalCell_ageAccelAssoc_forestPlots.png', type = 'png', dpi = 300)
```

```{r, echo = F, out.width = '75%', echo = F}
knitr::include_graphics("buccalCell_ageAccelAssoc_forestPlots.png")
```


## Combining PBMC and buccal cell plots
```{r, message=FALSE, warning=FALSE, eval=F}
# scatter plots 
titles1 = ggplot(data.frame(x = c(1:3), y = 1)) + annotate('text', label = 'Horvath', x = 1.1, y = 1, size = 3.5) + annotate('text', label = 'Hannum', x = 2, y = 1, size = 3.5) + annotate('text', label = 'PhenoAge', x = 2.92, y = 1, size = 3.5) + xlim(c(0.75,3.25)) + theme_void() 
titles2 = ggplot(data.frame(x = c(1:3), y = 1)) + annotate('text', label = 'Skin and blood', x = 1.18, y = 1, size = 3.5) + annotate('text', label = 'GrimAge', x = 2, y = 1, size = 3.5) + annotate('text', label = 'DNAmTL', x = 2.95, y = 1, size = 3.5)+ xlim(c(0.75,3.25)) + theme_void()

header1 = ggplot(data.frame(x = c(1:3), y = 1)) + annotate('text', label = 'PBMCs', x = 2, y = 1, size = 5) + xlim(c(0.75,3.25)) + theme_void() 
header2 = ggplot(data.frame(x = c(1:3), y = 1)) + annotate('text', label = 'Buccal cells', x = 2, y = 1, size = 5) + xlim(c(0.75,3.25)) + theme_void() 

scatter_DNAmAge_PBMC2 = scatter_DNAmAge_PBMC + theme(axis.title.x = element_blank()) + labs(y = 'DNAm age')
scatter_DNAmAgeHannum_PBMC2 = scatter_DNAmAgeHannum_PBMC + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
scatter_DNAmPhenoAge_PBMC2 = scatter_DNAmPhenoAge_PBMC + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
scatter_DNAmAgeSkinBloodClock_PBMC2 = scatter_DNAmAgeSkinBloodClock_PBMC + theme(axis.title.x = element_blank()) + labs(y = 'DNAm age')
scatter_DNAmGrimAge_PBMC2 = scatter_DNAmGrimAge_PBMC + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
scatter_DNAmTL_PBMC2 = scatter_DNAmTL_PBMC + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

scatter_DNAmAge_buccalCell2 = scatter_DNAmAge_buccalCell + labs(y = 'DNAm age') + theme(axis.title.x = element_blank())
scatter_DNAmAgeHannum_buccalCell2 = scatter_DNAmAgeHannum_buccalCell + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
scatter_DNAmPhenoAge_buccalCell2 = scatter_DNAmPhenoAge_buccalCell + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
scatter_DNAmAgeSkinBloodClock_buccalCell2 = scatter_DNAmAgeSkinBloodClock_buccalCell + labs(x = 'Chronological age', y = 'DNAm age')
scatter_DNAmGrimAge_buccalCell2 = scatter_DNAmGrimAge_buccalCell + theme(axis.title.y = element_blank()) + labs(x = 'Chronological age')
scatter_DNAmTL_buccalCell2 = scatter_DNAmTL_buccalCell + theme(axis.title.y = element_blank()) + labs(x = 'Chronological age')

g1 <- arrangeGrob(scatter_DNAmAge_PBMC2, scatter_DNAmAgeHannum_PBMC2, scatter_DNAmPhenoAge_PBMC2, nrow = 1)
g2 <- arrangeGrob(scatter_DNAmAgeSkinBloodClock_PBMC2, scatter_DNAmGrimAge_PBMC2, scatter_DNAmTL_PBMC2, nrow = 1)
g3 <- arrangeGrob(scatter_DNAmAge_buccalCell2,  scatter_DNAmAgeHannum_buccalCell2, scatter_DNAmPhenoAge_buccalCell2, nrow = 1)
g4 <- arrangeGrob(scatter_DNAmAgeSkinBloodClock_buccalCell2, scatter_DNAmGrimAge_buccalCell2, scatter_DNAmTL_buccalCell2, nrow = 1)

grid.arrange(header1,titles1,g1,titles2,g2,blank,header2,titles1,g3,titles2,g4,nrow=11, heights=unit(c(0.6,0.4,5.25,0.4,5.25,0.2,0.6,0.4,5.25,0.4,5.4), c("cm")))
quartz.save(file = 'all_clocks_scatterPlots_long.png', type = 'png', dpi = 300)

# forest plots
fp_buccalCell2 = fp_buccalCell + theme(axis.title.y = element_blank())
grid.arrange(fp_PBMCs, fp_buccalCell2, ncol = 2)
quartz.save(file = 'all_ageAccelAssoc_forestPlots.png', type = 'png', dpi = 300)
```


