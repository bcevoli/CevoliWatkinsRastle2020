# STUDY 1 - EXP. 2 Semantic Diversity & Ambiguity


###############
# Preparation #
###############

#Clean workspace
rm(list=ls(all=T))

#Set working directory
wd = paste(dirname(dirname(getwd())),"/", sep="")
setwd(wd)

#Import used libraries
library(MASS)
library(grid)
library(magrittr)
library(data.table)
library(dplyr)
library(lmerTest)
library(MuMIn)
library(lattice)
library(lme4)
library(sjPlot)
library(ggplot2)
library(gridExtra)
library(psych)
library(languageR)
library(car)
library(dplyr)
library(sjmisc)
library(emmeans)
library(arules)
library(tidyr)
library(yarrr)
library(car)
library(plyr)
library("readxl")
library(ggsignif)

inverse <- function(x) {-1000/x};

# Read in & simple preprocessing
##############################################################################################################################################
# Read in Lexical Variables
diversity = read.csv("S1_SemanticDiversity_LSA/Data/SemanticDiversity/diversity_WideFormat.csv")
diversity = diversity[,c("Word","SemD_P2_W","Length")]
setnames(diversity, c("SemD_P2_W"), c("SemD"))
#Look up the number of words (28,555)
length(unique(diversity[!(is.na(diversity$SemD)),]$Word))
#Add frequency estimates from SUBTLEX UK (Van Heuven et al., 2013)
subtlex = read.csv("Resources/Norms/Databases/SUBTLEX/SUBTLEX-UK.csv")
setnames(subtlex, c("Spelling","LogFreq.Zipf."), c("Word","Freq"))
lexicalVars = merge(diversity, subtlex[,c("Word","Freq")], by="Word")
str(lexicalVars)
summary(lexicalVars)

# Ambiguity Stimuli Dataset
ambiguityStimuli = read_excel("S1_SemanticDiversity_LSA/Exp2_Ambiguity&SemD/AmbiguityStimuliCollection.xlsx", sheet = "Stimuli")
colnames(ambiguityStimuli)

# BLP data
blp.ld.data = read.table("Resources/Data/BritishLexiconProject/blp-trials.txt/blp-trials.txt", sep="\t", header=TRUE)
blp.ld.data$lexicality = mapvalues(blp.ld.data$lexicality, from=c("N","W"), to=c(0,1))
blp.ld.data$lexicality  = as.numeric(as.character(blp.ld.data$lexicality))
#Add previous trial data
blp.ld.data$pLexicality = as.numeric(unlist(append(list(NA), head(blp.ld.data$lexicality, -1))))
#Select words only
blp.ld.data = subset(blp.ld.data, lexicality==1) 
#Subset only columns of interest (rt col: outliers and inaccurate responses set to NA)
blp.ld.data = blp.ld.data[c("spelling","participant","block","trial","accuracy","rt.raw","previous.rt","previous.accuracy","pLexicality")]
#Rename columns
setnames(blp.ld.data, c("spelling","participant","block","trial","accuracy","rt.raw", "previous.rt", "previous.accuracy"), c("Word","SubjID","Block","Trial","ACC","RT","pRT","pACC"))
#Look up the number of words (28,730)
length(unique(blp.ld.data$Word)) 
summary(blp.ld.data)

#ELP data
elp.ld.data = read.table("Resources/Data/EnglishLexiconProject/ELPDecisionData.csv", sep=",", header=TRUE)
#Rename columns
setnames(elp.ld.data, c("Stimulus","Participant", "Accuracy"), c("Word","SubjID", "ACC"))
#Set SubjID as factor
elp.ld.data$SubjID = as.factor(elp.ld.data$SubjID)
#Add previous trial data
elp.ld.data$pACC = as.numeric(unlist(append(list(NA), head(elp.ld.data$ACC, -1))))
elp.ld.data$pRT = as.numeric(unlist(append(list(NA), head(elp.ld.data$RT, -1))))
elp.ld.data$pLexicality = as.numeric(unlist(append(list(NA), head(elp.ld.data$Type, -1))))
#Select words only
elp.ld.data = subset(elp.ld.data, Type==1) 
#Look up the number of words (40,481)
length(unique(elp.ld.data$Word)) 
summary(elp.ld.data)
#remove strage negative RT
elp.ld.data = subset(elp.ld.data, RT>0)
#Look up the number of words (40,481)
length(unique(elp.ld.data$Word)) 
summary(elp.ld.data)

##############################################################################################################################################

# Replication & Simulation Analysis of Rodd et al. (2001) Experiment 1 (Visual Lexical Decision - Regression Design)
#####################################################################################################################

# Get stimuli & properties
###############################
roddStimuli = subset(ambiguityStimuli, Authors=="Rodd et al." & Experiment==1)[c("Word", "Meanings")]
roddStimuli = roddStimuli[c("Word", "Meanings")]
roddStimuli$Meanings = as.factor(roddStimuli$Meanings)
roddStimuli$Word = as.factor(roddStimuli$Word)
summary(roddStimuli)
###############################


# BLP - Reaction Time Analysis #
################################

# Get behavioural data from BLP
rodd.blp = merge(roddStimuli, blp.ld.data, by="Word", all.x = TRUE)
rodd.blp = merge(rodd.blp, subtlex[,c("Word","Freq")], by="Word", all.x = TRUE)
rodd.blp$Length = nchar(as.character(rodd.blp$Word))
#number of meanings/senses
wordsmythDict = read_excel("Resources/Norms/Dictionaries/WordsmythDictionary/Wordsmyth_Words_NumMeaning_NumSenses_PartsOfSpeechFreq.xls")
setnames(wordsmythDict, c("word", "NumMeanings", "NumSenses"), c("Word", "nMeanings_WMD", "nSenses_WMD"))
wordsmythDict = wordsmythDict[,c("Word","nMeanings_WMD","nSenses_WMD")]
rodd.blp = merge(rodd.blp, wordsmythDict, by="Word", all.x = TRUE)
#lexical neighborhood
blpStimuli = read.table("Resources/Data/BritishLexiconProject/blp-stimuli.txt/blp-stimuli.txt", sep="\t", header=TRUE)
setnames(blpStimuli, c("spelling", "coltheart.N"), c("Word", "coltN"))
rodd.blp = merge(rodd.blp, blpStimuli[c("Word","coltN")], by="Word", all.x = TRUE)
#Concreteness
concr = read.table("Resources/Norms/Measures/Concreteness_Brysbaert2014.csv", sep=",", header=TRUE)
rodd.blp = merge(rodd.blp, concr[c("Word","Conc.M")], by="Word", all.x = TRUE)
summary(rodd.blp)
str(rodd.blp)
rodd.blp$SubjID = as.factor(rodd.blp$SubjID)
length(unique(rodd.blp$Word)) #182 w

#Check accuracy
roddData.acc.subj = aggregate(ACC ~ SubjID, FUN=mean, rodd.blp)
summary(roddData.acc.subj) #no subjs with an error rate of more than 10%

#Check subject-level data
roddData.acc.subj = aggregate(ACC ~ SubjID, FUN=mean, rodd.blp)
summary(roddData.acc.subj) #no subjs with an error rate of more than 10%
roddData.rt.subj = aggregate(RT ~ SubjID, FUN=mean, rodd.blp)
summary(roddData.rt.subj) # no subjs with mean RT more than 1000

# Cleaning procedure 
#Remove incorrect responses
rodd.blp.rt = subset(rodd.blp, ACC==1)
#Remove responses longer than 1200ms and shorter than
rodd.blp.rt = subset(rodd.blp.rt, RT<1200)
summary(rodd.blp.rt)


#Check normality
hist(rodd.blp.rt$RT)
qqnorm(rodd.blp.rt$RT)
qqline(rodd.blp.rt$RT)
boxcox(rodd.blp.rt$RT ~ 1) #lambda ~ -1, inverse RT is recommended
hist(-1000/rodd.blp.rt$RT)
qqnorm(-1000/rodd.blp.rt$RT)
qqline(-1000/rodd.blp.rt$RT)

rodd.blp.rt$invRT = -1000/rodd.blp.rt$RT
summary(rodd.blp.rt)
rodd.blp.rt$Meanings <- relevel(rodd.blp.rt$Meanings, ref = "One")
rodd.exp1.blp.rt.m = lm(invRT ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
                   rodd.blp.rt)
summary(rodd.exp1.blp.rt.m)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.48109 -0.25341 -0.01956  0.23258  1.22741 
# 
# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)        -1.937325   0.009593 -201.948  < 2e-16 ***
#   MeaningsMany        0.023971   0.012143    1.974 0.048417 *  
#   scale(nSenses_WMD) -0.013501   0.006056   -2.229 0.025832 *  
#   scale(coltN)        0.015710   0.006415    2.449 0.014350 *  
#   scale(Freq)        -0.057227   0.005478  -10.447  < 2e-16 ***
#   scale(Length)       0.009910   0.006474    1.531 0.125868    
# scale(Conc.M)      -0.018483   0.005068   -3.647 0.000267 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.3899 on 6407 degrees of freedom
# (151 observations deleted due to missingness)
# Multiple R-squared:  0.02576,	Adjusted R-squared:  0.02484 
# F-statistic: 28.23 on 6 and 6407 DF,  p-value: < 2.2e-16

m.plot = lm(RT ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
                              rodd.blp.rt)
summary(m.plot)

library(effects)
rodd.exp1.blp.rt.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=m.plot, xlevels=list(nSenses_WMD=2)))
#rodd.exp1.blp.rt.senses.e[,c('fit','lower','upper')] <- inverse(rodd.exp1.blp.rt.senses.e[,c('fit','lower','upper')])
rodd.exp1.blp.rt.meanings.e = as.data.frame(effect(term="Meanings", mod=m.plot))
#rodd.exp1.blp.rt.meanings.e[,c('fit','lower','upper')] <- inverse(rodd.exp1.blp.rt.meanings.e[,c('fit','lower','upper')])
rodd.exp1.blp.rt.meanings.e$Meanings <- relevel(rodd.exp1.blp.rt.meanings.e$Meanings, ref = "One")

cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_BLP_RT_effects.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=4.5, width=9);
s = ggplot(rodd.exp1.blp.rt.senses.e, aes(x=nSenses_WMD, y=fit)) + 
  ylim(485,600)+
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=fit-se,ymax=fit+se),alpha=0.3, colour = NA) + 
  labs(x= "Number of Senses", y="Mean Reaction Time (ms)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") 
m= ggplot(rodd.exp1.blp.rt.meanings.e, aes(x=Meanings, y=fit, fill=Meanings)) + 
  #ylim(450,600) +
  coord_cartesian(ylim=c(485,600)) +
  geom_bar(stat="identity", position="dodge", color="white", alpha = 0.5) + 
  geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.1, position=position_dodge(width=0.9))+
  labs(x= "Meanings", y="Mean Reaction Time (ms)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") 
grid.arrange(s, m, nrow=1, ncol=2,
             top=textGrob("BLP Lexical Decision", gp=gpar(fontsize=20)))
dev.off()

rodd.blp.rt$Senses = discretize(rodd.blp.rt$nSenses_WMD, breaks = 2)
rodd.blp.rt$Senses = mapvalues(rodd.blp.rt$Senses, from = levels(rodd.blp.rt$Senses), to = c("Few", "Many"))
rodd.rt.item <- aggregate(RT ~ Word + Meanings + Senses + Freq + Length, FUN=mean, rodd.blp.rt)
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_BLP_RT.png.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=6, width=8);
par(mfrow=c(1,2))
pirateplot(formula = RT ~ Senses,
           data = rodd.rt.item,
           ylab = "Mean Reaction Time (ms)",
           xlab = "Senses",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
pirateplot(formula = RT ~ Meanings,
           data = rodd.rt.item,
           ylab = "Mean Reaction Time (ms)",
           xlab = "Meanings",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
dev.off();
################################

# ELP - Reaction Time Analysis #
################################

# Get behavioural data from ELP
rodd.elp = merge(roddStimuli, elp.ld.data, by="Word", all.x = TRUE)
rodd.elp = merge(rodd.elp, subtlex[,c("Word","Freq")], by="Word", all.x = TRUE)
rodd.elp$Length = nchar(as.character(rodd.elp$Word))
#number of meanings/senses
rodd.elp = merge(rodd.elp, wordsmythDict, by="Word", all.x = TRUE)
#lexical neighborhood
rodd.elp = merge(rodd.elp, blpStimuli[c("Word","coltN")], by="Word", all.x = TRUE)
#Concreteness
rodd.elp = merge(rodd.elp, concr[c("Word","Conc.M")], by="Word", all.x = TRUE)
summary(rodd.elp)
str(rodd.elp)
rodd.elp$SubjID = as.factor(rodd.elp$SubjID)
length(unique(rodd.elp$Word)) #182 w

#Check subject-level data
roddData.acc.subj = aggregate(ACC ~ SubjID, FUN=mean, rodd.elp)
summary(roddData.acc.subj) #no subjs with an error rate of more than 10%
#exclude subjects with mean RT more than 1000
subjToExclude = subset(roddData.rt.subj, RT>1000)$SubjID
rodd.elp = subset(rodd.elp, !(rodd.elp$SubjID %in% subjToExclude))

# Cleaning procedure 
#Remove incorrect responses
rodd.elp.rt = subset(rodd.elp, ACC==1)
#Remove responses longer than 1200ms and shorter than
rodd.elp.rt = subset(rodd.elp.rt, RT>200 & RT<1200)
summary(rodd.elp.rt)

#Check normality
hist(rodd.elp.rt$RT)
qqnorm(rodd.elp.rt$RT)
qqline(rodd.elp.rt$RT)
boxcox(rodd.elp.rt$RT ~ 1) #lambda ~ -1, inverse RT is recommended
hist(-1000/rodd.elp.rt$RT)
qqnorm(-1000/rodd.elp.rt$RT)
qqline(-1000/rodd.elp.rt$RT)

rodd.elp.rt$invRT = (-1000/rodd.elp.rt$RT)
summary(rodd.elp.rt)

rodd.elp.rt$Meanings <- relevel(rodd.elp.rt$Meanings, ref = "One")
rodd.elp.rt.m = lm(invRT ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
                   rodd.elp.rt)
summary(rodd.elp.rt.m)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.98017 -0.24985 -0.00569  0.26002  1.14970 
# 
# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)        -1.813385   0.011470 -158.096  < 2e-16 ***
#   MeaningsMany        0.036284   0.014510    2.501 0.012430 *  
#   scale(nSenses_WMD) -0.014073   0.007165   -1.964 0.049562 *  
#   scale(coltN)        0.013399   0.007641    1.754 0.079572 .  
# scale(Freq)        -0.056975   0.006553   -8.694  < 2e-16 ***
#   scale(Length)       0.028424   0.008110    3.505 0.000461 ***
#   scale(Conc.M)      -0.014086   0.006029   -2.336 0.019512 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.4136 on 5083 degrees of freedom
# (278 observations deleted due to missingness)
# Multiple R-squared:  0.02892,	Adjusted R-squared:  0.02777 
# F-statistic: 25.23 on 6 and 5083 DF,  p-value: < 2.2e-16


m.plot = lm(RT ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
            rodd.elp.rt)
summary(m.plot)

library(effects)
rodd.exp1.elp.rt.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=m.plot, xlevels=list(nSenses_WMD=2)))
#rodd.exp1.elp.rt.senses.e[,c('fit','lower','upper')] <- inverse(rodd.exp1.elp.rt.senses.e[,c('fit','lower','upper')])
rodd.exp1.elp.rt.meanings.e = as.data.frame(effect(term="Meanings", mod=m.plot))
#rodd.exp1.elp.rt.meanings.e[,c('fit','lower','upper')] <- inverse(rodd.exp1.elp.rt.meanings.e[,c('fit','lower','upper')])
rodd.exp1.elp.rt.meanings.e$Meanings <- relevel(rodd.exp1.elp.rt.meanings.e$Meanings, ref = "One")

cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_ELP_RT_effects.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=4.5, width=9);
s = ggplot(rodd.exp1.elp.rt.senses.e, aes(x=nSenses_WMD, y=fit)) + 
  ylim(550,650)+
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=fit-se,ymax=fit+se),alpha=0.3, colour = NA) + 
  labs(x= "Number of Senses", y="Mean Reaction Time (ms)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") 
m= ggplot(rodd.exp1.elp.rt.meanings.e, aes(x=Meanings, y=fit, fill=Meanings)) + 
  #ylim(450,600) +
  coord_cartesian(ylim=c(550,650)) +
  geom_bar(stat="identity", position="dodge", color="white", alpha = 0.5) + 
  geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.1, position=position_dodge(width=0.9))+
  labs(x= "Meanings", y="Mean Reaction Time (ms)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") 
grid.arrange(s, m, nrow=1, ncol=2,
             top=textGrob("ELP Lexical Decision", gp=gpar(fontsize=20)))
dev.off()


rodd.elp.rt$Senses = discretize(rodd.elp.rt$nSenses_WMD, breaks = 2)
rodd.elp.rt$Senses = mapvalues(rodd.elp.rt$Senses, from = levels(rodd.elp.rt$Senses), to = c("Few", "Many"))
rodd.rt.item <- aggregate(RT ~ Word + Meanings + Senses + Freq + Length, FUN=mean, rodd.elp.rt)
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_ELP_RT.png.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=6, width=8);
par(mfrow=c(1,2))
pirateplot(formula = RT ~ Senses,
           data = rodd.rt.item,
           ylab = "Mean Reaction Time (ms)",
           xlab = "Senses",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
pirateplot(formula = RT ~ Meanings,
           data = rodd.rt.item,
           ylab = "Mean Reaction Time (ms)",
           xlab = "Meanings",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")

dev.off();

################################

# Semantic Diversity Simulation #
#################################

# Get data 
rodd.exp1.semd = merge(roddStimuli, lexicalVars, by="Word", all.x = TRUE)
#number of meanings/senses
rodd.exp1.semd = merge(rodd.exp1.semd, wordsmythDict, by="Word", all.x = TRUE)
#lexical neighborhood
rodd.exp1.semd = merge(rodd.exp1.semd, blpStimuli[c("Word","coltN")], by="Word", all.x = TRUE)
#Concreteness
rodd.exp1.semd = merge(rodd.exp1.semd, concr[c("Word","Conc.M")], by="Word", all.x = TRUE)
summary(rodd.exp1.semd)

rodd.exp1.semd$Meanings <- relevel(rodd.exp1.semd$Meanings, ref = "One")
rodd.exp1.semd.m = lm(SemD ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
                      rodd.exp1.semd)
summary(rodd.exp1.semd.m)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.58984 -0.17932  0.00601  0.17353  0.67135 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         2.224456   0.038371  57.973  < 2e-16 ***
#   MeaningsMany        0.051402   0.048367   1.063 0.289471    
# scale(nSenses_WMD)  0.022481   0.023772   0.946 0.345701    
# scale(coltN)        0.008564   0.025571   0.335 0.738132    
# scale(Freq)         0.074274   0.022253   3.338 0.001047 ** 
#   scale(Length)       0.032783   0.026869   1.220 0.224197    
# scale(Conc.M)      -0.068919   0.020222  -3.408 0.000824 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.2526 on 163 degrees of freedom
# (12 observations deleted due to missingness)
# Multiple R-squared:  0.2091,	Adjusted R-squared:   0.18 
# F-statistic: 7.184 on 6 and 163 DF,  p-value: 8.19e-07
library(effects)
rodd.exp1.semd.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=rodd.exp1.semd.m, xlevels=list(nSenses_WMD=2)))
rodd.exp1.semd.meanings.e = as.data.frame(effect(term="Meanings", mod=rodd.exp1.semd.m))
rodd.exp1.semd.meanings.e$Meanings <- relevel(rodd.exp1.semd.meanings.e$Meanings, ref = "One")

cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_SemD_effects.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=4.5, width=9);
s = ggplot(rodd.exp1.semd.senses.e, aes(x=nSenses_WMD, y=fit)) + 
  ylim(0,3)+
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3, colour = NA) + 
  labs(x= "Number of Senses", y="Semantic Diversity") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") 
m= ggplot(rodd.exp1.semd.meanings.e, aes(x=Meanings, y=fit, fill=Meanings)) + 
  ylim(0,3)+
  geom_bar(stat="identity", position="dodge", color="white", alpha = 0.5) + 
  geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.2, position=position_dodge(width=0.9))+
  labs(x= "Meanings", y="Semantic Diversity") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none")
grid.arrange(s, m, nrow=1, ncol=2,
             top=textGrob("Semantic Diversity", gp=gpar(fontsize=20)))
dev.off()


rodd.exp1.semd$Senses = discretize(rodd.exp1.semd$nSenses_WMD, breaks = 2)
rodd.exp1.semd$Senses = mapvalues(rodd.exp1.semd$Senses, from = levels(rodd.exp1.semd$Senses), to = c("Few", "Many"))
rodd.semd.item <- aggregate(SemD ~ Word + Meanings + Senses + Freq + Length, FUN=mean, rodd.exp1.semd)
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_SemD.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=6, width=8);
par(mfrow=c(1,2))
pirateplot(formula = SemD ~ Senses,
           data = rodd.exp1.semd,
           ylab = "Semantic Diversity",
           xlab = "Senses",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
pirateplot(formula = SemD ~ Meanings,
           data = rodd.exp1.semd,
           ylab  = "Semantic Diversity",
           xlab = "Meanings",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")

dev.off();

#################################


# Replication & Simulation Analysis of Rodd et al. (2001) Experiment 2 (Visual Lexical Decision - Factorial Design) 
####################################################################################################################

# Get stimuli & properties
###############################
roddStimuli = subset(ambiguityStimuli, Authors=="Rodd et al." & Experiment==2)[c("Word", "AmbiguityType", "Meanings", "Senses")]
roddStimuli$AmbiguityType = as.factor(roddStimuli$AmbiguityType)
roddStimuli$Meanings = as.factor(roddStimuli$Meanings)
roddStimuli$Senses = as.factor(roddStimuli$Senses)
roddStimuli$Word = as.factor(roddStimuli$Word)
summary(roddStimuli)
###############################

#Plot Results from Rodd et al. (2001) #
#######################################
rt.m = c(587, 578, 586, 567)
rt.std = c(143, 135, 141, 129)
error.p = c(4.08, 1.77, 2.99, 1.63)
Senses = c("Few", "Many", "Few", "Many")
Meanings = c("Many", "Many", "One", "One")
rodd2001_results = data.frame(Senses, Meanings, rt.m, rt.std, error.p)   
rodd2001_results$ACC = abs(error.p-100)
rodd2001_results$Meanings <- factor(rodd2001_results$Meanings, levels = c("One", "Many"))
rodd2001_results$std.err = rodd2001_results$rt.std/sqrt(122)
rodd2001_results.senses <- aggregate(rt.m ~ Senses, FUN=mean, rodd2001_results)
rodd2001_results.senses$std.err <- aggregate(std.err ~ Senses, FUN=mean, rodd2001_results)$std.err
rodd2001_results.senses
rodd2001_results.meanings <- aggregate(rt.m ~ Meanings, FUN=mean, rodd2001_results)
rodd2001_results.meanings$std.err <- aggregate(std.err ~ Meanings, FUN=mean, rodd2001_results)$std.err
rodd2001_results.meanings

figFolder = paste(wd, "S1_SemanticDiversity_LSA/Exp2_Ambiguity&SemD/Figures/", sep="")
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Results_RT_Senses&Meanings.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=5, width=8);
s = ggplot(rodd2001_results.senses, aes(x=Senses, y=rt.m, fill=Senses)) +
  geom_bar(stat="identity", position=position_dodge(), color="white", alpha = 0.5) +
  geom_errorbar(aes(ymin=rt.m-std.err, ymax=rt.m+std.err), width=.15,
                position=position_dodge(.9)) + 
  coord_cartesian(ylim=c(500,650)) +
  ylab("Mean Reaction Time (ms)") +
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) +
  theme(legend.position = "none")
m = ggplot(rodd2001_results.meanings, aes(x=Meanings, y=rt.m, fill=Meanings)) +
  geom_bar(stat="identity", position=position_dodge(), color="white", alpha = 0.5) +
  geom_errorbar(aes(ymin=rt.m-std.err, ymax=rt.m+std.err), width=.15,
                position=position_dodge(.9)) + 
  coord_cartesian(ylim=c(500,650)) +
  ylab("Mean Reaction Time (ms)") +
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) +
  theme(legend.position = "none")
grid.arrange(s, m, nrow=1, ncol=2)
dev.off()

cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Results_RT_SensesXMeanings.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=6);
ggplot(rodd2001_results, aes(x=Senses, y=rt.m, fill=Meanings)) +
  geom_bar(stat="identity", position=position_dodge(),color="white", alpha = 0.5) +
  geom_errorbar(aes(ymin=rt.m-std.err, ymax=rt.m+std.err), width=.2,
                position=position_dodge(.9)) + 
  coord_cartesian(ylim=c(500,650)) +
  ylab("Mean Reaction Time (ms)")+ 
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) 
dev.off()

cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Results_Accuracy.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=5);
ggplot(rodd2001_results, aes(x=Senses, y=ACC, fill=Meanings)) +
  geom_bar(stat="identity", position=position_dodge(), color="white", alpha = 0.5) +
  coord_cartesian(ylim=c(80,100)) +
  ylab("Error Rate (%)") + 
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) 
dev.off()
#######################################

# BLP #
# Reaction Time Analysis #
##########################

# Get behavioural data from BLP
rodd.blp = merge(roddStimuli, blp.ld.data, by="Word", all.x = TRUE)
rodd.blp = merge(rodd.blp, subtlex[,c("Word","Freq")], by="Word", all.x = TRUE)
rodd.blp$Length = nchar(as.character(rodd.blp$Word))
summary(rodd.blp)
str(rodd.blp)
rodd.blp$SubjID = as.factor(rodd.blp$SubjID)
length(unique(rodd.blp$Word)) #128 w

#Check subject-level data
roddData.acc.subj = aggregate(ACC ~ SubjID, FUN=mean, rodd.blp)
summary(roddData.acc.subj) #no subjs with an error rate of more than 10%
roddData.rt.subj = aggregate(RT ~ SubjID, FUN=mean, rodd.blp)
summary(roddData.rt.subj) #no subjs with mean RT more than 1000ms

# Cleaning procedure 
#Remove incorrect responses
rodd.blp.rt = subset(rodd.blp, ACC==1)
#Remove responses longer than 1200ms and shorter than
rodd.blp.rt = subset(rodd.blp.rt, RT<1200)
summary(rodd.blp.rt)

#Check normality
hist(rodd.blp.rt$RT)
qqnorm(rodd.blp.rt$RT)
qqline(rodd.blp.rt$RT)
boxcox(rodd.blp.rt$RT ~ 1) #lambda ~ -1, inverse RT is recommended
hist(-1000/rodd.blp.rt$RT)
qqnorm(-1000/rodd.blp.rt$RT)
qqline(-1000/rodd.blp.rt$RT)

rodd.blp.rt$invRT = -1000/rodd.blp.rt$RT
summary(rodd.blp.rt)

#TRANSFORMED RT
#ByItem ANOVA
rodd.rt.item <- aggregate(invRT ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp.rt)
rodd.rt.item.aov = aov(invRT ~ Senses*Meanings + Freq + Length, data = rodd.rt.item)
summary(rodd.rt.item.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Senses            1 0.0831  0.0831  15.266 0.000154 ***
#   Meanings          1 0.0178  0.0178   3.278 0.072664 .  
# Freq              1 0.3886  0.3886  71.375 7.35e-14 ***
#   Length            1 0.0173  0.0173   3.177 0.077181 .  
# Senses:Meanings   1 0.0098  0.0098   1.806 0.181528    
# Residuals       122 0.6642  0.0054                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#BySubj ANOVA
rodd.rt.subj <- aggregate(invRT ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp.rt)
rodd.rt.subj.aov = aov(invRT ~  Senses*Meanings + Freq + Length, data = rodd.rt.subj)
summary(rodd.rt.subj.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Senses             1    3.1   3.071  18.864 1.43e-05 ***
#   Meanings           1    0.6   0.572   3.512   0.0610 .  
# Freq               1   13.7  13.655  83.875  < 2e-16 ***
#   Length             1    0.6   0.634   3.892   0.0486 *  
#   Senses:Meanings    1    0.3   0.345   2.121   0.1454    
# Residuals       4736  771.1   0.163                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

describeBy(rodd.blp.rt$RT, group = rodd.blp.rt$Senses) 
diff(as.matrix(aggregate(rodd.blp.rt$RT, by=list(rodd.blp.rt$Senses), FUN=mean)[2])) #15ms difference
describeBy(rodd.blp.rt$RT, group = rodd.blp.rt$Meanings) 
diff(as.matrix(aggregate(rodd.blp.rt$RT, by=list(rodd.blp.rt$Meanings), FUN=mean)[2])) #8ms difference

#Results
#Number of senses (F1(1,4726)=18.86, p<0.001, F2(1,121)=15.27, p<0.001) #15ms difference
#Number of meanings (F1(1,4726)=3.5, p=0.06, F2(1,121)=3.28, p=0.07)  #8ms difference


#UNTRANSFORMED RT
#ByItem ANOVA
rodd.rt.item <- aggregate(RT ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp.rt)
rodd.rt.item.aov = aov(RT ~ Senses*Meanings + Freq + Length, data = rodd.rt.item)
summary(rodd.rt.item.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Senses            1   8178    8178  14.114 0.000265 ***
#   Meanings          1   1248    1248   2.155 0.144725    
# Freq              1  37441   37441  64.617 6.63e-13 ***
#   Length            1   1149    1149   1.983 0.161643    
# Senses:Meanings   1   1774    1774   3.061 0.082705 .  
# Residuals       122  70691     579                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#BySubj ANOVA
rodd.rt.subj <- aggregate(RT ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp.rt)
rodd.rt.subj.aov = aov(RT ~  Senses*Meanings + Freq + Length, data = rodd.rt.subj)
summary(rodd.rt.subj.aov)
# Df   Sum Sq Mean Sq F value   Pr(>F)    
# Senses             1   305857  305857  17.173 3.47e-05 ***
#   Meanings           1    43321   43321   2.432   0.1189    
# Freq               1  1305787 1305787  73.318  < 2e-16 ***
#   Length             1    43121   43121   2.421   0.1198    
# Senses:Meanings    1    62006   62006   3.482   0.0621 .  
# Residuals       4736 84348211   17810                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#Plot
rodd.rt.item = aggregate(RT ~ Word + Senses + Meanings, FUN=mean, rodd.blp.rt)
rodd.rt.item$Meanings <- relevel(rodd.rt.item$Meanings, ref = "One")
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp2_BLP_RT.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=8);
#png(paste(figFolder,'Replication&Simulation_Rodd2001_BLP_RT.png',sep=''), 
#    width = 1200, height = 800, family='sans', type='cairo')
par(mfrow=c(1,2))
pirateplot(formula = RT ~ Senses,
           data = rodd.rt.item,
           ylab = "Mean Reaction Time (ms)",
           #xlab = "Condition",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
pirateplot(formula = RT ~ Meanings,
           data = rodd.rt.item,
           ylab = "Mean Reaction Time (ms)",
           #xlab = "Condition",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
dev.off();
##########################

# Accuracy Analysis #
#####################

#Get accuracy data
rodd.blp.acc = rodd.blp
summary(rodd.blp.acc)

#ByItem ANOVA
rodd.acc.item <- aggregate(ACC ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp)
rodd.acc.item.aov = aov(ACC ~ Senses*Meanings + Freq + Length, data = rodd.acc.item)
summary(rodd.acc.item.aov)
# Df  Sum Sq  Mean Sq F value   Pr(>F)    
# Senses            1 0.00727 0.007267   9.413 0.002655 ** 
#   Meanings          1 0.00002 0.000017   0.021 0.883772    
# Freq              1 0.01160 0.011598  15.022 0.000173 ***
#   Length            1 0.00110 0.001099   1.423 0.235152    
# Senses:Meanings   1 0.00027 0.000273   0.353 0.553356    
# Residuals       122 0.09419 0.000772                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#BySubj ANOVA
rodd.acc.subj <- aggregate(ACC ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp)
rodd.acc.subj.aov = aov(ACC ~ Senses*Meanings + Freq + Length, data = rodd.acc.subj)
summary(rodd.acc.subj.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Senses             1   0.30  0.3027  12.346 0.000446 ***
#   Meanings           1   0.00  0.0015   0.062 0.803376    
# Freq               1   0.43  0.4316  17.601 2.77e-05 ***
#   Length             1   0.05  0.0451   1.841 0.174873    
# Senses:Meanings    1   0.01  0.0105   0.427 0.513308    
# Residuals       4946 121.28  0.0245                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

describeBy(rodd.acc.item$ACC, group = rodd.acc.item$Senses) 
diff(as.matrix(aggregate(rodd.acc.item$ACC, by=list(rodd.acc.item$Senses), FUN=mean)[2])) #1.51% difference
describeBy(rodd.acc.item$ACC, group = rodd.acc.item$Meanings) 
diff(as.matrix(aggregate(rodd.acc.item$ACC, by=list(rodd.acc.item$Meanings), FUN=mean)[2])) #0.07% difference

#Results
#Number of senses (F1(1,4946)=12.35, p<0.001, F2(1,121)=9.41, p<0.01) #1.51% difference
#Number of meanings (F1(1,4946)=0.06, p=0.80, F2(1,121)=0.02, p=0.88)  
#####################

# ELP #
# Reaction Time Analysis #
##########################

# Get behavioural data from BLP
rodd.elp = merge(roddStimuli, elp.ld.data, by="Word", all.x = TRUE)
rodd.elp = merge(rodd.elp, subtlex[,c("Word","Freq")], by="Word", all.x = TRUE)
rodd.elp$Length = nchar(as.character(rodd.elp$Word))
summary(rodd.elp)
str(rodd.elp)
rodd.elp$SubjID = as.factor(rodd.elp$SubjID)
length(unique(rodd.elp$Word)) #128 w

#Check subject-level data
roddData.acc.subj = aggregate(ACC ~ SubjID, FUN=mean, rodd.elp)
summary(roddData.acc.subj) #no subjs with an error rate of more than 10%
roddData.rt.subj = aggregate(RT ~ SubjID, FUN=mean, rodd.elp)
summary(roddData.rt.subj)
#exclude subjects with mean RT more than 1000
subjToExclude = subset(roddData.rt.subj, RT>1000)$SubjID
rodd.elp = subset(rodd.elp, !(rodd.elp$SubjID %in% subjToExclude))

# Cleaning procedure 
#Remove incorrect responses
rodd.elp.rt = subset(rodd.elp, ACC==1)
#Remove responses longer than 1200ms and shorter than
rodd.elp.rt = subset(rodd.elp.rt, RT<1200)
summary(rodd.elp.rt)

#Check normality
hist(rodd.elp.rt$RT)
qqnorm(rodd.elp.rt$RT)
qqline(rodd.elp.rt$RT)
boxcox(rodd.elp.rt$RT ~ 1) #lambda ~ -1, inverse RT is recommended
hist(-1000/rodd.elp.rt$RT)
qqnorm(-1000/rodd.elp.rt$RT)
qqline(-1000/rodd.elp.rt$RT)

rodd.elp.rt$logRT = log(rodd.elp.rt$RT)
summary(rodd.elp.rt)

#TRANSFORMED RT
#ByItem ANOVA
rodd.rt.item <- aggregate(logRT ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp.rt)
rodd.rt.item.aov = aov(logRT ~ Senses*Meanings + Freq + Length, data = rodd.rt.item)
summary(rodd.rt.item.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Senses            1 0.0209 0.02087   6.334   0.0131 *  
#   Meanings          1 0.0073 0.00733   2.226   0.1383    
# Freq              1 0.0822 0.08219  24.952 1.98e-06 ***
#   Length            1 0.0092 0.00915   2.779   0.0981 .  
# Senses:Meanings   1 0.0058 0.00583   1.770   0.1859    
# Residuals       122 0.4019 0.00329                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#BySubj ANOVA
rodd.rt.subj <- aggregate(logRT ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp.rt)
rodd.rt.subj.aov = aov(logRT ~  Senses*Meanings + Freq + Length, data = rodd.rt.subj)
summary(rodd.rt.subj.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Senses             1   0.56  0.5563   7.230   0.0072 ** 
#   Meanings           1   0.23  0.2303   2.993   0.0837 .  
# Freq               1   2.30  2.3048  29.954 4.71e-08 ***
#   Length             1   0.31  0.3077   3.998   0.0456 *  
#   Senses:Meanings    1   0.17  0.1710   2.222   0.1361    
# Residuals       3764 289.62  0.0769                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

describeBy(rodd.elp.rt$RT, group = rodd.elp.rt$Senses) 
diff(as.matrix(aggregate(rodd.elp.rt$RT, by=list(rodd.elp.rt$Senses), FUN=mean)[2])) #16ms difference
describeBy(rodd.elp.rt$RT, group = rodd.elp.rt$Meanings) 
diff(as.matrix(aggregate(rodd.elp.rt$RT, by=list(rodd.elp.rt$Meanings), FUN=mean)[2])) #8ms difference

#Results
#Number of senses (F1(1,3763)=7.23, p<0.001, F2(1,121)=6.33, p<0.05) #14ms difference
#Number of meanings (F1(1,3764)=2.99, p=0.08, F2(1,121)=2.23, p=0.13)  #6ms difference


#UNTRANSFORMED RT
#ByItem ANOVA
rodd.rt.item <- aggregate(RT ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp.rt)
rodd.rt.item.aov = aov(RT ~ Senses*Meanings + Freq + Length, data = rodd.rt.item)
summary(rodd.rt.item.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Senses            1   9431    9431   8.474  0.00428 ** 
#   Meanings          1   1861    1861   1.673  0.19836    
# Freq              1  37025   37025  33.270 6.19e-08 ***
#   Length            1   6179    6179   5.553  0.02005 *  
#   Senses:Meanings   1    712     712   0.640  0.42526    
# Residuals       122 135768    1113                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#BySubj ANOVA
rodd.rt.subj <- aggregate(RT ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp.rt)
rodd.rt.subj.aov = aov(RT ~  Senses*Meanings + Freq + Length, data = rodd.rt.subj)
summary(rodd.rt.subj.aov)
# Df   Sum Sq Mean Sq F value   Pr(>F)    
# Senses             1   256316  256316  10.426  0.00125 ** 
#   Meanings           1    58743   58743   2.389  0.12225    
# Freq               1  1037529 1037529  42.201 9.31e-11 ***
#   Length             1   195199  195199   7.940  0.00486 ** 
#   Senses:Meanings    1    19647   19647   0.799  0.37141    
# Residuals       3764 92539575   24585                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#Plot
rodd.rt.item = aggregate(RT ~ Word + Senses + Meanings, FUN=mean, rodd.elp.rt)
rodd.rt.item$Meanings <- relevel(rodd.rt.item$Meanings, ref = "One")
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp2_ELP_RT.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=8);
#png(paste(figFolder,'Replication&Simulation_Rodd2001_BLP_RT.png',sep=''), 
#    width = 1200, height = 800, family='sans', type='cairo')
par(mfrow=c(1,2))
pirateplot(formula = RT ~ Senses,
           data = rodd.rt.item,
           ylab = "Mean Reaction Time (ms)",
           #xlab = "Condition",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
pirateplot(formula = RT ~ Meanings,
           data = rodd.rt.item,
           ylab = "Mean Reaction Time (ms)",
           #xlab = "Condition",
           #ylim = c(400,825),
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
dev.off();
##########################

# Accuracy Analysis #
#####################

#Get accuracy data
rodd.elp.acc = rodd.elp
summary(rodd.elp.acc)

#ByItem ANOVA
rodd.acc.item <- aggregate(ACC ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp)
rodd.acc.item.aov = aov(ACC ~ Senses*Meanings + Freq + Length, data = rodd.acc.item)
summary(rodd.acc.item.aov)
# Df  Sum Sq Mean Sq F value   Pr(>F)    
# Senses            1 0.00706 0.00706   4.426   0.0374 *  
#   Meanings          1 0.00024 0.00024   0.152   0.6973    
# Freq              1 0.05472 0.05472  34.319 4.06e-08 ***
#   Length            1 0.00180 0.00180   1.132   0.2895    
# Senses:Meanings   1 0.00008 0.00008   0.051   0.8218    
# Residuals       122 0.19451 0.00159                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#BySubj ANOVA
rodd.acc.subj <- aggregate(ACC ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp)
rodd.acc.subj.aov = aov(ACC ~ Senses*Meanings + Freq + Length, data = rodd.acc.subj)
summary(rodd.acc.subj.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Senses             1   0.23  0.2319   6.488   0.0109 *  
#   Meanings           1   0.01  0.0102   0.285   0.5936    
# Freq               1   1.81  1.8060  50.537 1.37e-12 ***
#   Length             1   0.06  0.0616   1.724   0.1892    
# Senses:Meanings    1   0.00  0.0027   0.077   0.7816    
# Residuals       4169 148.98  0.0357                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

describeBy(rodd.acc.item$ACC, group = rodd.acc.item$Senses) 
diff(as.matrix(aggregate(rodd.acc.item$ACC, by=list(rodd.acc.item$Senses), FUN=mean)[2])) #1.4% difference
describeBy(rodd.acc.item$ACC, group = rodd.acc.item$Meanings) 
diff(as.matrix(aggregate(rodd.acc.item$ACC, by=list(rodd.acc.item$Meanings), FUN=mean)[2])) #0.2% difference

#Results
#Number of senses (F1(1,4168)=6.4, p<0.05, F2(1,121)=4.43, p<0.05) #1.4% difference
#Number of meanings (F1(1,4168)=0.28, p=0.59, F2(1,121)=0.15, p=0.70)  
#####################


# Semantic Diversity Simulation #
#################################

# Get data 
rodd.semd = merge(roddStimuli, lexicalVars, by="Word", all.x = TRUE)
summary(rodd.semd)

#ByItem ANOVA
rodd.semd.item <- aggregate(SemD ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.semd)
rodd.semd.item.aov = aov(SemD ~ Senses*Meanings + Freq + Length, data = rodd.semd.item)
summary(rodd.semd.item.aov)
# Df Sum Sq Mean Sq F value Pr(>F)  
# Senses            1  0.000 0.00000   0.000 0.9953  
# Meanings          1  0.022 0.02176   0.385 0.5359  
# Freq              1  0.104 0.10364   1.835 0.1780  
# Length            1  0.159 0.15883   2.813 0.0961 .
# Senses:Meanings   1  0.001 0.00102   0.018 0.8935  
# Residuals       122  6.889 0.05647                 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Results
#SemD (F(1,121)<0.001, p=0.99)

describeBy(rodd.semd.item$SemD, group = rodd.semd.item$Senses)
diff(as.matrix(aggregate(rodd.semd.item$SemD, by=list(rodd.semd.item$Senses), FUN=mean)[2])) #0.0002 difference
describeBy(rodd.semd.item$SemD, group = rodd.semd.item$Meanings) 
diff(as.matrix(aggregate(rodd.semd.item$SemD, by=list(rodd.semd.item$Meanings), FUN=mean)[2])) #0.0260 difference

#Plot
rodd.semd.item$Meanings <- relevel(rodd.semd.item$Meanings, ref = "One")
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp2_SemD.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=8);
par(mfrow=c(1,2))
pirateplot(formula = SemD ~ Senses,
           data = rodd.semd.item,
           ylab = "Semantic Diversity",
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
pirateplot(formula = SemD ~ Meanings,
           data = rodd.semd.item,
           ylab = "Semantic Diversity",
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
dev.off()
#################################




# Replication Analysis of Armstrong & Plaut (2016)
##############################################################################################################################################

#Clean workspace
rm(list=setdiff(ls(), c("wd","figFolder","lexicalVars","ambiguityStimuli","blp.ld.data", "elp.ld.data")))

# Get stimuli & properties
###############################
armstrongStimuli = subset(ambiguityStimuli, Authors=="Armstrong & Plaut")[c("Word", "AmbiguityType")]
armstrongStimuli$AmbiguityType = as.factor(armstrongStimuli$AmbiguityType)
armstrongStimuli$Word = as.factor(armstrongStimuli$Word)
summary(armstrongStimuli)
#length
armstrongStimuli$Length = nchar(as.character(armstrongStimuli$Word))
#familiarity norms
famNorms = read.csv("Resources/Stimuli/ArmstrongPlaut2016/fam.itemMeans.txt", sep="\t", header=TRUE)
setnames(famNorms, c("word"), c("Word"))
library(stringr)
famNorms$Senses = str_split_fixed(famNorms$category,'g', 2)[,2]
famNorms$Senses = as.factor(famNorms$Senses)
levels(famNorms$Senses) <- c("Few","Many")
famNorms$Meanings = str_split_fixed(famNorms$category,'g', 2)[,1]
famNorms$Meanings = as.factor(famNorms$Meanings)
levels(famNorms$Meanings) <- c("Many","One")
armstrongStimuli = merge(armstrongStimuli, famNorms[,c("Word","Meanings","Senses","fSense","famfRes")], by="Word", all.x = TRUE) 
summary(armstrongStimuli)
meaningFreq = read_excel("Resources/Norms/Measures/Rice2019_MeaningFrequencyEstimates.xlsx", sheet = "MeaningFrequencyRatings")
setnames(meaningFreq,
c("word","eDom biggest","Log10 SUBTL word freq","Coltheart's N","Number of letters","Number of phonemes",
  "Number of syllables","Number of meanings","Number of senses"),
c("Word","eDom_biggest","LogFreq_Subtlex","colN","nLet","nPhon","nSyll","nMeanings","nSenses"))
armstrongStimuli = merge(armstrongStimuli,
                         meaningFreq[,c("Word","eDom_biggest","colN","nLet","nSyll","nMeanings",
                                        "nSenses","OLD","Imageability")],
                         by="Word", all.x = TRUE)
armstrongStimuli$nLet = as.numeric(as.character(armstrongStimuli$nLet))
summary(armstrongStimuli)
str(armstrongStimuli)
#word frequency
subtlex = read.csv("Resources/Norms/Databases/SUBTLEX/SUBTLEX-UK.csv")
setnames(subtlex, c("Spelling","LogFreq.Zipf."), c("Word","ZipfFreq_Subtlex"))
armstrongStimuli = merge(armstrongStimuli, subtlex[,c("Word","ZipfFreq_Subtlex")], by="Word", all.x = TRUE)
summary(armstrongStimuli)
#OLD20
blpStimuli = read.table("Resources/Data/BritishLexiconProject/blp-stimuli.txt/blp-stimuli.txt", sep="\t", header=TRUE)
setnames(blpStimuli, c("spelling", "coltheart.N"), c("Word", "coltN"))
blpStimuli = blpStimuli[c("Word","OLD20","nsyl", "coltN")]
armstrongStimuli = merge(armstrongStimuli, blpStimuli, by="Word", all.x = TRUE)
summary(armstrongStimuli)
#number of meanings/senses
wordsmythDict = read_excel("Resources/Norms/Dictionaries/WordsmythDictionary/Wordsmyth_Words_NumMeaning_NumSenses_PartsOfSpeechFreq.xls")
setnames(wordsmythDict, c("word", "NumMeanings", "NumSenses"), c("Word", "nMeanings_WMD", "nSenses_WMD"))
wordsmythDict = wordsmythDict[,c("Word","nMeanings_WMD","nSenses_WMD")]
armstrongStimuli = merge(armstrongStimuli, wordsmythDict, by="Word", all.x = TRUE)
summary(armstrongStimuli)
#number of phonemes
celex_epw<-read.table("Resources/Norms/Databases/CELEX/EPW/EPW.CD",sep="\\", fill=T, stringsAsFactors=F, quote="",comment.char = "")
setnames(celex_epw, c("V2", "V10"), c("Word", "nPhon"))
celex_epw = celex_epw[c("Word", "nPhon")]
armstrongStimuli = merge(armstrongStimuli, celex_epw, by="Word", all.x = TRUE)
summary(armstrongStimuli)
###############################


#Plot Results from Armstrong & Plaut (2016)
###########################################
armstrongResults = read.table("Resources/Stimuli/ArmstrongPlaut2016/wordStimuli_cleaned.ACC.txt", header=TRUE) #.grandaverage
rts = read.table("Resources/Stimuli/ArmstrongPlaut2016/wordStimuli_cleaned.correctRT.txt", header=TRUE) #.grandaverage
armstrongResults = merge(armstrongResults, rts, by=c("string", "stimClass", "cont", "nwDiff"))
setnames(armstrongResults, c("string","stimClass","acc","cont"), c("Word","AmbiguityType","ACC","Contrast"))
armstrongResults$AmbiguityType = revalue(armstrongResults$AmbiguityType, c("h"="Homonymy", "p"="Polysemy", "u"="Unambiguous", "y"="Hybrid"))
armstrongResults$AmbiguityType <- factor(armstrongResults$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy", "Hybrid"))
armstrongResults$nwDiff = revalue(armstrongResults$nwDiff, c("h"="Hard", "p"="Pseudohomophone", "e"="Very Hard"))
armstrongResults$nwDiff <- factor(armstrongResults$nwDiff, levels = c("Hard", "Very Hard", "Pseudohomophone"))
armstrongResults$Contrast = revalue(armstrongResults$Contrast, c("deg"="Degreded", "full"="Full"))
armstrongResults$Contrast <- factor(armstrongResults$Contrast, levels = c("Degreded", "Full"))
#meaning frequency estimates
mFreq = read.table("Resources/Stimuli/ArmstrongPlaut2016/wordStimuli_cleaned.grandaverage.ACC.txt", header=TRUE)
setnames(mFreq, c("string"), c("Word"))
armstrongResults = merge(armstrongResults, mFreq[,c("Word","biggest")], by="Word", all.x = TRUE) 
describeBy(armstrongResults$RT, group = armstrongResults$Condition)
#subset word with dominance meaning frequency less than 65 for visualisation purposes
plotSubset = subset(armstrongResults, biggest<65 & AmbiguityType!="Hybrid" | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous" )
#Plot across conditions
mean.rt.ambiguity <- aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(mean.rt.ambiguity, c("RT"), c("rt.mean"))
sd.rt.ambiguity = aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(sd.rt.ambiguity, c("RT"), c("rt.sd"))
armstrongResults.rt.ambiguity <- merge(mean.rt.ambiguity, sd.rt.ambiguity)
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_Results_RT.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=4.5, width=5);
armstrongResults.rt = ggplot(armstrongResults.rt.ambiguity,
       aes(x=AmbiguityType, y=rt.mean, fill=AmbiguityType, color=AmbiguityType)) +
  geom_bar(stat="identity", position=position_dodge(), color= "white", alpha = 0.5) +
  geom_errorbar(aes(ymin=rt.mean-rt.sd, ymax=rt.mean+rt.sd), width=.2,
                position=position_dodge(.9), color="black") + 
  ggtitle("Armstrong & Plaut (2016)") +
  theme(
    plot.title = element_text(size=20, hjust = 0.5),
    axis.text=element_text(size=15), axis.title=element_text(size=15), legend.position = "none") + 
  coord_cartesian(ylim=c(500,650)) +
  ylab("Mean Reaction Time (ms)") +
  xlab("Condition") +
  scale_fill_manual(values=c("#0C5BB0FF", "#EE0011FF","#15983DFF"))
grid.arrange(armstrongResults.rt, nrow=1, ncol=1)
dev.off()



###########################################

# BLP #
# Data screening #
##################

# Get behavioural data from BLP
armstrong.blp = merge(armstrongStimuli, blp.ld.data, by="Word", all.x = TRUE)
#exclude hybrids as not of interest
armstrong.blp = subset(armstrong.blp, AmbiguityType!="Hybrid")
armstrong.blp$AmbiguityType <- factor(armstrong.blp$AmbiguityType)
armstrong.blp$SubjID = as.factor(armstrong.blp$SubjID)
summary(armstrong.blp)
str(armstrong.blp)

# Cleaning procedure 
#Exlude outliers in speed-accuracy space using a Mahalanobis distance statistic
#Item-level
armstrongData.acc.item = aggregate(ACC ~ Word + AmbiguityType, FUN=mean, armstrong.blp)
armstrongData.rt.item = aggregate(RT ~ Word + AmbiguityType, FUN=mean, armstrong.blp)
summary(armstrongData.acc.item) 
summary(armstrongData.rt.item) 
armstrongData.item = merge(armstrongData.acc.item, armstrongData.rt.item)
itemsToexlude = setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Word","AmbiguityType","ACC","RT","mDist","cValue"))
nItemExluded = 0
for (aType in unique(armstrongData.item$AmbiguityType)){
  print(aType)
  s = subset(armstrongData.item, AmbiguityType==aType)
  print(head(s))
  s$mDist = mahalanobis(s[, 3:4], colMeans(s[, 3:4]), cov(s[, 3:4]))
  s$cValue = qchisq(1-.01,df=ncol(s[, 3:4]))
  print(nrow(s))
  print(nrow(subset(s, mDist>cValue)))
  nItemExluded = nItemExluded + nrow(subset(s, mDist>cValue))
  itemsToexlude = rbind(itemsToexlude, s)
}
nItemExluded #8 items
itemsToexlude = subset(itemsToexlude, mDist>cValue)
armstrong.blp.clean = armstrong.blp[!(armstrong.blp$Word %in% itemsToexlude$Word),]
#Subject-level
armstrongData.acc.subj = aggregate(ACC ~ SubjID + AmbiguityType, FUN=mean, armstrong.blp)
armstrongData.rt.subj = aggregate(RT ~ SubjID + AmbiguityType, FUN=mean, armstrong.blp)
summary(armstrongData.acc.subj) 
summary(armstrongData.rt.subj) 
armstrongData.subj= merge(armstrongData.acc.subj, armstrongData.rt.subj)
subjsToexlude = setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("SubjID","AmbiguityType","ACC","RT","mDist","cValue"))
nSubjExluded = 0
for (aType in unique(armstrongData.subj$AmbiguityType)){
  print(aType)
  s = subset(armstrongData.subj, AmbiguityType==aType)
  print(head(s))
  s$mDist = mahalanobis(s[, 3:4], colMeans(s[, 3:4]), cov(s[, 3:4]))
  s$cValue = qchisq(1-.01,df=ncol(s[, 3:4]))
  print(nrow(s))
  print(nrow(subset(s, mDist>cValue)))
  nSubjExluded = nSubjExluded + nrow(subset(s, mDist>cValue))
  subjsToexlude = rbind(subjsToexlude, s)
}
nSubjExluded #7 subjs
subjsToexlude = subset(subjsToexlude, mDist>cValue)
armstrong.blp.clean = armstrong.blp.clean[!(armstrong.blp.clean$SubjID %in% subjsToexlude$SubjID),]
((nrow(armstrong.blp)-nrow(armstrong.blp.clean))/nrow(armstrong.blp))*100 # 8.9% data removed
##################

# Reaction Time Analysis #
##########################

armstrong.blp.rt = armstrong.blp.clean[!(is.na(armstrong.blp.clean$RT)),]
nrow(armstrong.blp.rt)
#Remove incorrect responses as well as responses longer than 1200ms and shorter than 200ms
armstrong.blp.rt.clean = subset(armstrong.blp.rt, RT<2000 & RT>200)
nrow(armstrong.blp.rt.clean)
#Remove trial outliers within each block that exceeded the z-score associated with a two tailed p-value of .005
library(outliers)
trialClean = setNames(data.frame(matrix(ncol = ncol(armstrong.blp.rt.clean), nrow = 0)), colnames(armstrong.blp.rt.clean))
trialExcluded = 0
for (b in unique(armstrong.blp.rt.clean$Block)){
    s = subset(armstrong.blp.rt.clean, Block==b)
    trialExcluded = trialExcluded + nrow(s[!scores(s$RT, type="chisq", prob=1-.005),])
    s = s[!scores(s$RT, type="chisq", prob=1-.005),]
    trialClean = rbind(trialClean, s)
  }
nrow(trialClean)
armstrong.blp.rt.clean = trialClean
((nrow(armstrong.blp.rt)-nrow(armstrong.blp.rt.clean))/nrow(armstrong.blp.rt))*100 # 2.7% data removed
summary(armstrong.blp.rt.clean)

#Descriptives
describeBy(armstrong.blp.rt.clean$RT, group = armstrong.blp.rt.clean$AmbiguityType)
#Homonymy (M = 541ms, SD = 123ms)
#Polysemy (M = 532ms, SD = 111ms)
#Unambiguous (M = 540ms, SD = 117ms)

#Check normality
hist(armstrong.blp.rt.clean$RT)
qqnorm(armstrong.blp.rt.clean$RT)
qqline(armstrong.blp.rt.clean$RT)
boxcox(armstrong.blp.rt.clean$RT ~ 1) #lambda ~ -1, inverse RT is recommended
hist(-1000/armstrong.blp.rt.clean$RT)
qqnorm(-1000/armstrong.blp.rt.clean$RT)
qqline(-1000/armstrong.blp.rt.clean$RT)
armstrong.blp.rt.clean$invRT = -1000/armstrong.blp.rt.clean$RT

# #Set polysemes and unambiguous word to have a single completely dominant interpretation (biggest=100)
armstrong.blp.rt.clean[["eDom_biggest"]][is.na(armstrong.blp.rt.clean[["eDom_biggest"]])] = 100
describeBy(armstrong.blp.rt.clean$eDom_biggest, group = armstrong.blp.rt.clean$AmbiguityType)

#Modelling followeing Armstorng & Plaut (2016) procedure
armstrong.blp.rt.clean$AmbiguityType <- relevel(armstrong.blp.rt.clean$AmbiguityType, ref = "Unambiguous")
#get the inverse of dominant frequency to have higher values for balanced homonymous
armstrong.blp.rt.clean$invDom = 1000/(armstrong.blp.rt.clean$eDom_biggest)

# UNAMBIGUOUS - HOMONYMS
hSubset = subset(armstrong.blp.rt.clean, AmbiguityType=="Unambiguous" | AmbiguityType=="Homonymy")
armstrong.rt.m.h = lmer(invRT ~ 
                          #Homonymy
                          scale(invDom) +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) +  #lexicality of previous trial removed
                          #random effects
                          (1|SubjID) + (1|Word), 
                        hSubset)
summary(armstrong.rt.m.h)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)             -1.925e+00  2.165e-02  7.953e+01 -88.930  < 2e-16 ***
#   scale(invDom)            7.113e-03  5.331e-03  1.644e+02   1.334 0.183993    
# scale(ZipfFreq_Subtlex) -4.462e-02  6.131e-03  1.613e+02  -7.277 1.41e-11 ***
#   scale(OLD20)            -9.732e-03  7.002e-03  1.657e+02  -1.390 0.166421    
# scale(nsyl)              5.513e-03  1.023e-02  2.185e+03   0.539 0.590038    
# scale(nPhon)            -2.437e-03  9.017e-03  2.260e+04  -0.270 0.787002    
# scale(Length)            8.615e-03  6.880e-03  1.866e+02   1.252 0.212051    
# scale(famfRes)          -2.620e-02  5.987e-03  1.679e+02  -4.376 2.12e-05 ***
#   scale(pRT)               6.719e-02  2.025e-03  2.614e+04  33.183  < 2e-16 ***
#   scale(pACC)              7.589e-03  1.952e-03  2.611e+04   3.889 0.000101 ***
#   scale(Trial)            -1.248e-02  9.784e-03  2.192e+02  -1.275 0.203575    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.blp.rt.clean, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
armstrong.rt.m.p = lmer(invRT ~ 
                          #Polysemy
                          AmbiguityType +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) +  #lexicality of previous trial removed
                          #random effects
                          (1|SubjID) + (1|Word), 
                        pSubset)

summary(armstrong.rt.m.p)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)             -1.925e+00  2.049e-02  8.954e+01 -93.973  < 2e-16 ***
#   AmbiguityTypePolysemy   -2.503e-02  9.396e-03  1.721e+02  -2.664 0.008461 ** 
#   scale(ZipfFreq_Subtlex) -3.207e-02  5.428e-03  1.700e+02  -5.909 1.83e-08 ***
#   scale(OLD20)            -4.705e-03  6.547e-03  1.711e+02  -0.719 0.473360    
# scale(nsyl)              9.825e-03  1.537e-02  9.613e+03   0.639 0.522693    
# scale(nPhon)            -4.154e-03  1.470e-02  2.830e+04  -0.283 0.777530    
# scale(Length)            3.912e-03  6.622e-03  1.763e+02   0.591 0.555394    
# scale(famfRes)          -2.005e-02  5.527e-03  1.743e+02  -3.629 0.000374 ***
#   scale(pRT)               6.949e-02  1.869e-03  2.962e+04  37.186  < 2e-16 ***
#   scale(pACC)              6.131e-03  1.802e-03  2.958e+04   3.403 0.000667 ***
#   scale(Trial)            -2.532e-03  8.153e-03  2.093e+02  -0.311 0.756469    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Results
#Unambiguous - Homonymy (b<0.01,SE<0.01,t=1.33,p=0.18)
#Unambiguous - Polysemy (b=0.02,SE<0.01,t=-5.91,p<0.01)

#Print table with models summary
tab_model(armstrong.rt.m.p, armstrong.rt.m.h, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "Armstrong & Plaut (2016) Replication Analysis on BLP Lexical Decision Data")

#Plot
figFolder = paste(wd, "S1_SemanticDiversity_LSA/Exp2_Ambiguity&SemD/Figures/", sep="")
armstrong.blp.rt.clean$AmbiguityType <- factor(armstrong.blp.rt.clean$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
plotSubset = subset(armstrong.blp.rt.clean, eDom_biggest<65 | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous")
armstrong.blp.rt.mean <- aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.blp.rt.mean, c("RT"), c("rt.mean"))
armstrong.blp.rt.se = aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.blp.rt.se, c("RT"), c("rt.se"))
armstrong.blp.rt.ambiguity <- merge(armstrong.blp.rt.mean, armstrong.blp.rt.se)
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_BLP_RT.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=4.5, width=5);
armstrong.blp.rt.p = ggplot(armstrong.blp.rt.ambiguity,
       aes(x=AmbiguityType, y=rt.mean, fill=AmbiguityType, color=AmbiguityType)) +
  geom_bar(stat="identity", position=position_dodge(), color= "white", alpha = 0.5) +
  geom_errorbar(aes(ymin=rt.mean-rt.se, ymax=rt.mean+rt.se), width=.2,
                position=position_dodge(.9), color="black") + 
  ggtitle("BLP Lexical Decision") +
  theme(
    plot.title = element_text(size=20, hjust = 0.5),
    axis.text=element_text(size=15), axis.title=element_text(size=15), legend.position = "none") +  
  coord_cartesian(ylim=c(500,650)) +
  ylab("Mean Reaction Time (ms)") +
  xlab("Condition") +
  scale_fill_manual(values=c("#0C5BB0FF", "#EE0011FF","#15983DFF"))+
  geom_signif(y_position=c(560), xmin=c(1.8), xmax=c(3.2),
              annotation=c("**"), tip_length=0, size = 0.7, textsize = 12, color = "black")
grid.arrange(armstrong.blp.rt.p, nrow=1, ncol=1)
dev.off()

##########################

# Accuracy Analysis #
##########################

armstrong.blp.acc = armstrong.blp.clean
summary(armstrong.blp.acc)
#Descriptives
describeBy(armstrong.blp.acc$ACC, group = armstrong.blp.acc$AmbiguityType)
#Homonymy (M = 97%, SD = 0.17)
#Polysemy (M = 99%, SD = 0.12)
#Unambiguous (M = 98%, SD = 0.15)

#Set polysemes and unambiguous word to have a single completely dominant interpretation (biggest=100)
armstrong.blp.acc[["eDom_biggest"]][is.na(armstrong.blp.acc[["eDom_biggest"]])] = 100
describeBy(armstrong.blp.acc$eDom_biggest, group = armstrong.blp.acc$AmbiguityType)

#Modelling followeing Armstorng & Plaut (2016) procedure
armstrong.blp.acc$AmbiguityType <- relevel(armstrong.blp.acc$AmbiguityType, ref = "Unambiguous")
#get the inverse of dominant frequency to have higher values for balanced homonymous
armstrong.blp.acc$invDom = 1000/(armstrong.blp.acc$eDom_biggest)

# UNAMBIGUOUS - HOMONYMS
hSubset = subset(armstrong.blp.acc, AmbiguityType=="Unambiguous" | AmbiguityType=="Homonymy")
armstrong.acc.m.h = glmer(ACC ~ 
                          #Homonymy
                          scale(invDom) +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) +  #lexicality of previous trial removed
                          #random effects
                          (1|SubjID) + (1|Word), 
                          hSubset,
                          family="binomial", control=glmerControl(optimizer="bobyqa",calc.derivs=F))
summary(armstrong.acc.m.h)

# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              4.96212    0.19118  25.955  < 2e-16 ***
#   scale(invDom)            0.01184    0.13444   0.088 0.929850    
# scale(ZipfFreq_Subtlex)  0.47750    0.15612   3.059 0.002224 ** 
#   scale(OLD20)             0.09152    0.18247   0.502 0.615983    
# scale(nsyl)              0.00952    0.26243   0.036 0.971061    
# scale(nPhon)             0.04734    0.23292   0.203 0.838946    
# scale(Length)            0.07528    0.18412   0.409 0.682642    
# scale(famfRes)           0.27508    0.14604   1.884 0.059624 .  
# scale(pRT)               0.19561    0.05351   3.655 0.000257 ***
#   scale(pACC)              0.17099    0.03499   4.886 1.03e-06 ***
#   scale(Trial)            -0.07596    0.17313  -0.439 0.660868    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.blp.acc, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
armstrong.acc.m.p = glmer(ACC ~ 
                            #Polysemy
                            AmbiguityType +
                            #control variable
                            scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes) + 
                            #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                            #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                            scale(pRT) + scale(pACC) + scale(Trial) +  
                            #random effects
                            (1|SubjID) + (1|Word), 
                          pSubset,
                          family="binomial", control=glmerControl(optimizer="bobyqa",calc.derivs=F))
summary(armstrong.acc.m.p)
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               5.69664    1.61413   3.529 0.000417 ***
#   AmbiguityTypePolysemy     1.03718    0.37263   2.783 0.005379 ** 
#   scale(ZipfFreq_Subtlex)   0.70835    0.21237   3.335 0.000852 ***
#   scale(OLD20)              0.03958    0.26272   0.151 0.880259    
# scale(nsyl)              -3.22365  338.80531  -0.010 0.992408    
# scale(nPhon)              3.33023  341.21356   0.010 0.992213    
# scale(Length)             0.27796    0.27789   1.000 0.317190    
# scale(famfRes)            0.08063    0.20679   0.390 0.696606    
# scale(pRT)                0.42960    0.07722   5.563 2.64e-08 ***
#   scale(pACC)               0.15031    0.04236   3.549 0.000387 ***
#   scale(Trial)             -0.26423    0.22769  -1.161 0.245840    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Results
#Unambiguous - Homonymy (b=0.01,SE=0.13,t=0.08,p=92)
#Unambiguous - Polysemy (b=1.03,SE=0.37,t=2.78,p<0.01)

#Print table with models summary
tab_model(armstrong.acc.m.p, armstrong.acc.m.h, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "Armstrong & Plaut (2016) Replication Analysis on BLP Lexical Decision Data")
##########################

# ELP #
# Data screening #
##################

# Get behavioural data from ELP
armstrong.elp = merge(armstrongStimuli, elp.ld.data, by="Word", all.x = TRUE)
#exclude hybrids as not of interest
armstrong.elp = subset(armstrong.elp, AmbiguityType!="Hybrid")
armstrong.elp$AmbiguityType <- factor(armstrong.elp$AmbiguityType)
armstrong.elp$SubjID = as.factor(armstrong.elp$SubjID)
summary(armstrong.elp)
str(armstrong.elp)

# Cleaning procedure 
#Exlude outliers in speed-accuracy space using a Mahalanobis distance statistic
#Item-level
armstrongData.acc.item = aggregate(ACC ~ Word + AmbiguityType, FUN=mean, armstrong.elp)
armstrongData.rt.item = aggregate(RT ~ Word + AmbiguityType, FUN=mean, armstrong.elp)
summary(armstrongData.acc.item) 
summary(armstrongData.rt.item) 
armstrongData.item = merge(armstrongData.acc.item, armstrongData.rt.item)
itemsToexlude = setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Word","AmbiguityType","ACC","RT","mDist","cValue"))
nItemExluded = 0
for (aType in unique(armstrongData.item$AmbiguityType)){
  print(aType)
  s = subset(armstrongData.item, AmbiguityType==aType)
  print(head(s))
  s$mDist = mahalanobis(s[, 3:4], colMeans(s[, 3:4]), cov(s[, 3:4]))
  s$cValue = qchisq(1-.01,df=ncol(s[, 3:4]))
  print(nrow(s))
  print(nrow(subset(s, mDist>cValue)))
  nItemExluded = nItemExluded + nrow(subset(s, mDist>cValue))
  itemsToexlude = rbind(itemsToexlude, s)
}
nItemExluded #8 items
itemsToexlude = subset(itemsToexlude, mDist>cValue)
armstrong.elp.clean = armstrong.elp[!(armstrong.elp$Word %in% itemsToexlude$Word),]
#Subject-level
armstrongData.acc.subj = aggregate(ACC ~ SubjID + AmbiguityType, FUN=mean, armstrong.elp)
armstrongData.rt.subj = aggregate(RT ~ SubjID + AmbiguityType, FUN=mean, armstrong.elp)
summary(armstrongData.acc.subj) 
summary(armstrongData.rt.subj) 
armstrongData.subj= merge(armstrongData.acc.subj, armstrongData.rt.subj)
subjsToexlude = setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("SubjID","AmbiguityType","ACC","RT","mDist","cValue"))
nSubjExluded = 0
for (aType in unique(armstrongData.subj$AmbiguityType)){
  print(aType)
  s = subset(armstrongData.subj, AmbiguityType==aType)
  print(head(s))
  s$mDist = mahalanobis(s[, 3:4], colMeans(s[, 3:4]), cov(s[, 3:4]))
  s$cValue = qchisq(1-.01,df=ncol(s[, 3:4]))
  print(nrow(s))
  print(nrow(subset(s, mDist>cValue)))
  nSubjExluded = nSubjExluded + nrow(subset(s, mDist>cValue))
  subjsToexlude = rbind(subjsToexlude, s)
}
nSubjExluded #7 subjs
subjsToexlude = subset(subjsToexlude, mDist>cValue)
armstrong.elp.clean = armstrong.elp.clean[!(armstrong.elp.clean$SubjID %in% subjsToexlude$SubjID),]
((nrow(armstrong.elp)-nrow(armstrong.elp.clean))/nrow(armstrong.elp))*100 # 11.72% data removed
##################

# Reaction Time Analysis #
##########################

armstrong.elp.rt = armstrong.blp.clean[!(is.na(armstrong.elp.clean$RT)),]
nrow(armstrong.elp.rt)
#Remove incorrect responses as well as responses longer than 1200ms and shorter than 200ms
armstrong.elp.rt.clean = subset(armstrong.elp.rt, RT<2000 & RT>200)
nrow(armstrong.elp.rt.clean)
((nrow(armstrong.elp.rt)-nrow(armstrong.elp.rt.clean))/nrow(armstrong.elp.rt))*100 # 0.3% data removed
summary(armstrong.elp.rt.clean)

#Check normality
hist(armstrong.elp.rt.clean$RT)
qqnorm(armstrong.elp.rt.clean$RT)
qqline(armstrong.elp.rt.clean$RT)
boxcox(armstrong.elp.rt.clean$RT ~ 1) #lambda ~ -1, inverse RT is recommended
hist(-1000/armstrong.elp.rt.clean$RT)
qqnorm(-1000/armstrong.elp.rt.clean$RT)
qqline(-1000/armstrong.elp.rt.clean$RT)
armstrong.elp.rt.clean$invRT = -1000/armstrong.elp.rt.clean$RT

#Set polysemes and unambiguous word to have a single completely dominant interpretation (biggest=100)
armstrong.elp.rt.clean[["eDom_biggest"]][is.na(armstrong.elp.rt.clean[["eDom_biggest"]])] = 100
describeBy(armstrong.elp.rt.clean$eDom_biggest, group = armstrong.elp.rt.clean$AmbiguityType)

#Modelling followeing Armstorng & Plaut (2016) procedure
armstrong.elp.rt.clean$AmbiguityType <- relevel(armstrong.elp.rt.clean$AmbiguityType, ref = "Unambiguous")
#get the inverse of dominant frequency to have higher values for balanced homonymous
armstrong.elp.rt.clean$invDom = 1000/(armstrong.elp.rt.clean$eDom_biggest)
summary(armstrong.elp.rt.clean)

# UNAMBIGUOUS - HOMONYMS
hSubset = subset(armstrong.elp.rt.clean, AmbiguityType=="Unambiguous" | AmbiguityType=="Homonymy")
armstrong.rt.m.h = lmer(invRT ~ 
                          #Homonymy
                          scale(invDom) +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) +  
                          #random effects
                          (1|SubjID) + (1|Word), 
                        hSubset)
summary(armstrong.rt.m.h)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)             -1.900e+00  2.321e-02  7.993e+01 -81.843  < 2e-16 ***
#   scale(invDom)            8.273e-03  5.866e-03  1.670e+02   1.410   0.1603    
# scale(ZipfFreq_Subtlex) -4.486e-02  6.741e-03  1.640e+02  -6.655 4.06e-10 ***
#   scale(OLD20)            -1.375e-02  7.713e-03  1.683e+02  -1.783   0.0764 .  
# scale(nsyl)              4.022e-03  1.100e-02  2.101e+03   0.365   0.7148    
# scale(nPhon)            -2.223e-03  9.641e-03  2.351e+04  -0.231   0.8176    
# scale(Length)            1.273e-02  7.546e-03  1.887e+02   1.687   0.0933 .  
# scale(famfRes)          -2.787e-02  6.581e-03  1.704e+02  -4.235 3.73e-05 ***
#   scale(pRT)               7.870e-02  2.171e-03  2.685e+04  36.247  < 2e-16 ***
#   scale(pACC)              1.043e-02  2.096e-03  2.682e+04   4.979 6.43e-07 ***
#   scale(Trial)            -1.681e-02  1.070e-02  2.240e+02  -1.571   0.1177    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.elp.rt.clean, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
armstrong.rt.m.p = lmer(invRT ~ 
                          #Polysemy
                          AmbiguityType +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) + 
                          #random effects
                          (1|SubjID) + (1|Word), 
                        pSubset)

summary(armstrong.rt.m.p)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)             -1.902e+00  2.186e-02  8.877e+01 -86.983  < 2e-16 ***
#   AmbiguityTypePolysemy   -2.650e-02  9.826e-03  1.744e+02  -2.696 0.007695 ** 
#   scale(ZipfFreq_Subtlex) -3.009e-02  5.680e-03  1.719e+02  -5.298 3.55e-07 ***
#   scale(OLD20)            -6.504e-03  6.836e-03  1.732e+02  -0.951 0.342755    
# scale(nsyl)              8.260e-03  1.640e-02  1.012e+04   0.504 0.614564    
# scale(nPhon)            -4.562e-03  1.573e-02  2.871e+04  -0.290 0.771751    
# scale(Length)            5.953e-03  6.906e-03  1.788e+02   0.862 0.389808    
# scale(famfRes)          -2.200e-02  5.775e-03  1.766e+02  -3.809 0.000192 ***
#   scale(pRT)               8.070e-02  2.002e-03  3.032e+04  40.300  < 2e-16 ***
#   scale(pACC)              8.074e-03  1.934e-03  3.028e+04   4.175 2.99e-05 ***
#   scale(Trial)             3.825e-04  8.532e-03  2.108e+02   0.045 0.964280    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Results
#Unambiguous - Homonymy (b<0.01,SE<0.01,t=1.41,p=0.16)
#Unambiguous - Polysemy (b=-0.02,SE=0.01,t=-2.70,p<0.01)

#Print table with models summary
tab_model(armstrong.rt.m.p, armstrong.rt.m.h, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "Armstrong & Plaut (2016) Replication Analysis on ELP Lexical Decision Data")

#Plot
figFolder = paste(wd, "S1_SemanticDiversity_LSA/Exp2_Ambiguity&SemD/Figures/", sep="")
armstrong.elp.rt.clean$AmbiguityType <- factor(armstrong.elp.rt.clean$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
plotSubset = subset(armstrong.elp.rt.clean, eDom_biggest<65 | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous")
armstrong.elp.rt.mean <- aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.elp.rt.mean, c("RT"), c("rt.mean"))
armstrong.elp.rt.se = aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.elp.rt.se, c("RT"), c("rt.se"))
armstrong.elp.rt.ambiguity <- merge(armstrong.elp.rt.mean, armstrong.elp.rt.se)
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_ELP_RT.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=4.5, width=5);
armstrong.elp.rt.p = ggplot(armstrong.elp.rt.ambiguity,
                            aes(x=AmbiguityType, y=rt.mean, fill=AmbiguityType, color=AmbiguityType)) +
  geom_bar(stat="identity", position=position_dodge(), color= "white", alpha = 0.5) +
  geom_errorbar(aes(ymin=rt.mean-rt.se, ymax=rt.mean+rt.se), width=.2,
                position=position_dodge(.9), color="black") + 
  ggtitle("ELP Lexical Decision") +
  theme(
    plot.title = element_text(size=20, hjust = 0.5),
    axis.text=element_text(size=15), axis.title=element_text(size=15), legend.position = "none") +
  coord_cartesian(ylim=c(500,650)) +
  ylab("Mean Reaction Time (ms)") +
  xlab("Condition") +
  scale_fill_manual(values=c("#0C5BB0FF", "#EE0011FF","#15983DFF")) + 
  geom_signif(y_position=c(570), xmin=c(1.8), xmax=c(3.2),
              annotation=c("**"), tip_length=0, size = 0.7, textsize = 12, color = "black")
grid.arrange(armstrong.elp.rt.p, nrow=1, ncol=1)
dev.off()

library(RGraphics)
library(gridExtra)

#Plot replication agaist original data 
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_RT_Comparison.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=6, width=14);
t0 <- textGrob("", gp=gpar(fontsize=28), just="centre")
t1 <- textGrob("Original Data", gp=gpar(fontsize=28), just="centre")
t2 <- textGrob("Replication Analysis", gp=gpar(fontsize=28), just="centre")
l1 = linesGrob(y = c(1,1),x = c(0.16, .98),  gp = gpar(col = "black", lwd = 2.5)) 
l2 = linesGrob(y = c(1,1),x = c(0.01, .98),  gp = gpar(col = "black", lwd = 2.5))
grid.arrange(arrangeGrob(t0, t1, t2, ncol=3, widths=c(0.15/3, 1/3.4, 2/3)),
             arrangeGrob(l1, l2, ncol=2, widths=c(1.15/3, 2/3)),
             arrangeGrob(armstrongResults.rt, 
             armstrong.blp.rt.p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
             armstrong.elp.rt.p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
             widths = c(1/3, 1/3.4, 1/3.4), ncol=3), 
             nrow=3, heights = c(1/8, 0.5/8, 6.8/8))
dev.off()
##########################

# Accuracy Analysis #
##########################

armstrong.elp.acc = armstrong.elp.clean
summary(armstrong.elp.acc)

#Set polysemes and unambiguous word to have a single completely dominant interpretation (biggest=100)
armstrong.elp.acc[["eDom_biggest"]][is.na(armstrong.elp.acc[["eDom_biggest"]])] = 100
describeBy(armstrong.elp.acc$eDom_biggest, group = armstrong.elp.acc$AmbiguityType)

#Modelling followeing Armstorng & Plaut (2016) procedure
armstrong.elp.acc$AmbiguityType <- relevel(armstrong.elp.acc$AmbiguityType, ref = "Unambiguous")
#get the inverse of dominant frequency to have higher values for balanced homonymous
armstrong.elp.acc$invDom = 1000/(armstrong.elp.acc$eDom_biggest)
# UNAMBIGUOUS - HOMONYMS
hSubset = subset(armstrong.elp.acc, AmbiguityType=="Unambiguous" | AmbiguityType=="Homonymy")
armstrong.acc.m.h = glmer(ACC ~ 
                            #Homonymy
                            scale(invDom) +
                            #control variable
                            scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes) + 
                            #scale(LogFreq_Subtlex) + scale(OLD) + scale(nSyll) + scale(nPhon) + scale(nLet) + scale(famfRes) +
                            #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                            #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                            scale(pRT) + scale(pACC) + scale(Trial) +  #lexicality of previous trial removed
                            #random effects
                            (1|SubjID) + (1|Word), 
                          hSubset,
                          family="binomial", control=glmerControl(optimizer="bobyqa",calc.derivs=F))
summary(armstrong.acc.m.h)
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              14.69759    4.54865   3.231  0.00123 ** 
#   scale(invDom)            -0.11188    0.28499  -0.393  0.69463    
# scale(ZipfFreq_Subtlex)   0.89042    0.34367   2.591  0.00957 ** 
#   scale(OLD20)             -0.27243    0.38653  -0.705  0.48093    
# scale(nsyl)              -4.23698  446.88240  -0.009  0.99244    
# scale(nPhon)              3.81406  453.17809   0.008  0.99328    
# scale(Length)             1.23820    0.43062   2.875  0.00404 ** 
#   scale(famfRes)            0.95564    0.32978   2.898  0.00376 ** 
#   scale(pRT)                0.26875    0.08759   3.068  0.00215 ** 
#   scale(pACC)               0.14952    0.08437   1.772  0.07638 .  
# scale(Trial)             -0.52812    0.09445  -5.592 2.25e-08 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.elp.acc, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
armstrong.acc.m.p = glmer(ACC ~ 
                            #Polysemy
                            AmbiguityType +
                            #control variable
                            scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes) + 
                            #scale(LogFreq_Subtlex) + scale(OLD) + scale(nSyll) + scale(nPhon) + scale(nLet) + scale(famfRes) +
                            #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                            #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                            scale(pRT) + scale(pACC) + scale(Trial) +  #lexicality of previous trial removed
                            #random effects
                            (1|SubjID) + (1|Word), 
                          pSubset,
                          family="binomial", control=glmerControl(optimizer="bobyqa",calc.derivs=F))
summary(armstrong.acc.m.p)
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             11.17862    0.59294  18.853  < 2e-16 ***
#   AmbiguityTypePolysemy    0.34080    0.27889   1.222  0.22171    
# scale(ZipfFreq_Subtlex)  0.25029    0.16492   1.518  0.12910    
# scale(OLD20)            -0.03596    0.18999  -0.189  0.84988    
# scale(nsyl)              0.23747    0.29580   0.803  0.42209    
# scale(nPhon)            -0.02734    0.25859  -0.106  0.91580    
# scale(Length)            0.07049    0.19186   0.367  0.71331    
# scale(famfRes)           0.47420    0.16577   2.861  0.00423 ** 
#   scale(pRT)               0.28564    0.05317   5.372 7.78e-08 ***
#   scale(pACC)              0.20170    0.04706   4.286 1.82e-05 ***
#   scale(Trial)             0.00471    0.05298   0.089  0.92916    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Results
#Unambiguous - Homonymy (b=-0.11,SE=0.28,t=-0.39,p=0.69)
#Unambiguous - Polysemy (b=0.34,SE=0.27,t=1.22,p=0.22)

#Print table with models summary
tab_model(armstrong.acc.m.p, armstrong.acc.m.h, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "Armstrong & Plaut (2016) Replication Analysis on ELP Lexical Decision Data")

#Plot
armstrong.elp.acc$AmbiguityType <- factor(armstrong.elp.acc$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
plotSubset = subset(armstrong.elp.acc, eDom_biggest<65 | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous")
armstrong.elp.acc.mean <- aggregate(ACC ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.elp.acc.mean, c("ACC"), c("acc.mean"))
armstrong.elp.acc.se = aggregate(ACC ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.elp.acc.se, c("ACC"), c("acc.se"))
armstrong.elp.acc.ambiguity <- merge(armstrong.elp.acc.mean, armstrong.elp.acc.se)
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_ELP_ACC.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=5);
armstrong.elp.acc.p = ggplot(armstrong.elp.acc.ambiguity,
                             aes(x=AmbiguityType, y=acc.mean, fill=AmbiguityType, color=AmbiguityType)) +
  geom_bar(stat="identity", position=position_dodge(), color= "white", alpha = 0.5) +
  geom_errorbar(aes(ymin=acc.mean-acc.se, ymax=acc.mean+acc.se), width=.2,
                position=position_dodge(.9), color="black") + 
  #facet_wrap(~ Contrast + nwDiff)+ 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") + #
  coord_cartesian(ylim=c(0.85,1)) +
  ylab("Mean Accuracy") +
  xlab("Condition") +
  scale_fill_manual(values=c("#0C5BB0FF", "#EE0011FF","#15983DFF"))
grid.arrange(armstrong.elp.acc.p, nrow=1, ncol=1)
dev.off()
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_ELP_ACC&RT.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=4.5, width=9);
grid.arrange(armstrong.elp.acc.p, armstrong.elp.rt.p, nrow=1, ncol=2, 
             top=textGrob("ELP Lexical Decision", gp=gpar(fontsize=20)))
dev.off()
##########################

# Semantic Diversity Analysis #
###############################

# Get data 
armstrong.semd = merge(armstrongStimuli, lexicalVars[c("Word","SemD")], by="Word", all.x = TRUE)
armstrong.semd = subset(armstrong.semd, AmbiguityType!="Hybrid")
summary(armstrong.semd)

#Modelling
armstrong.semd$AmbiguityType <- relevel(armstrong.semd$AmbiguityType, ref = "Unambiguous")
#Set polysemes and unambiguous word to have a single completely dominant interpretation (biggest=100)
armstrong.semd[["eDom_biggest"]][is.na(armstrong.semd[["eDom_biggest"]])] = 100
#get the inverse of dominant frequency to have higher values for balanced homonymous
armstrong.semd$invDom = 1000/(armstrong.semd$eDom_biggest)

# UNAMBIGUOUS - HOMONYMS
hSubset = subset(armstrong.semd, AmbiguityType=="Unambiguous" | AmbiguityType=="Homonymy")
armstrong.semd.m.h = lm(SemD ~ scale(invDom) +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes), 
                        hSubset)
summary(armstrong.semd.m.h)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              2.178710   0.008868 245.669  < 2e-16 ***
#   scale(invDom)           -0.017966   0.009283  -1.935   0.0533 .  
# scale(ZipfFreq_Subtlex)  0.137839   0.010412  13.238  < 2e-16 ***
#   scale(OLD20)             0.016066   0.012003   1.339   0.1811    
# scale(nsyl)             -0.007529   0.042119  -0.179   0.8582    
# scale(nPhon)             0.020781   0.041605   0.499   0.6176    
# scale(Length)           -0.059520   0.012351  -4.819 1.75e-06 ***
#   scale(famfRes)          -0.075858   0.010705  -7.086 3.19e-12 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.2437 on 748 degrees of freedom
# (10 observations deleted due to missingness)
# Multiple R-squared:  0.2325,	Adjusted R-squared:  0.2253 
# F-statistic: 32.37 on 7 and 748 DF,  p-value: < 2.2e-16

# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.semd, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
armstrong.semd.m.p = lm(SemD ~ AmbiguityType +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(nPhon) + scale(Length) + scale(famfRes),  
                        pSubset)
summary(armstrong.semd.m.p)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              2.202894   0.012188 180.743  < 2e-16 ***
#   AmbiguityTypePolysemy   -0.005934   0.016019  -0.370    0.711    
# scale(ZipfFreq_Subtlex)  0.126433   0.009466  13.356  < 2e-16 ***
#   scale(OLD20)             0.012642   0.010975   1.152    0.250    
# scale(nsyl)              0.004063   0.063483   0.064    0.949    
# scale(nPhon)             0.016264   0.063304   0.257    0.797    
# scale(Length)           -0.044524   0.011185  -3.981 7.47e-05 ***
#   scale(famfRes)          -0.113147   0.009643 -11.734  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.2237 on 832 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.2199,	Adjusted R-squared:  0.2133 
# F-statistic:  33.5 on 7 and 832 DF,  p-value: < 2.2e-16

#Results
#Unambiguous - Homonymy (b=-0.02,SE<0.01,t=-1.93,p=0.53)
#Unambiguous - Polysemy (b<0.01,SE=0.01,t=-0.37,p=0.71)

#Print table with models summary
tab_model(armstrong.semd.m.p, armstrong.semd.m.h, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "Armstrong & Plaut (2016) Simulation")

#Plot
figFolder = paste(wd, "S1_SemanticDiversity/Exp2_Ambiguity&SemD/Figures/", sep="")
armstrong.semd$AmbiguityType <- factor(armstrong.semd$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
plotSubset = subset(armstrong.semd, eDom_biggest<65 | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous")
armstrong.semd.mean <- aggregate(SemD ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.semd.mean, c("SemD"), c("semd.mean"))
armstrong.semd.se = aggregate(SemD ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.semd.se, c("SemD"), c("semd.se"))
armstrong.semd.ambiguity <- merge(armstrong.semd.mean, armstrong.semd.se)
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_SemD.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=4.5, width=5);
armstrong.semd.p = ggplot(armstrong.semd.ambiguity,
                             aes(x=AmbiguityType, y=semd.mean, fill=AmbiguityType, color=AmbiguityType)) +
  geom_bar(stat="identity", position=position_dodge(), color= "white", alpha = 0.5) +
  geom_errorbar(aes(ymin=semd.mean-semd.se, ymax=semd.mean+semd.se), width=.2,
                position=position_dodge(.9), color="black") + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") + #
  coord_cartesian(ylim=c(0,3)) +
  ylab("Semantic Diversity") +
  xlab("Condition") +
  scale_fill_manual(values=c("#0C5BB0FF", "#EE0011FF","#15983DFF"))
grid.arrange(armstrong.semd.p, nrow=1, ncol=1)
dev.off()
###############################


##############################################################################################################################################

