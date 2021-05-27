
# A2 Simulation Analysis

# Replication analysis of Rodd et al. (2002) and Armstrong & Plaut (2016) on megastudies and 
# simulation analysis on the newly semantic diversity measures


###############
# Preparation #
###############

#Clean workspace
rm(list=ls(all=T))

#Set working directory
setwd(dirname(getwd()))
getwd()

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
library(tidyr)
library(yarrr)
library(car)
library(plyr)
library(ggsignif)
library(effects)
library(RGraphics)
library(gridExtra)

inverse <- function(x) {-1000/x};

# Read in & preprocessing #
###########################
# Read in Lexical Variables
lexicalVars = read.csv("stimuli/lexicalVariables.csv")
lexicalVars$Length = nchar(as.character(lexicalVars$Word))
setnames(lexicalVars, c("SemD_P2_W"), c("SemD"))
lexicalVars = lexicalVars[,c("Word","SemD","Freq","Length")]
summary(lexicalVars)
str(lexicalVars)

# BLP data
blp.data = read.table("data/***", sep="\t", header=TRUE) #to download BLP trial data, save in /data and add file name
blp.ld.data$lexicality = mapvalues(blp.ld.data$lexicality, from=c("N","W"), to=c(0,1))
blp.ld.data$lexicality  = as.numeric(as.character(blp.ld.data$lexicality))
#Add previous trial data
blp.ld.data$pLexicality = as.numeric(unlist(append(list(NA), head(blp.ld.data$lexicality, -1))))
#Select words only
blp.ld.data = subset(blp.ld.data, lexicality==1) 
#Subset only columns of interest (rt col: outliers and inaccurate responses set to NA)
blp.ld.data = blp.ld.data[c("spelling","participant","block","trial","accuracy","rt.raw","previous.rt","previous.accuracy","pLexicality")]
#Rename columns
setnames(blp.ld.data, c("spelling","participant","block","trial","accuracy","rt.raw", "previous.rt", "previous.accuracy"), 
         c("Word","SubjID","Block","Trial","ACC","RT","pRT","pACC"))
#Look up the number of words (28,730)
length(unique(blp.ld.data$Word)) 
summary(blp.ld.data)

#ELP data
elp.ld.data = read.table("data/**", sep=",", header=TRUE) #to download ELP trial data, save in /data and add file name
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
#remove strange negative RT
elp.ld.data = subset(elp.ld.data, RT>0)
#Look up the number of words (40,481)
length(unique(elp.ld.data$Word)) 
summary(elp.ld.data)
###########################

# Replication & Simulation Analysis of Rodd et al. (2001) Experiment 1 (Visual Lexical Decision - Regression Design)
#####################################################################################################################

# Get stimuli & properties
###############################
roddStimuli_Exp1 = read.table("stimuli/Rodd2002_Exp1_Stimuli.csv", sep=",", header=TRUE)
roddStimuli_Exp1 = merge(roddStimuli_Exp1, lexicalVars, by="Word", all.x = TRUE)
str(roddStimuli_Exp1)
summary(roddStimuli_Exp1)
###############################

# BLP - Reaction Time Analysis #
################################

# Get behavioural data from BLP
rodd.blp = merge(roddStimuli_Exp1, blp.ld.data, by="Word", all.x = TRUE)

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

#Plot
m.plot = lm(RT ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M), rodd.blp.rt)
rodd.exp1.blp.rt.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=m.plot, xlevels=list(nSenses_WMD=2)))
rodd.exp1.blp.rt.meanings.e = as.data.frame(effect(term="Meanings", mod=m.plot))
rodd.exp1.blp.rt.meanings.e$Meanings <- relevel(rodd.exp1.blp.rt.meanings.e$Meanings, ref = "One")
cairo_ps('figures\A2_SimulationAnalysis_Rodd2001_Exp1_BLP_RT_effects.eps',
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
################################

# ELP - Reaction Time Analysis #
################################

# Get behavioural data from ELP
rodd.elp = merge(roddStimuli_Exp1, elp.ld.data, by="Word", all.x = TRUE)
rodd.elp$SubjID = as.factor(rodd.elp$SubjID)
length(unique(rodd.elp$Word)) #182 w
summary(rodd.elp)
str(rodd.elp)

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


#Plot
m.plot = lm(RT ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
            rodd.elp.rt)
rodd.exp1.elp.rt.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=m.plot, xlevels=list(nSenses_WMD=2)))
rodd.exp1.elp.rt.meanings.e = as.data.frame(effect(term="Meanings", mod=m.plot))
rodd.exp1.elp.rt.meanings.e$Meanings <- relevel(rodd.exp1.elp.rt.meanings.e$Meanings, ref = "One")
cairo_ps('figures\A2_SimulationAnalysis_Rodd2001_Exp1_ELP_RT_effects.eps',
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
################################

# Semantic Diversity Simulation #
#################################
# Get data 
rodd.exp1.semd = roddStimuli_Exp1
summary(rodd.exp1.semd)
rodd.exp1.semd$Meanings <- relevel(rodd.exp1.semd$Meanings, ref = "One")
rodd.exp1.semd.m = lm(SemD ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
                      rodd.exp1.semd)
summary(rodd.exp1.semd.m)

#Plot
rodd.exp1.semd.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=rodd.exp1.semd.m, xlevels=list(nSenses_WMD=2)))
rodd.exp1.semd.meanings.e = as.data.frame(effect(term="Meanings", mod=rodd.exp1.semd.m))
rodd.exp1.semd.meanings.e$Meanings <- relevel(rodd.exp1.semd.meanings.e$Meanings, ref = "One")
cairo_ps('figures\A2_SimulationAnalysis_Rodd2001_Exp1_SemD_effects.eps',
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
#################################


# Replication & Simulation Analysis of Rodd et al. (2001) Experiment 2 (Visual Lexical Decision - Factorial Design) 
####################################################################################################################

# Get stimuli & properties
###############################
roddStimuli_Exp2 = read.table("stimuli/Rodd2002_Exp2_Stimuli.csv", sep=",", header=TRUE)
roddStimuli_Exp2 = merge(roddStimuli_Exp2, lexicalVars, by="Word", all.x = TRUE)
str(roddStimuli_Exp2)
summary(roddStimuli_Exp2)
###############################

# BLP #
# Reaction Time Analysis #
##########################

# Get behavioural data from BLP
rodd.blp = merge(roddStimuli_Exp2, blp.ld.data, by="Word", all.x = TRUE)
rodd.blp$SubjID = as.factor(rodd.blp$SubjID)
length(unique(rodd.blp$Word)) #128 w
summary(rodd.blp)
str(rodd.blp)

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

#BySubj ANOVA
rodd.rt.subj <- aggregate(invRT ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp.rt)
rodd.rt.subj.aov = aov(invRT ~  Senses*Meanings + Freq + Length, data = rodd.rt.subj)
summary(rodd.rt.subj.aov)

describeBy(rodd.blp.rt$RT, group = rodd.blp.rt$Senses) 
diff(as.matrix(aggregate(rodd.blp.rt$RT, by=list(rodd.blp.rt$Senses), FUN=mean)[2])) #15ms difference
describeBy(rodd.blp.rt$RT, group = rodd.blp.rt$Meanings) 
diff(as.matrix(aggregate(rodd.blp.rt$RT, by=list(rodd.blp.rt$Meanings), FUN=mean)[2])) #8ms difference

#UNTRANSFORMED RT
#ByItem ANOVA
rodd.rt.item <- aggregate(RT ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp.rt)
rodd.rt.item.aov = aov(RT ~ Senses*Meanings + Freq + Length, data = rodd.rt.item)
summary(rodd.rt.item.aov)

#BySubj ANOVA
rodd.rt.subj <- aggregate(RT ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.blp.rt)
rodd.rt.subj.aov = aov(RT ~  Senses*Meanings + Freq + Length, data = rodd.rt.subj)
summary(rodd.rt.subj.aov)


#Plot
rodd.rt.item = aggregate(RT ~ Word + Senses + Meanings, FUN=mean, rodd.blp.rt)
rodd.rt.item$Meanings <- relevel(rodd.rt.item$Meanings, ref = "One")
cairo_ps('figures\A2_SimulationAnalysis_Rodd2001_Exp2_BLP_RT.eps', 
         family='sans', onefile=T, antialias='default', height=6, width=8);
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


# ELP #
# Reaction Time Analysis #
##########################

# Get behavioural data from BLP
rodd.elp = merge(roddStimuli_Exp2, elp.ld.data, by="Word", all.x = TRUE)
rodd.elp$SubjID = as.factor(rodd.elp$SubjID)
length(unique(rodd.elp$Word)) #128 w
summary(rodd.elp)
str(rodd.elp)

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

#BySubj ANOVA
rodd.rt.subj <- aggregate(logRT ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp.rt)
rodd.rt.subj.aov = aov(logRT ~  Senses*Meanings + Freq + Length, data = rodd.rt.subj)
summary(rodd.rt.subj.aov)


describeBy(rodd.elp.rt$RT, group = rodd.elp.rt$Senses) 
diff(as.matrix(aggregate(rodd.elp.rt$RT, by=list(rodd.elp.rt$Senses), FUN=mean)[2])) #16ms difference
describeBy(rodd.elp.rt$RT, group = rodd.elp.rt$Meanings) 
diff(as.matrix(aggregate(rodd.elp.rt$RT, by=list(rodd.elp.rt$Meanings), FUN=mean)[2])) #8ms difference


#UNTRANSFORMED RT
#ByItem ANOVA
rodd.rt.item <- aggregate(RT ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp.rt)
rodd.rt.item.aov = aov(RT ~ Senses*Meanings + Freq + Length, data = rodd.rt.item)
summary(rodd.rt.item.aov)

#BySubj ANOVA
rodd.rt.subj <- aggregate(RT ~ SubjID + Senses + Meanings + Freq + Length, FUN=mean, rodd.elp.rt)
rodd.rt.subj.aov = aov(RT ~  Senses*Meanings + Freq + Length, data = rodd.rt.subj)
summary(rodd.rt.subj.aov)


#Plot
rodd.rt.item = aggregate(RT ~ Word + Senses + Meanings, FUN=mean, rodd.elp.rt)
rodd.rt.item$Meanings <- relevel(rodd.rt.item$Meanings, ref = "One")
cairo_ps('figures\A2_SimulationAnalysis_Rodd2001_Exp2_ELP_RT.eps', 
         family='sans', onefile=T, antialias='default', height=6, width=8);
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


# Semantic Diversity Simulation #
#################################

# Get data 
rodd.semd = roddStimuli_Exp2
summary(rodd.semd)

#ByItem ANOVA
rodd.semd.item <- aggregate(SemD ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.semd)
rodd.semd.item.aov = aov(SemD ~ Senses*Meanings + Freq + Length, data = rodd.semd.item)
summary(rodd.semd.item.aov)

describeBy(rodd.semd.item$SemD, group = rodd.semd.item$Senses)
diff(as.matrix(aggregate(rodd.semd.item$SemD, by=list(rodd.semd.item$Senses), FUN=mean)[2])) #0.0002 difference
describeBy(rodd.semd.item$SemD, group = rodd.semd.item$Meanings) 
diff(as.matrix(aggregate(rodd.semd.item$SemD, by=list(rodd.semd.item$Meanings), FUN=mean)[2])) #0.0260 difference

#Plot
rodd.semd.item$Meanings <- relevel(rodd.semd.item$Meanings, ref = "One")
cairo_ps('figures/A2_SimulationAnalysis_Rodd2001_Exp2_SemD.eps', 
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
rm(list=setdiff(ls(), c("lexicalVars","blp.ld.data", "elp.ld.data")))

# Get stimuli & properties
###############################
armstrongStimuli = read.table("stimuli/Armstrong2016_Stimuli.csv", sep=",", header=TRUE)
armstrongStimuli = merge(armstrongStimuli, lexicalVars, by="Word", all.x = TRUE)
str(armstrongStimuli)
summary(armstrongStimuli)
###############################

#Plot Results from Armstrong & Plaut (2016)
###########################################
armstrongResults = read.table("stimuli/Armstrong2016_Results", sep=",", header=TRUE)
cairo_ps('figures/A2_SimulationAnalysis_Armstrong&Plaut2016_Results_RT.eps', 
         family='sans', onefile=T, antialias='default', height=4.5, width=5);
armstrongResults.rt = ggplot(armstrongResults,
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
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) +  #lexicality of previous trial removed
                          #random effects
                          (1|SubjID) + (1|Word), 
                        hSubset)
summary(armstrong.rt.m.h)


# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.blp.rt.clean, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
armstrong.rt.m.p = lmer(invRT ~ 
                          #Polysemy
                          AmbiguityType +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) +  #lexicality of previous trial removed
                          #random effects
                          (1|SubjID) + (1|Word), 
                        pSubset)

summary(armstrong.rt.m.p)

#Plot
armstrong.blp.rt.clean$AmbiguityType <- factor(armstrong.blp.rt.clean$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
plotSubset = subset(armstrong.blp.rt.clean, eDom_biggest<65 | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous")
armstrong.blp.rt.mean <- aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.blp.rt.mean, c("RT"), c("rt.mean"))
armstrong.blp.rt.se = aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.blp.rt.se, c("RT"), c("rt.se"))
armstrong.blp.rt.ambiguity <- merge(armstrong.blp.rt.mean, armstrong.blp.rt.se)
cairo_ps('figures/A2_SimulationAnalysis_Armstrong&Plaut2016_BLP_RT.eps',
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
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) +  
                          #random effects
                          (1|SubjID) + (1|Word), 
                        hSubset)
summary(armstrong.rt.m.h)


# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.elp.rt.clean, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
armstrong.rt.m.p = lmer(invRT ~ 
                          #Polysemy
                          AmbiguityType +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes) + 
                          #previous trial data to eliminate auto-correlations that violate model assumptions (for discussion,
                          #see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
                          scale(pRT) + scale(pACC) + scale(Trial) + 
                          #random effects
                          (1|SubjID) + (1|Word), 
                        pSubset)

summary(armstrong.rt.m.p)


#Plot
armstrong.elp.rt.clean$AmbiguityType <- factor(armstrong.elp.rt.clean$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
plotSubset = subset(armstrong.elp.rt.clean, eDom_biggest<65 | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous")
armstrong.elp.rt.mean <- aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.elp.rt.mean, c("RT"), c("rt.mean"))
armstrong.elp.rt.se = aggregate(RT ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.elp.rt.se, c("RT"), c("rt.se"))
armstrong.elp.rt.ambiguity <- merge(armstrong.elp.rt.mean, armstrong.elp.rt.se)
cairo_ps('figures/A2_SimulationAnalysis_Armstrong&Plaut2016_ELP_RT.eps',
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

#Plot replication against original data 
cairo_ps('figures\A2_SimulationAnalysis_Armstrong&Plaut2016_RT_Comparison.eps',
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

# Semantic Diversity Analysis #
###############################

# Get data 
armstrong.semd = armstrongStimuli
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
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes), 
                        hSubset)
summary(armstrong.semd.m.h)

# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.semd, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
armstrong.semd.m.p = lm(SemD ~ AmbiguityType +
                          #control variable
                          scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes),  
                        pSubset)
summary(armstrong.semd.m.p)

#Plot
armstrong.semd$AmbiguityType <- factor(armstrong.semd$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
plotSubset = subset(armstrong.semd, eDom_biggest<65 | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous")
armstrong.semd.mean <- aggregate(SemD ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.semd.mean, c("SemD"), c("semd.mean"))
armstrong.semd.se = aggregate(SemD ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.semd.se, c("SemD"), c("semd.se"))
armstrong.semd.ambiguity <- merge(armstrong.semd.mean, armstrong.semd.se)
cairo_ps('figures\A2_SimulationAnalysis_Armstrong&Plaut2016_SemD.eps', 
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



#Additional analysis - Revisions
#############################################################################################################################################

# Read in Replicability Data
############################

repData = read.csv("data/SemanticDiversity/diversity_ReplicabilityData.csv")
#Add frequency estimates from SUBTLEX UK (Van Heuven et al., 2013)
subtlex = read.csv("Resources/Norms/Databases/SUBTLEX/SUBTLEX-UK.csv")
setnames(subtlex, c("Spelling","LogFreq.Zipf."), c("Word","Freq"))
repData = merge(repData, subtlex[,c("Word","Freq")], by="Word")
repData = merge(repData, lexicalVars[,c("Word","Length")], by="Word")
str(repData)
summary(repData)

repData$SemD_bnc = repData$SemD
repData$SemD_w1000 = repData$SemD
corpora = repData[,c("Word","SemD_bnc","SemD_uw","SemD_wp")]
corpora$rep = "Corpora"
corpora$rep <- as.factor(corpora$rep)
summary(corpora)
corporaLong = gather(corpora, condition, value, SemD_bnc:SemD_wp, factor_key=TRUE)
summary(corporaLong)
context = repData[,c("Word","SemD_w1000","SemD_w100")]
context$rep = "Context Length"
context$rep <- as.factor(context$rep)
summary(context)
contextLong = gather(context, condition, value, SemD_w1000:SemD_w100, factor_key=TRUE)
summary(contextLong)
repDataLong = rbind(corporaLong, contextLong)
summary(repDataLong)
figFolder = paste(wd, "figures/", sep="")

#plot
roddStimuli$Study = "Rodd et al. (2002)"
roddStimuli$Study <- as.factor(roddStimuli$Study)
summary(roddStimuli)
armstrongStimuli = armstrongStimuli[!duplicated(armstrongStimuli[ , c("Word")]),]
armstrongStimuli$Study = "Armstrong & Plaut (2016)"
armstrongStimuli$Study <- as.factor(armstrongStimuli$Study)
summary(armstrongStimuli)
stimuli = rbind(roddStimuli[,c("Word", "AmbiguityType", "Study")],armstrongStimuli[,c("Word", "AmbiguityType", "Study")])
summary(stimuli)
semdAcross = merge(stimuli, repDataLong, by="Word", all.x = TRUE)
semdAcross = subset(semdAcross, AmbiguityType!="Hybrid")
summary(semdAcross)


semdAcross$AmbiguityType <- factor(semdAcross$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
semdAcross$condition = revalue(semdAcross$condition, c("SemD_w100"="Short (100w)", 
                                                       "SemD_w1000"="Long (1000w)",
                                                       "SemD_bnc"="BNC",
                                                       "SemD_uw"="ukWaC",
                                                       "SemD_wp"="Wackypedia"))
semdAcross$Corpora = semdAcross$condition
semdAcross$ContextLength = semdAcross$condition
summary(semdAcross)

semdAcross = semdAcross[complete.cases(semdAcross), ]


cairo_ps(paste(figFolder,'Replication&Simulation_Replicability_AcrossCorpora&ContextLength.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=8, width=9);
f1 = ggplot(subset(semdAcross, rep=="Corpora"), aes(Corpora, value, fill=AmbiguityType)) +
  geom_boxplot(outlier.shape=NA)+
  facet_wrap(~Study) + 
  geom_point(position=position_jitterdodge(), alpha=1/10) +
  labs(y = "Semantic Diversity")+ 
  theme(axis.text=element_text(size=15), 
        axis.title=element_text(size=15), 
        strip.text.x = element_text(size=15),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
        legend.position="none") +
  ylim(1,3.5)
f2 = ggplot(subset(semdAcross, rep=="Context Length"), aes(ContextLength, value, fill=AmbiguityType)) +
  geom_boxplot(outlier.shape=NA)+
  facet_wrap(~Study)+ 
  geom_point(position=position_jitterdodge(), alpha=1/10) +
  labs(y = "Semantic Diversity") + 
  theme(axis.text=element_text(size=15), 
        axis.title=element_text(size=15), 
        strip.text.x = element_text(size=15),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
        legend.position="bottom")+
  ylim(1,3.5)
grid.arrange(f1, f2, nrow = 2, 
             heights = c(0.463,0.537),
             top=textGrob("Semantic Diverisity by Lexical Ambiguity across \n Corpora and Context Length", gp=gpar(fontsize=20)))
dev.off()


cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp2_SemD_wp.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=8);
par(mfrow=c(2,1))
pirateplot(formula = value ~ AmbiguityType + Corpora + Study,
           data = subset(semdAcross, rep=="Corpora"),
           ylab = "Semantic Diversity",
           ylim = c(0,3.5),
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

pirateplot(formula = value ~ AmbiguityType + ContextLength + Study,
           data = subset(semdAcross, rep=="Context Length"),
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
############################

# Rodd et al. (2002) EXP 1 #
############################

# Get data 
rodd.exp1.rep = merge(roddStimuli, repData, by="Word", all.x = TRUE)
#number of meanings/senses
rodd.exp1.rep = merge(rodd.exp1.rep, wordsmythDict, by="Word", all.x = TRUE)
#lexical neighborhood
rodd.exp1.rep = merge(rodd.exp1.rep, blpStimuli[c("Word","coltN")], by="Word", all.x = TRUE)
#Concreteness
rodd.exp1.rep = merge(rodd.exp1.rep, concr[c("Word","Conc.M")], by="Word", all.x = TRUE)
summary(rodd.exp1.rep)

rodd.exp1.rep$Meanings <- relevel(rodd.exp1.rep$Meanings, ref = "One")
rodd.exp1.rep.ww.m = lm(SemD_w100 ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
                        rodd.exp1.rep)
summary(rodd.exp1.rep.ww.m)

rodd.exp1.rep.wp.m = lm(SemD_wp ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
                        rodd.exp1.rep)
summary(rodd.exp1.rep.wp.m)

rodd.exp1.rep.uw.m = lm(SemD_uw ~ Meanings + scale(nSenses_WMD) + scale(coltN) + scale(Freq) + scale(Length) + scale(Conc.M),  
                        rodd.exp1.rep)
summary(rodd.exp1.rep.uw.m)

rodd.exp1.semd.rep.ww.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=rodd.exp1.rep.ww.m, xlevels=list(nSenses_WMD=2)))
rodd.exp1.semd.rep.ww.meanings.e = as.data.frame(effect(term="Meanings", mod=rodd.exp1.rep.ww.m))
rodd.exp1.semd.rep.ww.meanings.e$Meanings <- relevel(rodd.exp1.semd.rep.ww.meanings.e$Meanings, ref = "One")

cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_SemD_effects_w100.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=4.5, width=9);
s = ggplot(rodd.exp1.semd.rep.ww.senses.e, aes(x=nSenses_WMD, y=fit)) + 
  ylim(0,3.5)+
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3, colour = NA) + 
  labs(x= "Number of Senses", y="Semantic Diversity") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") 
m= ggplot(rodd.exp1.semd.rep.ww.meanings.e, aes(x=Meanings, y=fit, fill=Meanings)) + 
  ylim(0,3.5)+
  geom_bar(stat="identity", position="dodge", color="white", alpha = 0.5) + 
  geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.2, position=position_dodge(width=0.9))+
  labs(x= "Meanings", y="Semantic Diversity") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none")
grid.arrange(s, m, nrow=1, ncol=2,
             top=textGrob("Semantic Diversity", gp=gpar(fontsize=20)))
dev.off()


rodd.exp1.semd.rep.wp.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=rodd.exp1.rep.wp.m, xlevels=list(nSenses_WMD=2)))
rodd.exp1.semd.rep.wp.meanings.e = as.data.frame(effect(term="Meanings", mod=rodd.exp1.rep.wp.m))
rodd.exp1.semd.rep.wp.meanings.e$Meanings <- relevel(rodd.exp1.semd.rep.wp.meanings.e$Meanings, ref = "One")

cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_SemD_effects_wp.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=4.5, width=9);
s = ggplot(rodd.exp1.semd.rep.wp.senses.e, aes(x=nSenses_WMD, y=fit)) + 
  ylim(0,3.5)+
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3, colour = NA) + 
  labs(x= "Number of Senses", y="Semantic Diversity") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") 
m= ggplot(rodd.exp1.semd.rep.wp.meanings.e, aes(x=Meanings, y=fit, fill=Meanings)) + 
  ylim(0,3.5)+
  geom_bar(stat="identity", position="dodge", color="white", alpha = 0.5) + 
  geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.2, position=position_dodge(width=0.9))+
  labs(x= "Meanings", y="Semantic Diversity") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none")
grid.arrange(s, m, nrow=1, ncol=2,
             top=textGrob("Semantic Diversity", gp=gpar(fontsize=20)))
dev.off()


rodd.exp1.semd.rep.uw.senses.e = as.data.frame(effect(term="scale(nSenses_WMD)", mod=rodd.exp1.rep.uw.m, xlevels=list(nSenses_WMD=2)))
rodd.exp1.semd.rep.uw.meanings.e = as.data.frame(effect(term="Meanings", mod=rodd.exp1.rep.uw.m))
rodd.exp1.semd.rep.uw.meanings.e$Meanings <- relevel(rodd.exp1.semd.rep.uw.meanings.e$Meanings, ref = "One")

cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp1_SemD_effects_uw.eps',sep=''),
         family='sans', onefile=T, antialias='default', height=4.5, width=9);
s = ggplot(rodd.exp1.semd.rep.uw.senses.e, aes(x=nSenses_WMD, y=fit)) + 
  ylim(0,3.5)+
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3, colour = NA) + 
  labs(x= "Number of Senses", y="Semantic Diversity") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") 
m= ggplot(rodd.exp1.semd.rep.uw.meanings.e, aes(x=Meanings, y=fit, fill=Meanings)) + 
  ylim(0,3.5)+
  geom_bar(stat="identity", position="dodge", color="white", alpha = 0.5) + 
  geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.2, position=position_dodge(width=0.9))+
  labs(x= "Meanings", y="Semantic Diversity") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2F86FFFF", "#DE0012FF")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none")
grid.arrange(s, m, nrow=1, ncol=2,
             top=textGrob("Semantic Diversity", gp=gpar(fontsize=20)))
dev.off()
#############################

# Rodd et al. (2002) EXP 2 #
############################

# Get data 
rodd.semd.rep = merge(roddStimuli, repData, by="Word", all.x = TRUE)
summary(rodd.semd.rep)

setnames(rodd.semd.rep, c("SemD"), c("SemD_w1000"))
rodd.semd.rep$Diff = rodd.semd.rep$SemD_wp - rodd.semd.rep$SemD_w1000
par(mfrow=c(1,1))
rodd.semd.rep$AmbiguityType <- factor(rodd.semd.rep$AmbiguityType, levels = c("Unambiguous", "Polysemy","Homonymy"))
rodd.semd.rep = subset(rodd.semd.rep, AmbiguityType!="Hybrid")


####### Word Window 
#ByItem ANOVA
rodd.semd.rep.item <- aggregate(SemD_w100 ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.semd.rep)
rodd.semd.rep.item.aov = aov(SemD_w100 ~ Senses*Meanings + Freq + Length, data = rodd.semd.rep.item)
summary(rodd.semd.rep.item.aov)

describeBy(rodd.semd.rep.item$SemD_w100, group = rodd.semd.rep.item$Senses)
diff(as.matrix(aggregate(rodd.semd.rep.item$SemD_w100, by=list(rodd.semd.rep.item$Senses), FUN=mean)[2])) #<0.01 difference
describeBy(rodd.semd.rep.item$SemD_w100, group = rodd.semd.rep.item$Meanings) 
diff(as.matrix(aggregate(rodd.semd.rep.item$SemD_w100, by=list(rodd.semd.rep.item$Meanings), FUN=mean)[2])) #0.05 difference

#Plot
rodd.semd.rep.item$Meanings <- relevel(rodd.semd.rep.item$Meanings, ref = "One")
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp2_SemD_ContextLengthReplicability.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=8);
par(mfrow=c(1,2))
pirateplot(formula = SemD_w100 ~ Senses,
           data = rodd.semd.rep.item,
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
pirateplot(formula = SemD_w100 ~ Meanings,
           data = rodd.semd.rep.item,
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

####### Wackypedia
#ByItem ANOVA
rodd.semd.rep.item <- aggregate(SemD_wp ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.semd.rep)
rodd.semd.rep.item.aov = aov(SemD_wp ~ Senses*Meanings + Freq + Length, data = rodd.semd.rep.item)
summary(rodd.semd.rep.item.aov)

describeBy(rodd.semd.rep.item$SemD_wp, group = rodd.semd.rep.item$Senses)
diff(as.matrix(aggregate(rodd.semd.rep.item$SemD_wp, by=list(rodd.semd.rep.item$Senses), FUN=mean)[2])) #0.12 difference
describeBy(rodd.semd.rep.item$SemD_wp, group = rodd.semd.rep.item$Meanings) 
diff(as.matrix(aggregate(rodd.semd.rep.item$SemD_wp, by=list(rodd.semd.rep.item$Meanings), FUN=mean)[2])) #0.04 difference

#Plot
rodd.semd.rep.item$Meanings <- relevel(rodd.semd.rep.item$Meanings, ref = "One")
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp2_SemD_wp.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=8);
par(mfrow=c(1,2))
pirateplot(formula = SemD_wp ~ Senses,
           data = rodd.semd.rep.item,
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
pirateplot(formula = SemD_wp ~ Meanings,
           data = rodd.semd.rep.item,
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

####### UkWac 
#ByItem ANOVA
rodd.semd.rep.item <- aggregate(SemD_uw ~ Word + Senses + Meanings + Freq + Length, FUN=mean, rodd.semd.rep)
rodd.semd.rep.uw.aov = aov(SemD_uw ~ Senses*Meanings + Freq + Length, data = rodd.semd.rep.item)
summary(rodd.semd.rep.uw.aov)

describeBy(rodd.semd.rep.item$SemD_uw, group = rodd.semd.rep.item$Senses)
diff(as.matrix(aggregate(rodd.semd.rep.item$SemD_uw, by=list(rodd.semd.rep.item$Senses), FUN=mean)[2])) #0.04 difference
describeBy(rodd.semd.rep.item$SemD_uw, group = rodd.semd.rep.item$Meanings) 
diff(as.matrix(aggregate(rodd.semd.rep.item$SemD_uw, by=list(rodd.semd.rep.item$Meanings), FUN=mean)[2])) #0.03 difference

#Plot
rodd.semd.rep.item$Meanings <- relevel(rodd.semd.rep.item$Meanings, ref = "One")
cairo_ps(paste(figFolder,'Replication&Simulation_Rodd2001_Exp2_SemD_uw.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=6, width=8);
par(mfrow=c(1,2), title="UkWac")
pirateplot(formula = SemD_uw ~ Senses,
           data = rodd.semd.rep.item,
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
pirateplot(formula = SemD_uw ~ Meanings,
           data = rodd.semd.rep.item,
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
############################

# Armstrong & Plaut (2016) #
############################


# Get data 
armstrongStimuli = armstrongStimuli[!duplicated(armstrongStimuli[ , c("Word")]),]
armstrong.semd.rep = merge(armstrongStimuli, repData[c("Word","SemD","SemD_w100","SemD_wp","SemD_uw")], by="Word", all.y = TRUE)
armstrong.semd.rep = subset(armstrong.semd.rep, AmbiguityType!="Hybrid")
summary(armstrong.semd.rep)

#Modelling
armstrong.semd.rep$AmbiguityType <- relevel(armstrong.semd.rep$AmbiguityType, ref = "Unambiguous")
#Set polysemes and unambiguous word to have a single completely dominant interpretation (biggest=100)
armstrong.semd.rep[["eDom_biggest"]][is.na(armstrong.semd.rep[["eDom_biggest"]])] = 100
#get the inverse of dominant frequency to have higher values for balanced homonymous
armstrong.semd.rep$invDom = 1000/(armstrong.semd.rep$eDom_biggest)

# UNAMBIGUOUS - HOMONYMS
hSubset = subset(armstrong.semd.rep, AmbiguityType=="Unambiguous" | AmbiguityType=="Homonymy")
#####Word window
armstrong.semd.ww.m.h = lm(SemD_w100 ~ scale(invDom) +
                             #control variable
                             scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes), 
                           hSubset)
summary(armstrong.semd.ww.m.h)

#####Wackypedia
armstrong.semd.wp.m.h = lm(SemD_wp ~ scale(invDom) +
                             #control variable
                             scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes), 
                           hSubset)
summary(armstrong.semd.wp.m.h)

#####Word window
armstrong.semd.uw.m.h = lm(SemD_uw ~ scale(invDom) +
                             #control variable
                             scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes), 
                           hSubset)
summary(armstrong.semd.uw.m.h)

# UNAMBIGUOUS - POLYSEMY
pSubset = subset(armstrong.semd.rep, AmbiguityType=="Unambiguous" | AmbiguityType=="Polysemy")
###Word window
armstrong.semd.ww.m.p = lm(SemD_w100 ~ AmbiguityType +
                             #control variable
                             scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes),  
                           pSubset)
summary(armstrong.semd.ww.m.p)


###Wackypedia
armstrong.semd.wp.m.p = lm(SemD_wp ~ AmbiguityType +
                             #control variable
                             scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes),  
                           pSubset)
summary(armstrong.semd.wp.m.p)


###UkWac
armstrong.semd.uw.m.p = lm(SemD_uw ~ AmbiguityType +
                             #control variable
                             scale(ZipfFreq_Subtlex) + scale(OLD20) + scale(nsyl) + scale(Length) + scale(famfRes),  
                           pSubset)
summary(armstrong.semd.uw.m.p)



#Plot
armstrong.semd.rep$AmbiguityType <- factor(armstrong.semd.rep$AmbiguityType, levels = c("Homonymy", "Unambiguous", "Polysemy"))
plotSubset = subset(armstrong.semd.rep, eDom_biggest<65 | AmbiguityType=="Polysemy" | AmbiguityType=="Unambiguous")
###Context Length
armstrong.semd.mean <- aggregate(SemD_w100 ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.semd.mean, c("SemD_w100"), c("semd.mean"))
armstrong.semd.se = aggregate(SemD_w100 ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.semd.se, c("SemD_w100"), c("semd.se"))
armstrong.semd.ambiguity <- merge(armstrong.semd.mean, armstrong.semd.se)
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_SemD_ContextLengthReplicability.eps',sep=''), 
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
###Wakypedia
armstrong.semd.mean <- aggregate(SemD_wp ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.semd.mean, c("SemD_wp"), c("semd.mean"))
armstrong.semd.se = aggregate(SemD_wp ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.semd.se, c("SemD_wp"), c("semd.se"))
armstrong.semd.ambiguity <- merge(armstrong.semd.mean, armstrong.semd.se)
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_SemD_wp.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=4.5, width=5);
armstrong.semd.wp.p = ggplot(armstrong.semd.ambiguity,
                             aes(x=AmbiguityType, y=semd.mean, fill=AmbiguityType, color=AmbiguityType)) +
  geom_bar(stat="identity", position=position_dodge(), color= "white", alpha = 0.5) +
  geom_errorbar(aes(ymin=semd.mean-semd.se, ymax=semd.mean+semd.se), width=.2,
                position=position_dodge(.9), color="black") + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") + #
  coord_cartesian(ylim=c(0,3)) +
  ylab("Semantic Diversity") +
  xlab("Condition") +
  ggtitle("Wackypedia") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("#0C5BB0FF", "#EE0011FF","#15983DFF"))
grid.arrange(armstrong.semd.wp.p, nrow=1, ncol=1)
dev.off()
###UkWac
armstrong.semd.mean <- aggregate(SemD_uw ~ AmbiguityType, plotSubset, function(x) c(mean = mean(x)))
setnames(armstrong.semd.mean, c("SemD_uw"), c("semd.mean"))
armstrong.semd.se = aggregate(SemD_uw ~ AmbiguityType, plotSubset, function(x) c(se = sd(x)/sqrt(length(x))))
setnames(armstrong.semd.se, c("SemD_uw"), c("semd.se"))
armstrong.semd.ambiguity <- merge(armstrong.semd.mean, armstrong.semd.se)
cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_SemD_uw.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=4.5, width=5);
armstrong.semd.uw.p = ggplot(armstrong.semd.ambiguity,
                             aes(x=AmbiguityType, y=semd.mean, fill=AmbiguityType, color=AmbiguityType)) +
  geom_bar(stat="identity", position=position_dodge(), color= "white", alpha = 0.5) +
  geom_errorbar(aes(ymin=semd.mean-semd.se, ymax=semd.mean+semd.se), width=.2,
                position=position_dodge(.9), color="black") + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15),legend.position = "none") + #
  coord_cartesian(ylim=c(0,3)) +
  ylab("Semantic Diversity") +
  xlab("Condition") +
  ggtitle("UkWac") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("#0C5BB0FF", "#EE0011FF","#15983DFF"))
grid.arrange(armstrong.semd.uw.p, nrow=1, ncol=1)
dev.off()

cairo_ps(paste(figFolder,'Replication&Simulation_Armstrong&Plaut2016_SemD_CorporaReplicability.eps',sep=''), 
         family='sans', onefile=T, antialias='default', height=4.5, width=10);
grid.arrange(armstrong.semd.wp.p, armstrong.semd.uw.p, nrow=1, ncol=2)
dev.off()

par(mfrow=c(2,2))
pirateplot(formula = SemD_w100 ~ AmbiguityType,
           data = armstrong.semd.rep,
           ylab = "Semantic Diversity",
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           main = "BNC (1,000-WordWindow)",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
pirateplot(formula = SemD_w100 ~ AmbiguityType,
           data = armstrong.semd.rep,
           ylab = "Semantic Diversity",
           theme = 2, 
           inf.method = "se",# Start with theme 2
           inf.f.col = "black",
           inf.b.col = "black",
           main = "BNC (100-WordWindow)",
           #inf.f.o = 0, # Turn off inf fill
           #inf.b.o = 0, # Turn off inf border
           point.o = .2,   # Turn up points
           bar.f.o = .5, # Turn up bars
           bean.f.o = .4, # Light bean filling
           bean.b.o = .2, # Light bean border
           point.col = "black")
pirateplot(formula = SemD_wp ~ AmbiguityType,
           data = armstrong.semd.rep,
           ylab = "Semantic Diversity",
           theme = 2, 
           main = "Wackypedia",
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
pirateplot(formula = SemD_uw ~ AmbiguityType,
           data = armstrong.semd.rep,
           ylab = "Semantic Diversity",
           main = "UkWac",
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


###############################

##############################################################################################################################################

