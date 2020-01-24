
# A1 - Validation Analysis

# Statistical analysis of the newly computed semantic diversity measures on behavioural data from the following megastudies 
# British Lexicon Project (BLP; Keuleers et al., 2012)
# English Lexicon Project (ELP; Balota et al., 2007)


# Preparation #
###############

#Clean workspace
rm(list=ls(all=T))

#Set working directory
setwd(dirname(getwd()))
getwd()

#Import used libraries
library(MASS)
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
library(car)
library(dplyr)
library(plyr)
library(sjmisc)
library(effects)

inverse <- function(x) {-1000/x};

# Read in Lexical Variables
lexicalVars = read.csv("stimuli/lexicalVariables.csv")
lexicalVars$Length = nchar(as.character(lexicalVars$Word))
setnames(lexicalVars, c("SemD_P2_W"), c("SemD"))
lexicalVars = lexicalVars[,c("Word","SemD","Freq","Length")]
summary(lexicalVars)
str(lexicalVars)
###############

###############################
# BLP - Lexical Decision Data #
###############################

# Read in & simple preprocessing
################################
blp.data = read.table("data/***", sep="\t", header=TRUE) #to download BLP trial data, save in /data and add file name
str(blp.data)
setnames(blp.data, c("spelling","participant","trial","accuracy","rt.raw", "previous.rt", "previous.accuracy"), 
         c("Word","SubjID","Trial","ACC","RT","pRT","pACC")) #RT outliers already removed (see Keuleers et al., 2012)
nrow(blp.data) # 2,240,940 observations
length(unique(blp.data$Word)) # 55,865 word types
#Set SubjID as factor
blp.data$SubjID = as.factor(blp.data$SubjID)
#Add previous trial data to include in the modelling to eliminate auto-correlations that violate model assumptions 
#(for discussion, see Baayen & Milin, 2010; Barr et al., 2013; Bolker et al., 2009).
blp.data$lexicality = mapvalues(blp.data$lexicality, from=c("N","W"), to=c(0,1))
blp.data$lexicality  = as.numeric(as.character(blp.data$lexicality))
blp.data$pLexicality = as.numeric(unlist(append(list(NA), head(blp.data$lexicality, -1))))
#Select words only
blp.data = subset(blp.data, lexicality==1)
nrow(blp.data) # 1,120,470 observations
length(unique(blp.data$Word)) # 28,730 word types
#Subset only columns of interest (rt col: outliers and inaccurate responses set to NA)
blp.data = blp.data[c("Word","SubjID","Trial","ACC","RT","pRT","pACC","pLexicality")]
#Add lexical variables 
blp.data = merge(blp.data, lexicalVars, by="Word") # inner join
summary(blp.data)

# Descriptives
nrow(blp.data) # 453,044 observations
length(unique(blp.data$Word)) # number of word types (11,618)
length(unique(blp.data$SubjID)) # number of subjects (78)
blp.data %>% group_by(SubjID) %>% summarise_each(list(n_distinct), Trial) #number of trials (5852 odd and 5791 even subjs)
describe(blp.data$RT) #(M = 609ms, SD = 188ms)
describe(blp.data$ACC) #(M = 86%, SD = 0.35)
summary(blp.data)
################################

# Reaction Time Analysis #
##########################

# Remove inaccurate responses
blp.data.rt = subset(blp.data, ACC==1)
blp.data.rt = subset(blp.data.rt, RT>0)
summary(blp.data.rt)

# Check normality of residuals following Box-Cox procedure (transform data by X^lambda; Box and Cox, 1964) 
hist(blp.data.rt$RT)
qqnorm(blp.data.rt$RT)
qqline(blp.data.rt$RT)
boxcox(blp.data.rt$RT ~ 1) #lambda ~ 0, log transfomation of RT is recommended
hist(log(blp.data.rt$RT))
qqnorm(log(blp.data.rt$RT))
qqline(log(blp.data.rt$RT))

# Outliers trimming
#Generate base model
rt.m0 = lmer(log(RT) ~ 1 + scale(Trial) + (1|SubjID) + (1|Word), blp.data.rt)
summary(rt.m0)
hist(resid(rt.m0))
#Remove if residuals > 2.5SDs
blp.data.rt.clean = blp.data.rt[abs(scale(resid(rt.m0)))<2.5,]
((nrow(blp.data.rt)-nrow(blp.data.rt.clean))/nrow(blp.data.rt))*100 # 2.7% data removed
dotplot(ranef(rt.m0, condVar=T))
summary(blp.data.rt.clean)

# Modelling
blp.ld.rt.m <- lmer(log(RT) ~ scale(SemD)*scale(Freq) + 
                      scale(Length) +
                      scale(pRT) + scale(pACC) + scale(pLexicality) + 
                      scale(Trial) +
                     (1|SubjID) + (1|Word), 
                    blp.data.rt.clean)
summary(blp.ld.rt.m)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)              6.386e+00  1.253e-02  7.745e+01 509.751  < 2e-16 ***
#   scale(SemD)             -2.059e-03  7.502e-04  1.053e+04  -2.745  0.00606 ** 
#   scale(Freq)             -7.656e-02  7.915e-04  1.101e+04 -96.723  < 2e-16 ***
#   scale(Length)            1.148e-03  7.383e-04  1.068e+04   1.556  0.11982    
# scale(pRT)               3.299e-02  3.274e-04  3.713e+05 100.764  < 2e-16 ***
#   scale(pACC)              1.653e-04  3.250e-04  3.713e+05   0.509  0.61092    
# scale(pLexicality)      -6.336e-04  4.477e-04  3.712e+05  -1.415  0.15699    
# scale(Trial)            -4.254e-03  1.414e-03  1.057e+04  -3.007  0.00264 ** 
#   scale(SemD):scale(Freq)  2.796e-03  7.092e-04  1.102e+04   3.943 8.11e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Print table with model summary
tab_model(blp.ld.rt.m, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "BLP Lexical Decision Data")

#Plot effects
blp.rt.FxS.e = as.data.frame(effect(term="scale(SemD):scale(Freq)", mod=blp.ld.rt.m, xlevels=list(SemD=10, Freq=6)))
blp.rt.FxS.e[,c('fit','lower','upper')] <- exp(blp.rt.FxS.e[,c('fit','lower','upper')])
blp.rt.FxS.p = ggplot(blp.rt.FxS.e, aes(x=SemD, y=fit, group=Freq, color=Freq, fill=Freq)) + 
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3, colour = NA) + 
  labs(x= "Semantic Diversity", y="Reaction Time (ms)", color="Frequency", fill="Frequency") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  scale_color_continuous(high = "#132B43", low = "#56B1F7") + 
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.title = element_text(size=30), legend.position = "none") +
  ggtitle("BLP Lexical Decision")
blp.rt.FxS.p
ggsave("figures/A1_ValidationAnalysis_BLP_RT_Interaction_FreqxSemD.png", blp.rt.FxS.p, 
       width = 7, height = 9, dpi = 300, type = "cairo");

##########################

# Accurary Analysis #
#####################

blp.data.acc = blp.data
blp.acc.m <- glmer(ACC ~ scale(SemD)*scale(Freq) +
                     scale(Length) +
                     scale(pRT) + scale(pACC) + scale(pLexicality) + 
                     scale(Trial) +
                    (1|SubjID) + (1|Word), 
                    blp.data.acc, 
                    family="binomial", control=glmerControl(optimizer="bobyqa",calc.derivs=F))
summary(blp.acc.m)
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)              2.972217   0.077558  38.322  < 2e-16 ***
#   scale(SemD)              0.068960   0.016001   4.310 1.63e-05 ***
#   scale(Freq)              1.500007   0.017907  83.768  < 2e-16 ***
#   scale(Length)            0.657733   0.015702  41.889  < 2e-16 ***
#   scale(pRT)               0.143198   0.006780  21.122  < 2e-16 ***
#   scale(pACC)              0.022375   0.005302   4.220 2.44e-05 ***
#   scale(pLexicality)      -0.012288   0.011595  -1.060    0.289    
#   scale(Trial)             0.008717   0.028359   0.307    0.759    
#   scale(SemD):scale(Freq) -0.077865   0.016049  -4.852 1.22e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Print table with model summary
tab_model(blp.acc.m, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "BLP Lexical Decision Data")
###################




###############################
# ELP - Lexical Decision Data #
###############################

#Clean workspace
rm(list=setdiff(ls(), c("wd","figFolder","lexicalVars","blp.rt.FxS.p")))

# Read in & simple preprocessing
################################
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
#remove strage negative RT
elp.ld.data = subset(elp.ld.data, RT>0)
#Add lexical variables
elp.ld.data = merge(elp.ld.data, lexicalVars, by="Word") # inner join

# Descriptives
length(unique(elp.ld.data$Word)) #number of words (15,688)
length(unique(elp.ld.data$SubjID)) #number of subjects (785)
elp.ld.data %>% group_by(SubjID) %>% summarise_each(list(n_distinct), Trial) #number of trials (variable)
describe(elp.ld.data$RT) #(M = 770ms, SD = 380ms)
describe(elp.ld.data$ACC) #(M = 88%, SD = 0.33)
summary(elp.ld.data)
################################

# Reaction Time Analysis #
##########################

# Remove inaccurate responses
elp.ld.rt = subset(elp.ld.data, ACC==1)
summary(elp.ld.rt)

# Check normality of residuals following Box-Cox procedure (transform data by X^lambda; Box and Cox, 1964) 
hist(elp.ld.rt$RT)
qqnorm(elp.ld.rt$RT)
qqline(elp.ld.rt$RT)
boxcox(elp.ld.rt$RT ~ 1) #lambda ~ -1, log RT is recommended
hist(log(elp.ld.rt$RT))
qqnorm(log(elp.ld.rt$RT))
qqline(log(elp.ld.rt$RT))

# Outliers trimming
#Generate base model
rt.m0 = lmer(log(RT) ~ 1 + scale(Trial) + (1|SubjID) + (1|Word), elp.ld.rt)
summary(rt.m0)
hist(resid(rt.m0))
#Remove if residuals > 2.5SDs
elp.ld.rt.clean = elp.ld.rt[abs(scale(resid(rt.m0)))<2.5,]
((nrow(elp.ld.rt)-nrow(elp.ld.rt.clean))/nrow(elp.ld.rt))*100 # 2.4% data removed
dotplot(ranef(rt.m0, condVar=T))
summary(elp.ld.rt.clean)

# Modelling
elp.ld.rt.m <- lmer(log(RT) ~ scale(SemD)*scale(Freq) + 
                      scale(Length) +
                      scale(pRT) + scale(pACC) + scale(pLexicality) + 
                      scale(Trial) +
                      (1|SubjID) + (1|Word), 
                    elp.ld.rt.clean)
summary(elp.ld.rt.m)
# Fixed effects:
#   Estimate Std. Error         df  t value Pr(>|t|)    
# (Intercept)              6.551e+00  6.082e-03  7.986e+02 1077.050  < 2e-16 ***
#   scale(SemD)             -7.615e-03  7.350e-04  1.761e+04  -10.359  < 2e-16 ***
#   scale(Freq)             -7.491e-02  7.793e-04  1.806e+04  -96.127  < 2e-16 ***
#   scale(Length)            5.596e-02  7.524e-04  1.793e+04   74.373  < 2e-16 ***
#   scale(pRT)               5.394e-02  3.761e-04  5.155e+05  143.432  < 2e-16 ***
#   scale(pACC)              1.450e-03  3.329e-04  5.150e+05    4.355 1.33e-05 ***
#   scale(pLexicality)      -1.922e-02  3.292e-04  5.156e+05  -58.377  < 2e-16 ***
#   scale(Trial)            -1.112e-02  3.267e-04  5.150e+05  -34.022  < 2e-16 ***
#   scale(SemD):scale(Freq)  1.449e-03  6.962e-04  1.794e+04    2.081   0.0374 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#Print table with model summary
tab_model(elp.ld.rt.m, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "ELP Lexical Decision Data")

#Plot Effects
elp.rt.FxS.e = as.data.frame(effect(term="scale(SemD):scale(Freq)", mod=elp.ld.rt.m, xlevels=list(SemD=10, Freq=6)))
elp.rt.FxS.e[,c('fit','lower','upper')] <- exp(elp.rt.FxS.e[,c('fit','lower','upper')])
elp.ld.rt.FxS.p = ggplot(elp.rt.FxS.e, aes(x=SemD, y=fit, group=Freq, color=Freq, fill=Freq)) + 
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3, colour = NA) + 
  labs(x= "Semantic Diversity", y="Reaction Time (ms)", color="Frequency", fill="Frequency") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  scale_color_continuous(high = "#132B43", low = "#56B1F7") + 
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.title = element_text(size=30), legend.position = "none") +
  ggtitle("ELP Lexical Decision")
elp.ld.rt.FxS.p
ggsave("figures/A1_ValidationAnalysis_ELP_LD_RT_Interaction_FreqxSemD.png", elp.ld.rt.FxS.p, 
       width = 7, height = 9, dpi = 300, type = "cairo");
##########################

# Accurary Analysis #
#####################

elp.ld.data.acc = elp.ld.data

elp.ld.acc.m <- glmer(ACC ~ scale(SemD)*scale(Freq) +
                     scale(Length) +
                     scale(pRT) + scale(pACC) + scale(pLexicality) + 
                     scale(Trial) +
                     (1|SubjID) + (1|Word), 
                   elp.ld.data.acc, 
                   family="binomial", control=glmerControl(optimizer="bobyqa",calc.derivs=F))
summary(elp.ld.acc.m)
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              2.906398   0.031254  92.992  < 2e-16 ***
#   scale(SemD)              0.118817   0.011325  10.491  < 2e-16 ***
#   scale(Freq)              1.001845   0.012395  80.826  < 2e-16 ***
#   scale(Length)            0.407736   0.011390  35.797  < 2e-16 ***
#   scale(pRT)               0.091688   0.005196  17.644  < 2e-16 ***
#   scale(pACC)              0.017532   0.004453   3.937 8.25e-05 ***
#   scale(pLexicality)       0.140174   0.004687  29.908  < 2e-16 ***
#   scale(Trial)            -0.037114   0.004639  -8.000 1.24e-15 ***
#   scale(SemD):scale(Freq) -0.032294   0.011088  -2.913  0.00358 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Print table with model summary
tab_model(elp.ld.acc.m, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "ELP Lexical Decision Data")
#####################



#####################
# ELP - Naming Data #
#####################

#Clean workspace
rm(list=setdiff(ls(), c("wd","figFolder","lexicalVars","blp.rt.FxS.p","elp.ld.rt.FxS.p")))


# Read in & simple preprocessing
################################
# ELP NMG trial-level data
elpNMGTrials = read.table("data/**", sep=",", header=TRUE) #to download ELP trial data, save in /data and add file name
#Rename columns
setnames(elpNMGTrials, c("Stimulus","Participant", "Accuracy", "NMG_RT"), c("Word","SubjID", "ACC", "RT"))
str(elpNMGTrials)
elpNMGTrials$Trial = as.numeric(as.character(elpNMGTrials$Trial))
elpNMGTrials$RT = as.numeric(as.character(elpNMGTrials$RT))
#Add previous trial data
elpNMGTrials$ACC[elpNMGTrials$ACC > 1] = 0
elpNMGTrials$pACC = as.numeric(unlist(append(list(NA), head(elpNMGTrials$ACC, -1))))
elpNMGTrials$pRT = as.numeric(unlist(append(list(NA), head(elpNMGTrials$RT, -1))))
summary(elpNMGTrials)
str(elpNMGTrials)
#Subset only columns of interest
elpNMGTrials = elpNMGTrials[ , -which(names(elpNMGTrials) %in% c("CodingRT"))]
#Look up the number of words (40,449)
length(unique(elpNMGTrials$Word)) 
#Add lexical variables
elp.nmg.data = merge(elpNMGTrials, lexicalVars, by="Word") # inner join

# Descriptives
length(unique(elp.nmg.data$Word)) #number of words (15,688)
length(unique(elp.nmg.data$SubjID)) #number of subjects (405)
elp.nmg.data %>% group_by(SubjID) %>% summarise_each(list(n_distinct), Trial) #number of trials (variable)
elp.nmg.data$Trial = as.numeric(as.character(elp.nmg.data$Trial))
elp.nmg.data$RT = as.numeric(as.character(elp.nmg.data$RT))
describe(elp.nmg.data$RT) #(M = 764ms, SD = 416ms)
describe(elp.nmg.data$ACC) #(M = 94%, SD = 0.25)
summary(elp.nmg.data)
str(elp.nmg.data)
################################


# Reaction Time Analysis #
##########################

# Remove inaccurate responses
elp.nmg.rt = subset(elp.nmg.data, ACC==1)
elp.nmg.rt = elp.nmg.rt[!(is.na(elp.nmg.rt$RT)),]
summary(elp.nmg.rt)

# Check normality of residuals following Box-Cox procedure (transform data by X^lambda; Box and Cox, 1964) 
hist(elp.nmg.rt$RT)
qqnorm(elp.nmg.rt$RT)
qqline(elp.nmg.rt$RT)
boxcox(elp.nmg.rt$RT ~ 1) #lambda ~ 0, log RT is recommended
hist(log(elp.nmg.rt$RT))
qqnorm(log(elp.nmg.rt$RT))
qqline(log(elp.nmg.rt$RT))

# Outliers trimming
#Generate base model
rt.m0 = lmer(log(RT) ~ 1 + scale(Trial) + (1|SubjID) + (1|Word), elp.nmg.rt)
summary(rt.m0)
hist(resid(rt.m0))
#Remove if residuals > 2.5SDs
elp.nmg.rt.clean = elp.nmg.rt[abs(scale(resid(rt.m0)))<2.5,]
(nrow(elp.nmg.rt)-nrow(elp.nmg.rt.clean))/nrow(elp.nmg.rt) # 2.1% data removed
dotplot(ranef(rt.m0, condVar=T))
summary(elp.nmg.rt.clean)

# Modelling
elp.nmg.rt.m <- lmer(log(RT) ~ scale(SemD)*scale(Freq) + 
                      scale(Length) +
                      scale(pRT) + scale(pACC) +  
                      scale(Trial) +
                      (1|SubjID) + (1|Word), 
                    elp.nmg.rt.clean)
summary(elp.nmg.rt.m)
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)              6.533e+00  6.939e-03  4.104e+02 941.491  < 2e-16 ***
#   scale(SemD)             -3.626e-03  7.096e-04  1.802e+04  -5.111 3.24e-07 ***
#   scale(Freq)             -4.955e-02  7.491e-04  1.814e+04 -66.149  < 2e-16 ***
#   scale(Length)            5.245e-02  7.207e-04  1.811e+04  72.773  < 2e-16 ***
#   scale(pRT)               2.652e-02  3.206e-04  4.173e+05  82.725  < 2e-16 ***
#   scale(pACC)              8.705e-04  3.062e-04  4.172e+05   2.843  0.00447 ** 
#   scale(Trial)            -8.679e-03  2.759e-04  4.176e+05 -31.460  < 2e-16 ***
#   scale(SemD):scale(Freq)  2.092e-03  6.711e-04  1.821e+04   3.117  0.00183 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Print table with model summary
tab_model(elp.nmg.rt.m, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "ELP Naming Data")

#Plot Effects
elp.rt.FxS.e = as.data.frame(effect(term="scale(SemD):scale(Freq)", mod=elp.nmg.rt.m, xlevels=list(SemD=10, Freq=6)))
elp.rt.FxS.e[,c('fit','lower','upper')] <- exp(elp.rt.FxS.e[,c('fit','lower','upper')])
elp.nmg.rt.FxS.p = ggplot(elp.rt.FxS.e, aes(x=SemD, y=fit, group=Freq, color=Freq, fill=Freq)) + 
  geom_point() + 
  geom_line(size=1.4) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3, colour = NA) + 
  labs(x= "Semantic Diversity", y="Reaction Time (ms)", color="Frequency", fill="Frequency") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  scale_color_continuous(high = "#132B43", low = "#56B1F7") + 
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25), plot.title = element_text(size=30), 
        legend.title=element_text(size=20), legend.text=element_text(size=12), legend.position = c(0.72, 0.93), legend.direction = "horizontal") +
  ggtitle("ELP Naming")
elp.nmg.rt.FxS.p
ggsave("A1_ValidationAnalysis_ELP_NMG_RT_Interaction_FreqxSemD.png", elp.nmg.rt.FxS.p, 
       width = 7, height = 9, dpi = 300, type = "cairo");
##########################


# Accurary Analysis #
#####################

elp.nmg.data.acc = elp.nmg.data

elp.nmg.acc.m <- glmer(ACC ~ scale(SemD)*scale(Freq) +
                        scale(Length) +
                        scale(pRT) + scale(pACC) + 
                        scale(Trial) +
                        (1|SubjID) + (1|Word), 
                      elp.nmg.data.acc, 
                      family="binomial", control=glmerControl(optimizer="bobyqa",calc.derivs=F))
summary(elp.nmg.acc.m)
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              3.688522   0.039079  94.386  < 2e-16 ***
#   scale(SemD)              0.020663   0.013050   1.583   0.1133    
# scale(Freq)              0.744456   0.014117  52.735  < 2e-16 ***
#   scale(Length)           -0.006147   0.012678  -0.485   0.6278    
# scale(pRT)              -0.107995   0.006240 -17.307  < 2e-16 ***
#   scale(pACC)              0.050434   0.006413   7.864 3.73e-15 ***
#   scale(Trial)             0.407896   0.023313  17.497  < 2e-16 ***
#   scale(SemD):scale(Freq) -0.030396   0.012611  -2.410   0.0159 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Print table with model summary
tab_model(elp.nmg.acc.m, p.style = "both", show.intercept = FALSE, show.se = TRUE, show.ci = FALSE,
          show.stat = TRUE, show.df = TRUE, string.se = "SE", string.stat = "t", 
          string.est = "b", title = "ELP Naming Data")
#####################


