# Vowel Harmony experiment in Finnish
# single presentation
# Lexical Decision Experiment (Experiment 3)

# libraries (we may not use all of them)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(multcomp)
# library(tidybayes)
library(jtools)
library(ggpubr)
library(brms)

library(multcomp)
library(ggmcmc)
library(effects)
library(see)
library(dplyr)
library(plyr)


library(BayesFactor)
library(bayestestR)

# Lexical Decision Experiment (Experiment 3)

# open file LDT (set working directory to the right folder)
LDT_VH=read.csv("LME_SinglePres_HARMONY_data.csv", header=T, sep=";")
head(LDT_VH)

# Analyses on the correct RT data for word targets
byTrial <- LDT_VH %>% 
  filter(lexicality=="nw")

head(byTrial)

# Averages RTs and accuracy
tapply(byTrial$Trcorrect,list(byTrial$DISHAR),mean, na.rm=T)
tapply(byTrial$accuracy,list(byTrial$DISHAR),mean, na.rm=T)


# definining the factors
byTrial$subject=as.factor(byTrial$subject)
byTrial$item=as.factor(byTrial$item)


# Linear mixed effect model (RTs)
LDTnonword_LME = lmer(-1000/Trcorrect ~ DISHARc + (1+DISHARc|item) + (1+DISHARc|subject), data = byTrial, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(LDTnonword_LME)



# Linear mixed effect model (accuracy)
LDTnonword_ACC_LME = glmer(accuracy ~ DISHARc + (1+DISHARc|item) + (1|subject), data = byTrial, family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(LDTnonword_ACC_LME)


#BAYES FACTORS
#BF01 a vs b

byTrialab <- byTrial %>% 
  filter(Trcorrect !="")


byTrialab$iTrcorrect=-1000/(byTrialab$Trcorrect)

byTrialab$DISHAR=as.factor(byTrialab$DISHAR)
byTrialab$subject=as.factor(byTrialab$subject)
byTrialab$item=as.factor(byTrialab$item)



bf_FULL = generalTestBF(iTrcorrect ~ DISHAR + subject*DISHAR+item,
                        data = byTrialab, whichRandom = c("subject", "item") )%>% sort()
bf_FULL

#BF01
sort(bf_FULL)[9]/sort(bf_FULL)[7]















