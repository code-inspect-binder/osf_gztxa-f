# Vowel Harmony experiments in Finnish
# using masked priming (Perea, Hyönä, and Marcet)
# Pre-registered study

# Lexical Decision Experiment (Experiment 1)
# Naming Experiment (Experiment 2)

# libraries (we may not need all of them)
library(tidyverse)
library(lmerTest)
library(jtools)
library(ggpubr)
library(ggmcmc)
library(effects)
library(see)
library(dplyr)
library(plyr)
library(BayesFactor)


# Lexical Decision Experiment (Experiment 1)

# open file LDT (set working directory to the right folder)
LDT_VH=read.csv("LDT_masked_priming_VH.csv", header=T, sep=";")
head(LDT_VH)

# Analyses on the correct RT data for word targets
byTrial <- LDT_VH %>% 
  filter(lexicality=="WORD")

head(byTrial)

# Averages RTs and accuracy
# abcde--order of conditions in Table 1 (i.e., ID=a , VSNH=b...)
tapply(byTrial$RT,list(byTrial$prime),mean, na.rm=T)
tapply(byTrial$RT,list(byTrial$primec),mean, na.rm=T)
tapply(byTrial$accuracy,list(byTrial$prime),mean, na.rm=T)
tapply(byTrial$accuracy,list(byTrial$primec),mean, na.rm=T)


# definining the factors
byTrial$primec=as.factor(byTrial$primec)
byTrial$subject=as.factor(byTrial$subject)
byTrial$item=as.factor(byTrial$item)

# List of pre-registered contrasts
c1 = rbind(c(1, -1, 0, 0, 0), c(0, 1, -1, 0, 0), c(0, 0, 1, -1, 0), c(0, 0, 0, 1, -1))
cc1=ginv(c1)

# Linear mixed effect model (RTs)
# -1000/RT as pre-registered
LDTword_LME = lmer(-1000/RT ~ primec + (1|item) + (1+primec|subject), data = byTrial, contrasts = list(primec = cc1), control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(LDTword_LME)



# Exploratory: Linear mixed effect model (RTs) with AuthorTest_Finnish
LDTword_LME_ARF = lmer(-1000/RT ~ primec*AR_Finnish + (1|item) + (1|subject), data = byTrial, contrasts = list(primec = cc1), control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(LDTword_LME_ARF)
plot(allEffects(LDTword_LME_ARF), x.var = "AR_Finnish")

# Exploratory: Linear mixed effect model (RTs) with AuthorTest
LDTword_LME_ARE = lmer(-1000/RT ~ primec*AR_English + (1|item) + (1|subject), data = byTrial, contrasts = list(primec = cc1), control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(LDTword_LME_ARE)
plot(allEffects(LDTword_LME_ARE), x.var = "AR_English")


# Density Plot (a vs. b; ID vs. VS_dysharmonious)
byTrialab = subset(byTrial, primec == "a" | primec == "b")
p = ggplot(byTrialab, aes(x = RT, color = prime)) +
  geom_density()
p + scale_fill_brewer(palette="Dark2") + theme_apa()


#BAYES FACTORS
#BF10 a vs b
byTrialab <- byTrial %>% 
  filter(primec=="a" | primec=="b") %>% 
  filter(RT !="")

byTrialab$primec=as.factor(byTrialab$primec)
byTrialab$subject=as.factor(byTrialab$subject)
byTrialab$item=as.factor(byTrialab$item)

bf_FULLab = generalTestBF(RTtransf ~ primec + subject*primec+item,
                        data = byTrialab, whichRandom = c("subject", "item") )%>% sort()
bf_FULLab

#BF10
#Note that this value has an associated "% error"
sort(bf_FULLab)[9]/sort(bf_FULLab)[8]


#BF10 b vs c
byTrialbc <- byTrial %>% 
  filter(primec=="b" | primec=="c") %>% 
  filter(RT !="")

byTrialbc$primec=as.factor(byTrialbc$primec)
byTrialbc$subject=as.factor(byTrialbc$subject)
byTrialbc$item=as.factor(byTrialbc$item)

bf_FULLbc = generalTestBF(RTtransf ~ primec + subject*primec+item,
                          data = byTrialbc, whichRandom = c("subject", "item") )%>% sort()
bf_FULLbc

#BF10
sort(bf_FULLbc)[9]/sort(bf_FULLbc)[8]

#BF10 c vs d
byTrialcd <- byTrial %>% 
  filter(primec=="c" | primec=="d") %>% 
  filter(RT !="")

byTrialcd$primec=as.factor(byTrialcd$primec)
byTrialcd$subject=as.factor(byTrialcd$subject)
byTrialcd$item=as.factor(byTrialcd$item)

bf_FULLcd = generalTestBF(RTtransf ~ primec + subject*primec+item,
                          data = byTrialcd, whichRandom = c("subject", "item") )%>% sort()
bf_FULLcd

#BF10
sort(bf_FULLcd)[8]/sort(bf_FULLcd)[9]


#BF10 d vs e
byTrialde <- byTrial %>% 
  filter(primec=="d" | primec=="e") %>% 
  filter(RT !="")

byTrialde$primec=as.factor(byTrialde$primec)
byTrialde$subject=as.factor(byTrialde$subject)
byTrialde$item=as.factor(byTrialde$item)

bf_FULLde = generalTestBF(RTtransf ~ primec + subject*primec+item,
                          data = byTrialde, whichRandom = c("subject", "item") )%>% sort()
bf_FULLde

#BF10
sort(bf_FULLde)[8]/sort(bf_FULLde)[9]



# Generalized linear mixed effect model (accuracy)
LDTword_GLME_accuracy = glmer(accuracy ~ primec + (1|item) + (1|subject), data = byTrial, contrasts = list(primec = cc1), family = binomial)
summary(LDTword_GLME_accuracy)



# Analyses on the NONWORD targets LDT
byTrialnw <- LDT_VH %>% 
  filter(lexicality=="NW")

# Averages RTs and accuracy
tapply(byTrialnw$RT,list(byTrialnw$primec),mean, na.rm=T)
tapply(byTrialnw$accuracy,list(byTrialnw$primec),mean, na.rm=T)

# definining the factors
byTrialnw$primec=as.factor(byTrialnw$primec)
byTrialnw$subject=as.factor(byTrialnw$subject)
byTrialnw$item=as.factor(byTrialnw$item)

# Linear mixed effect model (RTs) nonwords
LDTword_LMEnw = lmer(-1000/RT ~ primec + (1|item) + (1+primec|subject), data = byTrialnw, contrasts = list(primec = cc1), control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(LDTword_LMEnw)

# Generalized linear mixed effect model (accuracy) nonwords
LDTword_GLME_accuracynw = glmer(accuracy ~ primec + (1|item) + (1+primec|subject), data = byTrialnw, contrasts = list(primec = cc1), family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(LDTword_GLME_accuracynw)





# NAMING experiment (Experiment 2)
# open file (set working directory to the right folder)
N_VH=read.csv("Naming_masked_priming_VH.csv", header=T, sep=";")
head(N_VH)


# Analyses on the correct RT data for word targets (unnecessary the filter "RT < 200" here, already done in .csv)
nbyTrial <- N_VH %>% 
  filter(RT > 200)

# Averages RTs and accuracy
tapply(nbyTrial$RT,list(nbyTrial$prime),mean, na.rm=T)
tapply(nbyTrial$RT,list(nbyTrial$primec),mean, na.rm=T)


# definining the factors
nbyTrial$primec=as.factor(nbyTrial$primec)
nbyTrial$subject=as.factor(nbyTrial$subject)
nbyTrial$item=as.factor(nbyTrial$item)

# List of contrasts
c1 = rbind(c(1, -1, 0, 0, 0), c(0, 1, -1, 0, 0), c(0, 0, 1, -1, 0), c(0, 0, 0, 1, -1))
cc1=ginv(c1)

# Linear mixed effect model (RTs)
Nword_LME = lmer(-1000/RT ~ primec + (1|item) + (1+primec|subject), data = nbyTrial, contrasts = list(primec = cc1), control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(Nword_LME)


# Linear mixed effect model (RTs) with ART_English
Nword_LME_AR = lmer(-1000/RT ~ primec*AR_English + (1|item) + (1+primec|subject), data = nbyTrial, contrasts = list(primec = cc1), control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(Nword_LME_AR)
plot(allEffects(Nword_LME_AR), x.var = "AR_English")



#BAYES FACTORS
#BF10 a vs b
byTrialab <- nbyTrial %>% 
  filter(primec=="a" | primec=="b") %>% 
  filter(RT !="")

byTrialab$primec=as.factor(byTrialab$primec)
byTrialab$subject=as.factor(byTrialab$subject)
byTrialab$item=as.factor(byTrialab$item)

bf_FULLab = generalTestBF(RTtransf ~ primec + subject*primec+item,
                          data = byTrialab, whichRandom = c("subject", "item") )%>% sort()
bf_FULLab
#BF10
sort(bf_FULLab)[8]/sort(bf_FULLab)[9]


#BF10 b vs c
byTrialbc <- nbyTrial %>% 
  filter(primec=="b" | primec=="c") %>% 
  filter(RT !="")

byTrialbc$primec=as.factor(byTrialbc$primec)
byTrialbc$subject=as.factor(byTrialbc$subject)
byTrialbc$item=as.factor(byTrialbc$item)

bf_FULLbc = generalTestBF(RTtransf ~ primec + subject*primec+item,
                          data = byTrialbc, whichRandom = c("subject", "item") )%>% sort()
bf_FULLbc
#BF10
sort(bf_FULLbc)[9]/sort(bf_FULLbc)[8]


#BF10 c vs d
byTrialcd <- nbyTrial %>% 
  filter(primec=="c" | primec=="d") %>% 
  filter(RT !="")

byTrialcd$primec=as.factor(byTrialcd$primec)
byTrialcd$subject=as.factor(byTrialcd$subject)
byTrialcd$item=as.factor(byTrialcd$item)

bf_FULLcd = generalTestBF(RTtransf ~ primec + subject*primec+item,
                          data = byTrialcd, whichRandom = c("subject", "item") )%>% sort()
bf_FULLcd
#BF10
sort(bf_FULLcd)[8]/sort(bf_FULLcd)[9]


#BF10 d vs e
byTrialde <- nbyTrial %>% 
  filter(primec=="d" | primec=="e") %>% 
  filter(RT !="")

byTrialde$primec=as.factor(byTrialde$primec)
byTrialde$subject=as.factor(byTrialde$subject)
byTrialde$item=as.factor(byTrialde$item)

bf_FULLde = generalTestBF(RTtransf ~ primec + subject*primec+item,
                          data = byTrialde, whichRandom = c("subject", "item") )%>% sort()
bf_FULLde
#BF10
sort(bf_FULLde)[8]/sort(bf_FULLde)[9]










