
################Covid-19 Copy Number (Random) Using specific Mean and SD###############
#                             Mohammad Nayeem Hasan                                   #
#######################################################################################

library(MASS)
require(foreign)
require(ggplot2)
require(maptools)
library(tidyverse)
library(betareg)
library(car)
library(gapminder)
library(dplyr)

setwd("E:\\ResearchProject\\Aminul\\Copy Number")


VLData <- read.csv("MainData.csv")
VLData$Viral.Load.Cy5.

library(nortest)
ad.test(VLData$Viral.Load.Cy5.)


library(fitdistrplus)
library(fitdistrplus)
library(logspline)

my_data <- c(702442.0,  1702933.2,   643533.4,   539461.6,   889656.6,  1714387.8,  1714387.8,
             1576932.6, 12000000.0, 276336.5,  1817479.2,   321499.7,   392190.1,  2962939.2,
             378444.7,   614079.1,   981293.4,  2665119.6, 1439477.4, 12000000.0, 10770000.0, 
             12000000.0, 12000000.0,  456989.5, 12000000.0,  1909116.0,   174228.2, 1221840.0,
             551243.3,   563025.0,   354881.3, 12000000.0, 12000000.0,   619970.0,  1072930.2,
             476625.7, 347026.8, 10500000.0,  2688028.8,   386299.2,   720114.6,  3363850.2, 11220000.0,
             11490000.0,   372553.9, 606224.6,   547316.0,   660564.6,   213500.6,  1210385.4,   
             445207.8,   380408.3,   405935.4,   529643.5, 712260.1,   515898.1,   555170.5,   866747.4,
             345063.2,   207609.8,   667096.9,   798019.8,   889656.6, 12000000.0,  2848393.2,   
             586588.4,  2195481.0,   559097.8,   244918.6, 12000000.0,   368626.6,   504116.4,
             3146212.8,  1702933.2,   584624.8,   889656.6,  1714387.8)

descdist(my_data)

fit.weibull <- fitdist(my_data, "weibull", method = "mle")
summary(fit.weibull)
plot(fit.weibull)
fit.weibull$aic

fit.lognorm <- fitdist(my_data, "lnorm")
summary(fit.lognorm)
plot(fit.lognorm)
fit.lognorm$aic

fnbinom <-fitdist(my_data,"nbinom")
summary(fnbinom)
plot(fnbinom)
fnbinom$aic

fnom<-fitdist(my_data,"norm")
summary(fnom)
plot(fnom)
fnom$aic

RF <- rlnorm(n = 2400, meanlog = log(10), sdlog = log(2))
RF

options(scipen = 999) 

F <- round(rlnorm(n = 1200, meanlog = log(149), sdlog = log(95)),2)
F

E <- rlnorm(n = 80)*100
E

################Descriptive#############
library(ggpubr)
library(pastecs)
options(scipen = 999)

Q <- c(4.00,
       3.81,
       4.71,
       4.00,
       3.20,
       3.81,
       3.33,
       4.71,
       2.00,
       3.33,
       2.96,
       4.21,
       4.71,
       6.67,
       3.20
)
stat.desc(Q)
rnorm(120,3.91, 1.07)


abs(rnorm(120,149, 95))


round(runif(100, min=100, max=200),2)


options(scipen = 999) 
Data <- read.csv("AnalysisData.csv")
attach(Data)
Data$ViralLoad_N <- (Data$ViralLoad - mean(Data$ViralLoad)) / sd(Data$ViralLoad)
Data$FlowRate_N <- (Data$FlowRate - mean(Data$FlowRate)) / sd(Data$FlowRate)
Data$VLinStool_N <- (Data$VLinStool - mean(Data$VLinStool)) / sd(Data$VLinStool)
Data$StoolProduc_N <- (Data$StoolProduc - mean(Data$StoolProduc)) / sd(Data$StoolProduc)
Data$Shedvirus_N <- (Data$Shedvirus - mean(Data$Shedvirus)) / sd(Data$Shedvirus)
Data$ppl_N <- (Data$ppl - mean(Data$ppl)) / sd(Data$ppl)

fit <- lm(Data$ppl  ~ Data$ViralLoad + Data$FlowRate + Data$VLinStool + Data$StoolProduc + Data$Shedvirus, data = Data)
summary(fit)
fit$coefficients
round(confint(fit),2)
mean(abs(predict(fit)))

mean(Data$ViralLoad)


VLData <- read.csv("Short Analysis File_Red.csv")

mean(VLData$Ct.value..IC., na.rm=T)
sd(VLData$Ct.value..IC., na.rm=T)

mean(VLData$Ct.value..IC..1, na.rm=T)
sd(VLData$Ct.value..IC..1, na.rm=T)
t.test(VLData$Ct.value..IC., VLData$Ct.value..IC..1)



mean(VLData$Ct.Value..ORF1ab., na.rm=T)
sd(VLData$Ct.Value..ORF1ab., na.rm=T)

mean(VLData$Ct.Value..ORF1ab..1, na.rm=T)
sd(VLData$Ct.Value..ORF1ab..1, na.rm=T)
t.test(VLData$Ct.Value..ORF1ab., VLData$Ct.Value..ORF1ab..1, na.rm=T)


mean(VLData$Ct.Value..N., na.rm=T)
sd(VLData$Ct.Value..N., na.rm=T)

mean(VLData$Ct.Value..N..1, na.rm=T)
sd(VLData$Ct.Value..N..1, na.rm=T)
t.test(VLData$Ct.Value..N., VLData$Ct.Value..N..1, na.rm=T)


########Yellow
VLData <- read.csv("Short Analysis File_Yellow.csv")

mean(VLData$Ct.value..IC., na.rm=T)
sd(VLData$Ct.value..IC., na.rm=T)

mean(VLData$Ct.value..IC..1, na.rm=T)
sd(VLData$Ct.value..IC..1, na.rm=T)
t.test(VLData$Ct.value..IC., VLData$Ct.value..IC..1)



mean(VLData$Ct.Value..ORF1ab., na.rm=T)
sd(VLData$Ct.Value..ORF1ab., na.rm=T)

mean(VLData$Ct.Value..ORF1ab..1, na.rm=T)
sd(VLData$Ct.Value..ORF1ab..1, na.rm=T)
t.test(VLData$Ct.Value..ORF1ab., VLData$Ct.Value..ORF1ab..1, na.rm=T)


mean(VLData$Ct.Value..N., na.rm=T)
sd(VLData$Ct.Value..N., na.rm=T)

mean(VLData$Ct.Value..N..1, na.rm=T)
sd(VLData$Ct.Value..N..1, na.rm=T)
t.test(VLData$Ct.Value..N., VLData$Ct.Value..N..1, na.rm=T)



VLData <- read.csv("Short Analysis File.csv")

mean(VLData$Ct.value..IC., na.rm=T)
sd(VLData$Ct.value..IC., na.rm=T)

mean(VLData$Ct.value..IC..1, na.rm=T)
sd(VLData$Ct.value..IC..1, na.rm=T)
t.test(VLData$Ct.value..IC., VLData$Ct.value..IC..1)



mean(VLData$Ct.Value..ORF1ab., na.rm=T)
sd(VLData$Ct.Value..ORF1ab., na.rm=T)

mean(VLData$Ct.Value..ORF1ab..1, na.rm=T)
sd(VLData$Ct.Value..ORF1ab..1, na.rm=T)
t.test(VLData$Ct.Value..ORF1ab., VLData$Ct.Value..ORF1ab..1, na.rm=T)


mean(VLData$Ct.Value..N., na.rm=T)
sd(VLData$Ct.Value..N., na.rm=T)

mean(VLData$Ct.Value..N..1, na.rm=T)
sd(VLData$Ct.Value..N..1, na.rm=T)
t.test(VLData$Ct.Value..N., VLData$Ct.Value..N..1, na.rm=T)





VLData <- read.csv("Short Analysis File_Red_myme.csv")

mean(VLData$Ct.value..IC., na.rm=T)
sd(VLData$Ct.value..IC., na.rm=T)

mean(VLData$Ct.value..IC..1, na.rm=T)
sd(VLData$Ct.value..IC..1, na.rm=T)
t.test(VLData$Ct.value..IC., VLData$Ct.value..IC..1)



mean(VLData$Ct.Value..ORF1ab., na.rm=T)
sd(VLData$Ct.Value..ORF1ab., na.rm=T)

mean(VLData$Ct.Value..ORF1ab..1, na.rm=T)
sd(VLData$Ct.Value..ORF1ab..1, na.rm=T)
t.test(VLData$Ct.Value..ORF1ab., VLData$Ct.Value..ORF1ab..1, na.rm=T)


mean(VLData$Ct.Value..N., na.rm=T)
sd(VLData$Ct.Value..N., na.rm=T)

mean(VLData$Ct.Value..N..1, na.rm=T)
sd(VLData$Ct.Value..N..1, na.rm=T)
t.test(VLData$Ct.Value..N., VLData$Ct.Value..N..1, na.rm=T)