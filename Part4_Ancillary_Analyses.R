#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Part 4: Ancillary Analyses - Carcass persistence   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Curated by Adam F. Parlin, PhD; Bradley J. Cosentino, PhD; James P. Gibbs, PhD
#email(s): afparlin@esf.edu; compecophys@gmail.com
#          Cosentino@hws.edu;
#          jpgibbs@esf.edu


#~~~~~~~~~~~#
# Packages  #
#~~~~~~~~~~~#
library(survival) #survival_3.3-1

#~~~~~~~~~~~~~~~~~~~~#
# Working Directory  #
#~~~~~~~~~~~~~~~~~~~~#
wd<-setwd("Set/Work/Directory/Here/to/load")
load(paste0(wd, "/All_Data.Rdata")) 

# Use the 'surv' dataframe from All_Data.RData
head(surv)

surv$SurvObj<-with(surv,Surv(time = Time, event = State, type="right"))

#Right Censored
ScCa.fit<-survfit(SurvObj~1, data=surv)
ScCa_Morph.fit<-survfit(SurvObj~morph, data=surv)
summary(ScCa.fit)
print(ScCa.fit,print.rmean=TRUE)
print(ScCa_Morph.fit, print.rmean=T)

plot(ScCa.fit)
plot(ScCa_Morph.fit)

#morph differences
ScCa.cox<-coxph(SurvObj~morph, data=surv)
cox.zph(ScCa.cox)
summary(ScCa.cox)
exp(confint(ScCa.cox))

survdiff(SurvObj~morph, data=surv)

