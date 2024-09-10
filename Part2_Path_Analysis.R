#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Part 2: Path Analysis         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Curated by Adam F. Parlin, PhD; Bradley J. Cosentino, PhD; James P. Gibbs, PhD
#email(s): afparlin@esf.edu; compecophys@gmail.com
#          Cosentino@hws.edu;
#          jpgibbs@esf.edu

#~~~~~~~~~~~~~~~#
#   Packages    #
#~~~~~~~~~~~~~~~#

library(dagitty) # dagitty_0.3-4
library(brglm2)  # brglm2_0.9.2
library(arm)     # arm_1.12-2; for binned residual plots

#~~~~~~~~~~~~~~~~~~~~#
# Working Directory  #
#~~~~~~~~~~~~~~~~~~~~#
wd<-setwd("Set/Work/Directory/Here/to/load")
load(paste0(wd, "/All_Data.Rdata")) 

#Use the 'data' dataframe from All_Data.RData
head(data)
str(data)

data$AADT.log <- log(data$AADT)
data$frag.log <- log(data$frag)
data$humanDensity.log <- log(data$humanDensity)

data.g <- data[-which(data$morph == "melanic"),] 
data.m <- data[-which(data$morph == "grey"),] 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# DAG Model Logistic Regressions #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#####GRAY MORPH DAG AND ADJUSTMENT SETS#####
dag.g <- dagitty("dag{Distance -> Human; Distance -> Frag;
                      Human -> Building; Human -> Frag; Human -> Gray; 
                      Human -> Speed; Human -> Traffic; Human -> Cross;
                      Building -> Forest; Building -> Cross;
                      Forest -> Frag; Forest -> Cross; Forest -> Gray;
                      Frag -> Split; Frag -> Gray;
                      Gray -> Mortality;
                      Speed -> Mortality
                      Traffic -> Mortality;
                      Cross -> Mortality;
                      Split -> Mortality;
                }")
#~~~~~~~~~~~~~~~~~~~#
# Adjustment sets   #
#~~~~~~~~~~~~~~~~~~~#

# Speed #

adjustmentSets(dag.g, exposure="Speed", outcome = "Mortality", effect = "direct") #human

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Gray mortality ~ speed + human ; Adjustment 1 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mg.s <- glm(y~new_speed + humanDensity.log, 
    data = data.g, 
    family=binomial(link="logit"))
binnedplot(fitted(mg.s), residuals(mg.s, type="response"))

mg.s.br <- update(mg.s, method = "brglm_fit")
summary(mg.s.br)$coefficients
binnedplot(fitted(mg.s.br), residuals(mg.s.br, type="response"))

#predicted probability - Adjustment 1
mean_log_humanDensity <- mean(data.g$humanDensity.log)
speed_values <- seq(min(data.g$new_speed), max(data.g$new_speed), length.out = 100)
gm1_df <- data.frame(new_speed = speed_values, humanDensity.log = mean_log_humanDensity)

gm1_df <- cbind.data.frame(gm1_df, predict(mg.s.br, newdata = gm1_df, type = "response", se.fit=TRUE))
gm1_df$lwr <- gm1_df$fit - 1.96*gm1_df$se.fit
gm1_df$upr <- gm1_df$fit + 1.96*gm1_df$se.fit

plot(gm1_df$new_speed, gm1_df$fit, type = "l",
     xlab = "Speed (mph)", ylab = "Predicted DOR Probability",
     main = "", ylim=c(0,1),
     lwd = 2)
shade(rbind(gm1_df$lwr, gm1_df$upr), gm1_df$new_speed, col=col.alpha("Sky blue", alpha=0.6))
lines(gm1_df$fit ~ gm1_df$new_speed, lwd=2)
points(data.g$y ~ data.g$new_speed)

# Traffic #

adjustmentSets(dag.g, exposure="Traffic", outcome = "Mortality", effect = "direct") #human

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Gray mortality ~ traffic + human ; Adjustment 2 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

mg.t <- glm(y~AADT.log + humanDensity.log, 
            data = data.g, 
            family=binomial(link="logit"))
binnedplot(fitted(mg.t), residuals(mg.t, type="response"))

mg.t.br <- update(mg.t, method = "brglm_fit")
summary(mg.t.br)$coefficients
binnedplot(fitted(mg.t.br), residuals(mg.t.br, type="response"))

#Predicted probability - Adjustment 2 - with quartiles
mean_log_humanDensity <- mean(data.g$humanDensity.log)
AADT_values <- seq(min(data.g$AADT.log), max(data.g$AADT.log), length.out = 100)
gm2_df <- data.frame(AADT.log = AADT_values, humanDensity.log = mean_log_humanDensity)

gm2_df <- cbind.data.frame(gm2_df, predict(mg.t.br, newdata = gm2_df, type = "response", se.fit=TRUE))
gm2_df$lwr <- gm2_df$fit - 1.96*gm2_df$se.fit
gm2_df$upr <- gm2_df$fit + 1.96*gm2_df$se.fit

plot(gm2_df$AADT.log, gm2_df$fit, type = "l",
     xlab = "log(Traffic)", ylab = "Predicted DOR Probability",
     main = "", ylim=c(0,1),
     lwd = 2)
shade(rbind(gm2_df$lwr, gm2_df$upr), gm2_df$AADT.log, col=col.alpha("Sky Blue", alpha=0.6))
lines(gm2_df$fit ~ gm2_df$AADT.log, lwd=2)
lines(gm2_df$lwr ~ gm2_df$AADT.log, lwd=3)
lines(gm2_df$upr ~ gm2_df$AADT.log, lwd=3)
points(data.g$y ~ data.g$AADT.log)

# Adjustment 3 #
adjustmentSets(dag.g, exposure="Cross", outcome = "Mortality", effect = "direct") #forest, human

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Gray mortality ~ crossings + human + forest ; Adjustment 3  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mg.c <- glm(y~crossing + humanDensity.log + forest_cover, 
            data = data.g, 
            family=binomial(link="logit"))
binnedplot(fitted(mg.c), residuals(mg.c, type="response"))

mg.c.br <- update(mg.c, method = "brglm_fit")
summary(mg.c.br)$coefficients
binnedplot(fitted(mg.c.br), residuals(mg.c.br, type="response"))

expo(mg.c, type = "correction*")

#Predicted probability - Adjustment 3 
mean_forest_cover <- mean(data.g$forest_cover)
mean_log_humanDensity <- mean(data.g$humanDensity.log)
crossings <- seq(min(data.g$crossing), max(data.g$crossing), length.out = 100)

gm3_df <- data.frame(crossing = crossings, forest_cover = mean_forest_cover, humanDensity.log = mean_log_humanDensity)

gm3_df <- cbind.data.frame(gm3_df, predict(mg.c.br, newdata = gm3_df, type = "response", se.fit=TRUE))
gm3_df$lwr <- gm3_df$fit - 1.96*gm3_df$se.fit
gm3_df$upr <- gm3_df$fit + 1.96*gm3_df$se.fit


#gm3_df$predicted_prob <- predict(mg.c.br,newdata = gm3_df, type = "response")

plot(gm3_df$crossing , gm3_df$fit, type = "l",
     xlab = "Crossing (0=absent;1=present)", ylab = "Predicted Probability",
     main = "Predicted Probability vs. Crossing (adjustment 3)",ylim=c(0,0.5),
     col = "blue", lwd = 2)

# Adjustment 4 #
adjustmentSets(dag.g, exposure="Split", outcome = "Mortality", effect = "direct") #frag

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Gray mortality ~ split + fragmentation  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mg.h <- glm(y~habsplit + frag.log, 
            data = data.g, 
            family=binomial(link="logit"))
binnedplot(fitted(mg.h), residuals(mg.h, type="response"))

mg.h.br <- update(mg.h, method = "brglm_fit")
summary(mg.h.br)$coefficients
binnedplot(fitted(mg.h.br), residuals(mg.h.br, type="response"))

expo(mg.h.br, type = "correction*")

#Predicted probability - Adjustment 4
mean_log_frag <- mean(data.g$frag.log)

habsplit<-seq(min(data.g$habsplit), max(data.g$habsplit), length.out = 100)

gm4_df <- data.frame(habsplit = habsplit, frag.log = mean_log_frag)

gm4_df <- cbind(gm4_df, predict(mg.h.br, newdata=gm4_df, type = "response", se.fit=TRUE))
gm4_df$lwr <- gm4_df$fit - 1.96*gm4_df$se.fit
gm4_df$upr <- gm4_df$fit + 1.96*gm4_df$se.fit

plot(gm4_df$habsplit , gm4_df$fit, type = "l",
     xlab = "Habitat Split", ylab = "Predicted Probability",
     main = "Predicted Probability vs. Habitat Split (adjustment 4)",ylim=c(0,0.5),
     col = "blue", lwd = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Melanic logistic regression #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dag.m <- dagitty("dag {Distance -> Human; Distance -> Frag; Distance -> Melanic
                      Human -> Building; Human -> Frag; 
                      Human -> Speed; Human -> Traffic; Human -> Cross;
                      Building -> Forest; Building -> Cross;
                      Forest -> Frag; Forest -> Cross; Forest -> Melanic;
                      Frag -> Split; Frag -> Melanic;
                      Melanic -> Mortality;
                      Speed -> Mortality
                      Traffic -> Mortality;
                      Cross -> Mortality;
                      Split -> Mortality;
                }")

#~~~~~~~~~~~~~~~~~~~#
# Adjustment sets   #
#~~~~~~~~~~~~~~~~~~~#

# Adjustment 1 #
adjustmentSets(dag.m, exposure="Speed", outcome = "Mortality", effect = "direct") #human

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Melanic mortality ~ speed + human ; Adjustment 1 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mm.s <- glm(y~new_speed + humanDensity.log, 
            data = data.m, 
            family=binomial(link="logit"))
binnedplot(fitted(mm.s), residuals(mm.s, type="response"))

mm.s.br <- update(mm.s, method = "brglm_fit")
summary(mm.s.br)$coefficients
binnedplot(fitted(mm.s.br), residuals(mm.s.br, type="response"))

expo(mm.s.br, type = "correction*") #doesn't seem to work on brglm_fit???

#Predicted probability - Adjustment 1
mean_log_humanDensity.m <- mean(data.m$humanDensity.log)

speed_values.m <- seq(min(data.m$new_speed), max(data.m$new_speed), length.out = 100)

mm1_df <- data.frame(new_speed = speed_values.m, humanDensity.log = mean_log_humanDensity.m)

mm1_df <- cbind(mm1_df, predict(mm.s.br, newdata = mm1_df, type = "response", se.fit=TRUE))

mm1_df$lwr <- mm1_df$fit - 1.96*mm1_df$se.fit
mm1_df$upr <- mm1_df$fit + 1.96*mm1_df$se.fit

plot(mm1_df$new_speed, mm1_df$fit, type = "l",
     xlab = "Speed (mph)", ylab = "Predicted Probability",
     main = "Predicted Probability vs. Speed (adjustment 1)", ylim=c(0,0.5),
     col = "blue", lwd = 2)

# Adjustment 2 #
adjustmentSets(dag.m, exposure="Traffic", outcome = "Mortality", effect = "direct") #human

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Melanic mortality ~ traffic + human ; Adjustment 2 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mm.t <- glm(y~AADT.log + humanDensity.log, 
            data = data.m, 
            family=binomial(link="logit"))
binnedplot(fitted(mm.t), residuals(mm.t, type="response"))

mm.t.br <- update(mm.t, method = "brglm_fit")
summary(mm.t.br)$coefficients
binnedplot(fitted(mm.t.br), residuals(mm.t.br, type="response"))

expo(mm.t.br, type = "correction*")

#Predicted Probability
mean_log_humanDensity.m <- mean(data.m$humanDensity.log)

AADT_log.m <- seq(min(log(data.m$AADT)), max(log(data.m$AADT)), length.out = 100)

mm2_df <- data.frame(AADT.log = AADT_log.m, humanDensity.log = mean_log_humanDensity.m)

mm2_df <- cbind(mm2_df, predict(mm.t.br, newdata = mm2_df, type = "response", se.fit=TRUE))
mm2_df$lwr <- mm2_df$fit - 1.96*mm2_df$se.fit
mm2_df$upr <- mm2_df$fit + 1.96*mm2_df$se.fit

plot(mm2_df$AADT.log , mm2_df$fit, type = "l",
     xlab = "Traffic (log transformed)", ylab = "Predicted Probability",
     main = "Predicted Probability vs. Traffic (adjustment 2)",ylim=c(0,0.5),
     col = "blue", lwd = 2)

# Adjustment 3 #
adjustmentSets(dag.m, exposure="Cross", outcome = "Mortality", effect = "direct") #forest, human

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Melanic mortality ~ crossings + human + forest ; Adjustment 3 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mm.c <- glm(y~crossing + humanDensity.log + forest_cover, 
            data = data.m, 
            family=binomial(link="logit"))
binnedplot(fitted(mm.c), residuals(mm.c, type="response"))

mm.c.br <- update(mm.c, method = "brglm_fit")
summary(mm.c.br)$coefficients
binnedplot(fitted(mm.c.br), residuals(mm.c.br, type="response"))

expo(mm.c.br, type = "correction*")

#Predicted probability - Adjustment 3
mean_forest_cover.m <- mean(data.m$forest_cover)
mean_log_humanDensity.m <- mean(data.m$humanDensity.log)

crossings.m <- seq(min(data.m$crossing), max(data.m$crossing), length.out = 100)

mm3_df <- data.frame(crossing = crossings.m, forest_cover = mean_forest_cover.m, humanDensity.log = mean_log_humanDensity.m)

mm3_df <- cbind(mm3_df, predict(mm.c.br, newdata = mm3_df, type = "response", se.fit=TRUE))
mm3_df$lwr <- mm3_df$fit - 1.96*mm3_df$se.fit
mm3_df$upr <- mm3_df$fit + 1.96*mm2_df$se.fit

plot(mm3_df$crossing , mm3_df$fit, type = "l",
     xlab = "Crossing (0=absent;1=present)", ylab = "Predicted Probability",
     main = "Predicted Probability vs. Crossing (adjustment 3)",ylim=c(0,0.5),
     col = "blue", lwd = 2)


# Adjustment 4 #
adjustmentSets(dag.m, exposure="Split", outcome = "Mortality", effect = "direct") #frag

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Melanic mortality ~ split + fragmentation ; Adjustment 4 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

mm.h <- glm(y~habsplit + frag.log, 
            data = data.m, 
            family=binomial(link="logit"))
binnedplot(fitted(mm.h), residuals(mm.h, type="response"))

mm.h.br <- update(mm.h, method = "brglm_fit")
summary(mm.h.br)$coefficients
binnedplot(fitted(mm.h.br), residuals(mm.h.br, type="response"))

expo(mm.h.br, type = "correction*")

#Predicted Probability
mean_log_frag.m <- mean(data.m$frag.log)

habsplit.m<-seq(min(data.m$habsplit), max(data.m$habsplit), length.out = 100)

mm4_df <- data.frame(habsplit = habsplit.m, frag.log = mean_log_frag.m)

mm4_df <-cbind(mm4_df, predict(mm.h.br, newdata = mm4_df, type = "response", se.fit=TRUE))
mm4_df$lwr <- mm4_df$fit - 1.96*mm4_df$se.fit
mm4_df$upr <- mm4_df$fit + 1.96*mm4_df$se.fit

plot(mm4_df$habsplit , mm4_df$fit, type = "l",
     xlab = "Habitat Split", ylab = "Predicted Probability",
     main = "Predicted Probability vs. Habitat Split (adjustment 4)",ylim=c(0,0.5),
     col = "blue", lwd = 2)
