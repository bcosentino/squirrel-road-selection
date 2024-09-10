#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Part 3: iNaturalist Analyses  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Curated by Adam F. Parlin, PhD; Bradley J. Cosentino, PhD; James P. Gibbs, PhD
#email(s): afparlin@esf.edu; compecophys@gmail.com
#          Cosentino@hws.edu;
#          jpgibbs@esf.edu

#~~~~~~~~~~~~~~~#
#   Packages    #
#~~~~~~~~~~~~~~~#

library(rcompanion) # rcompanion_2.4.15; cramerV test
library(ggplot2)    # ggplot2_3.5.0
library(gridExtra)  # gridExtra_2.3 

#~~~~~~~~~~~~~~~~~~~~#
# Working Directory  #
#~~~~~~~~~~~~~~~~~~~~#
wd<-setwd("Set/Work/Directory/Here/to/load")
load(paste0(wd, "/All_Data.Rdata")) 

#Use the 'squirrels' dataframe frm the All_Data.RData

#checking totals
table(squirrels$morph)
table(squirrels$surface)
table(squirrels$status)

#Subset on black or gray with agreement. 
ScCa_all<-subset(squirrels, morph =="black" | morph =="gray")
ScCa_all<-subset(ScCa_all, status != "no agreement")

str(ScCa_all)

#Subset for photographed on road only
ScCa_road<-subset(ScCa_all, surface == "road")

str(ScCa_road)

#~~~~~~~~~~~~~~~~~~~~~~~#
# Whole Analysis - All  # 
#~~~~~~~~~~~~~~~~~~~~~~~#

#Total of observations
sum(table(ScCa_all$morph))

#Chisq Test
ScCa_chisq<-chisq.test(table(ScCa_all$morph, ScCa_all$status))

#Print Results
ScCa_chisq

#Observed vs. expected
ScCa_chisq$expected
ScCa_chisq$observed

#Effect size
cramerV(table(ScCa_all$morph, ScCa_all$status)) 

#~~~~~~~~~~~~~~~~~#
# Subset Analysis - DOR #
#~~~~~~~~~~~~~~~~~#

#DOR comparison, grey vs. melanic on roads. 
road_chisq<-chisq.test(table(ScCa_road$morph, ScCa_road$status))
road_chisq
road_chisq$expected
road_chisq$observed

#Effect size
cramerV(table(ScCa_road$morph, ScCa_road$status)) 

#Confidence intervals for DOR

# Define the data - see output from ScCa_chisq$observed
data <- data.frame(
  morph = c("Melanic", "Gray"),
  alive = c(12531, 94444),
  dead = c(174, 1581)
)

# Calculate the proportion dead
data$proportion_dead <- data$dead / (data$alive + data$dead)

# Add confidence intervals
ci <- binom::binom.confint(x = data$dead, (data$alive + data$dead), 
                           conf.level = 0.95, method="agresti-coull")
data$lower_ci <- ci$lower
data$upper_ci <- ci$upper

# Create the first plot
p1 <- ggplot(data, aes(x = morph, y = proportion_dead, color = morph)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0) +
  scale_color_manual(values = c("Melanic" = "black", "Gray" = "gray")) +
  ylim(0, 0.02) +
  labs(x = "Morph", y = "Proportion Dead on Road") +
  ggtitle("All squirrels") +
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 13), 
        axis.title = element_text(size = 16), 
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold", color = "black", 
                                  margin = margin(t = 10, b = 10)))

# Define the data - see output from road_chisq$observed
data2 <- data.frame(
  morph = c("Melanic", "Gray"),
  alive = c(103, 463),
  dead = c(82, 825)
)

# Calculate the proportion dead
data2$proportion_dead <- data2$dead / (data2$alive + data2$dead)

# Add confidence intervals
ci2 <- binom::binom.confint(x = data2$dead, (data2$alive + data2$dead), 
                            conf.level = 0.95, method="agresti-coull")
data2$lower_ci <- ci2$lower
data2$upper_ci <- ci2$upper

# Create the second plot
p2 <- ggplot(data2, aes(x = morph, y = proportion_dead, color = morph)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0) +
  scale_color_manual(values = c("Melanic" = "black", "Gray" = "gray")) +
  ylim(0, 1) +
  labs(x = "Morph", y = "Proportion Dead on Road") +
  ggtitle("Squirrels on roads") +
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 13), 
        axis.title = element_text(size = 16), 
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold", color = "black", 
                                  margin = margin(t = 10, b = 10)))

# Arrange the plots in two columns
grid.arrange(p1, p2, ncol = 2)



