#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Part 1: DOR v. Live Analyses    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Curated by Adam F. Parlin, PhD; Bradley J. Cosentino, PhD; James P. Gibbs, PhD
#email(s): afparlin@esf.edu; compecophys@gmail.com
#          Cosentino@hws.edu;
#          jpgibbs@esf.edu


#~~~~~~~~~~~~~~#
# Packages     #
#~~~~~~~~~~~~~~#

#libraries
library(sf)        #sf_1.0-9 
library(raster)    #raster_3.5-15  
library(ggplot2)   #ggplot2_3.5.0 
library(RANN)      #RANN_2.6.1 
library(jagsUI)    #jagsUI_1.5.2 
library(lubridate) #lubridate_1.9.0 
library(dplyr)     #dplyr_1.0.9

#~~~~~~~~~~~~~~~~~~~~#
# Working Directory  #
#~~~~~~~~~~~~~~~~~~~~#
wd<-setwd("Set/Work/Directory/Here/to/load")
load(paste0(wd, "/All_Data.Rdata")) 

#load routes
#Use 'long' GIS Shapefile from All_Data.RData
#Use 'short' GIS Shapefile from All_Data.RData

#load camera locations
# Use 'locs' dataframe from All_Data.RData
locs<- st_as_sf(locs,                                                                  
                crs = 4326,
                coords = c("long", "lat")) #WGS84 crs

#plot camera locations and routes
plot(locs$geometry)
plot(long$geometry, col="blue", add=T)
plot(short$geometry, col="red", add=T)

#buffer routes by 3 km, and add buffer to plot
sf_use_s2(TRUE) #set TRUE for st_buffer, FALSE when trying to extract vertices and polygons
route_buff<-st_buffer(long$geometry, dist=3000)
plot(route_buff, add=T)

#crop camera and point count sites to within 3-km buffer
locs.sub <- st_intersection(locs, route_buff)

#replot at scale of buffer
plot(route_buff)
plot(locs.sub$geometry, add=T)
plot(long$geometry, add=T, col="blue")
plot(short$geometry, col="red", add=T)

#save subset of sites
sites.sub <- as.data.frame(locs.sub)
save(sites.sub, file = paste0(wd.gis, "/sites_subset.Rdata"))

#####ESTIMATE CLINE FROM ROADKILL DATA#####

#load data
#Use 'road' dataframe from from All_Data.RData.

#quantify distance to city center for each mortality
city.ctr <- data.frame(id = "city.ctr", lon = -76.14373, lat = 43.04068)
city.ctr <- st_as_sf(city.ctr, 
                     crs = 4326,
                     coords = c("lon", "lat"))
pts <- data.frame(id = 1:nrow(road), lon = road$long, lat =road$lat, morph = road$morph )
pts <- st_as_sf(pts, 
                crs = 4326,
                coords = c("lon", "lat")) #WGS84 crs
road$dist.cent.km <- as.numeric(st_distance(pts, city.ctr)/1000) #distance from city center, km

#subset only dor observations
road.dor <- subset(road, morph != "reference")
road.dor$morph <- as.factor(road.dor$morph)
road.dor$dist.cent.km.s <- (road.dor$dist.cent.km - mean(road.dor$dist.cent.km))/sd(road.dor$dist.cent.km)

#bundle data
jags.data.live.dor <- list(y = ifelse(road.dor$morph == "melanic", 1, 0),
                           n = nrow(road.dor), 
                           distance = road.dor$dist.cent.km.s,
                           Xdist = seq(min(road.dor$dist.cent.km.s), max(road.dor$dist.cent.km.s), length.out = 100))
str(jags.data.live.dor)

#specify model in BUGS language
sink("clineEstimation_road.txt")
cat("
    model {
    
    #Priors 
    
    beta0.rd <- logit(mean.rd.pm)
    mean.rd.pm ~ dunif(0, 1)
    beta1.rd ~ dnorm(0, 0.1)     #pm slope for distance

    #Likelihood
    for(i in 1:n) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- beta0.rd + beta1.rd*distance[i]
    }
    
    #derived predictions
    for(k in 1:100) {
      cline.pred.rd[k] <- (exp(beta0.rd + beta1.rd*Xdist[k])) / (1 + (exp(beta0.rd + beta1.rd*Xdist[k])))
    }
    
    }
    ", fill=TRUE)
sink()

#Initial values
inits <- function() list(mean.pm.rd = runif(1), beta1.rd = rnorm(1))

#Parameters monitored
params <- c("beta0.rd", "beta1.rd", "cline.pred.rd", "mean.rd.pm")

##MCMC
ni <- 11000 ; nt <- 10 ;  nb <- 1000 ; nc <- 3

m.dor <- jags(jags.data.live.dor, inits, params, "clineEstimation_road.txt", 
              n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=20000,
              parallel=TRUE, n.cores=3)

##Model output----

#traceplots
traceplot(m.dor, parameters=c("beta0.rd", "beta1.rd"))

#parameter estimates with Rhat
library(MCMCvis)
MCMCsummary(m.dor, 
            params = c("beta0.rd", "beta1.rd", "mean.rd.pm"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

#prediction plots
X.dist.km <- jags.data.live.dor$Xdist*sd(road$dist.cent.km) + mean(road$dist.cent.km)

par(mfrow=c(1,1))
plot(m.dor$mean$cline.pred.rd ~ X.dist.km, type="l", lwd=2,
     main="Cline inferred from road mortalities", xlab="Distance from city center (km)", ylab="Proportion melanic",
     ylim=c(0,1), cex.lab=1.25)
polygon(x = c(X.dist.km, rev(X.dist.km)),
        y = c(m.dor$q97.5$cline.pred.rd, rev(m.dor$q2.5$cline.pred.rd)),
        col =  adjustcolor("red", alpha.f = 0.30), border = NA)
points(ifelse(road$morph == "melanic", 1, 0) ~ road$dist.cent.km)

#####ESTIMATE CLINE FOR LIVING SQUIRRELS#####

#Site info
#Use 'sites' dataframe from All_Data.RData 
site.coords <- sites[,c("site", "lat", "long")]
sites$site <- as.factor(sites$site)

#only retain sites within buffer of roadkill routes
sites <- sites[sites$site %in% sites.sub$site,]
sites$site <- as.factor(sites$site)

#Quantify distance to city center
pts <- data.frame(id = 1:nrow(sites), lon = sites$long, lat =sites$lat)
pts <- st_as_sf(pts, 
                crs = 4326,
                coords = c("lon", "lat")) #WGS84 crs
sites$dist.cent.km <- as.numeric(st_distance(pts, city.ctr))/1000 #distance from city center, km

#load occupancy data
#Use 'y.occ' dataframe from All_Data.RData

#subset to sites near roadkill routes
y.occ <- y.occ[which(y.occ$site %in% sites$site),]

#add rows for point count only sites and fill with NA
y.occ$site <- as.character(y.occ$site)
y.occ[nrow(y.occ)+sum(!(sites$site %in% y.occ$site)),] <- NA
y.occ[which(is.na(y.occ$site)), 1] <- as.character(sites[which(!(sites$site %in% y.occ$site)), "site"])
y.occ$site <- as.factor(y.occ$site)
y.occ <- y.occ[order(y.occ$site),] #order by site name

#site names with camera observations
occ.sites <- as.factor(y.occ[,1]) 

#create separate data frames for gray and melanic
y.occ.g <- cbind(y.occ$site, y.occ %>% dplyr::select(starts_with('g')))
y.occ.m <- cbind(y.occ$site, y.occ %>% dplyr::select(starts_with('m')))

#drop site names
y.occ.g <- y.occ.g[,-1]
y.occ.m <- y.occ.m[,-1]

#make sure y.occ is a numeric matrix
y.occ.g <- as.matrix(as.data.frame(lapply(y.occ.g, function(x) {as.numeric(levels(x))[x]})))
y.occ.m <- as.matrix(as.data.frame(lapply(y.occ.m, function(x) {as.numeric(levels(x))[x]})))

#detections for species
y.occ.both <- y.occ.g + y.occ.m
y.occ.both[y.occ.both>0] <- 1

#dates of deployment for camera observations
occ.dates <-  as.Date(occ.dates, "%m/%d/%y")
occ.dates.jul <- yday(occ.dates) #julian date
occ.start.date <- occ.dates[1] #need to extract start date to sink up with point count observations
occ.dates.seq <- occ.dates - occ.start.date #sequence of dates with 0 = first camera date

#count surveys per site
occ.surveys.per.site <- as.numeric(rowSums(!(is.na(y.occ.g[,2:ncol(y.occ.g)]))))  #distribution of number of surveys per site
summary(occ.surveys.per.site[which(!occ.surveys.per.site == 0)]) #summary excluding sites with no cams

#load point count data
#Use 'ct' dataframe from All_Data.RData

#subset to sites near roadkill routes
ct <- ct[which(ct$site %in% sites$site),]

#site names with point counts
ct.sites <- as.factor(ct$site)

#make sure site names for occ and ct match
identical(occ.sites, ct.sites)

#dates of point count observations
ct$survey1.date <- as.Date(ct$survey1.date,"%m/%d/%y")
ct$survey2.date <- as.Date(ct$survey2.date,"%m/%d/%y")
ct$survey3.date <- as.Date(ct$survey3.date,"%m/%d/%y")
ct$survey4.date <- as.Date(ct$survey4.date,"%m/%d/%y")
ct$survey5.date <- as.Date(ct$survey5.date,"%m/%d/%y")
ct$survey6.date <- as.Date(ct$survey6.date,"%m/%d/%y")
ct$survey7.date <- as.Date(ct$survey7.date,"%m/%d/%y")
ct.dates.seq <- cbind(ct$survey1.date-occ.start.date, ct$survey2.date-occ.start.date,   #sequence of dates with 0 = first camera date
                      ct$survey3.date-occ.start.date, ct$survey4.date-occ.start.date, 
                      ct$survey5.date-occ.start.date, ct$survey6.date-occ.start.date,
                      ct$survey7.date-occ.start.date)
ct.dates    <- cbind(as.character(ct$survey1.date), as.character(ct$survey2.date),   #sequence of dates with 0 = first camera date
                     as.character(ct$survey3.date), as.character(ct$survey4.date), 
                     as.character(ct$survey5.date), as.character(ct$survey6.date),
                     as.character(ct$survey7.date))
date.ct <- ct.dates.seq

ct.surveys.per.site <- rowSums(!(is.na(ct.dates.seq)))  #distribution of number of surveys per site

#Organize data for integrated model in JAGS

#organize survey data
y.ct.g <- cbind(ct$survey1.g, ct$survey2.g, ct$survey3.g, ct$survey4.g, ct$survey5.g, ct$survey6.g, ct$survey7.g) #observed gray counts from pt counts
y.ct.m <- cbind(ct$survey1.m, ct$survey2.m, ct$survey3.m, ct$survey4.m, ct$survey5.m, ct$survey6.m, ct$survey7.m) #observed melanic counts from pt counts

#survey design
nsites <- nrow(y.occ.m)       #number of sites with cam observations
nsurveys.occ <- ncol(y.occ.m) #number of surveys for occupancy data
nsurveys.ct <- ncol(y.ct.m)  #number of surveys for count data

#detection covariate: date
date.occ <- t(replicate(nsites, occ.dates.seq))  #create matrix of date by site, standardized
date.occ[is.na(y.occ.g)]  <- NA

#standardize date
date.mean <- mean(c(date.occ, ct.dates.seq), na.rm=T)
date.sd <- sd(c(date.occ, ct.dates.seq), na.rm=T)

date.occ.s <- (date.occ-date.mean)/date.sd
date.ct.s <- (ct.dates.seq-date.mean)/date.sd

#distance
distance <- sites$dist.cent.km
dist.city.s <- (distance - mean(distance))/sd(distance) #standardize distance

#impute missing dates with the mean
date.occ.s[is.na(date.occ.s)] <- 0
date.ct.s[is.na(date.ct.s)] <- 0
#see https://groups.google.com/g/hmecology/c/6xkFO5KwuhQ/m/cuX96R3AAQAJ on missing values for covariates when no surveys were done

sites$site <- factor(sites$site)
identical(sites$site, occ.sites) #make sure sites are ordered the same - site info vs. occ.sites
identical(sites$site, ct.sites) #make sure sites are ordered the same - site info vs. ct.sites

#weather data
#Use 'weather' dataframe from All_Data.RData
weather$DATE <- as.Date(weather$DATE, "%m/%d/%y")
weather$DATE <- as.character(weather$DATE)
temp.avg <- weather$TAVG #September 16, 2021 to October 18, 2022, inclusive of cam and pt ct survey dates
temp.avg.s <- (temp.avg - mean(temp.avg))/sd(temp.avg)
ct.dates <- matrix(ct.dates, dim(ct.dates))

temp.avg.occ <- t(replicate(nsites, temp.avg.s[1:393]))
temp.avg.ct <- ct.dates
temp.avg.ct[,1] <- temp.avg.s[match(ct.dates[,1], weather$DATE)]
temp.avg.ct[,2] <- temp.avg.s[match(ct.dates[,2], weather$DATE)]
temp.avg.ct[,3] <- temp.avg.s[match(ct.dates[,3], weather$DATE)]
temp.avg.ct[,4] <- temp.avg.s[match(ct.dates[,4], weather$DATE)]
temp.avg.ct[,5] <- temp.avg.s[match(ct.dates[,5], weather$DATE)]
temp.avg.ct[,6] <- temp.avg.s[match(ct.dates[,6], weather$DATE)]
temp.avg.ct[,7] <- temp.avg.s[match(ct.dates[,7], weather$DATE)]
temp.avg.ct <- matrix(as.numeric(temp.avg.ct), ncol = ncol(temp.avg.ct))

temp.avg.occ[is.na(y.occ.g)]  <- NA
temp.avg.occ[is.na(temp.avg.occ)] <- 0
temp.avg.ct[is.na(temp.avg.ct)] <- 0

#bundle data
jags.data.live.live <- list(yg.occ = y.occ.g, ym.occ =y.occ.m, 
                            yg.ct = y.ct.g, ym.ct =y.ct.m, 
                            nsites = nsites,
                            nsurveys.occ = nsurveys.occ, nsurveys.ct = nsurveys.ct, 
                            distance = dist.city.s,
                            tmp.occ = temp.avg.occ, tmp.ct = temp.avg.ct,
                            Xdist = seq(min(dist.city.s), max(dist.city.s), length.out = 100),
                            Xtemp = seq(min(c(temp.avg.occ, temp.avg.ct)), max(c(temp.avg.occ, temp.avg.ct)), length.out=100))
str(jags.data.live.live)

#specify model in JAGS language
sink("cline_integratedModel.txt")
cat("
    model {
    
    #Priors
    
    for(m in 1:2) {                   #m = morphs, 1 = melanic 2 = gray
    alpha0[m] <- logit(mean.p[m])     #detection intercept
    mean.p[m] ~ dunif(0, 0.2)           #mean detection probabilty at value = 0 for covariates (mean daily temperature)
    alpha1[m] ~ dnorm(0, 0.1)       #detection slope for temperature
    alpha2[m] ~ dnorm(0, 0.1)       #detection slope for temperature (quadratic term)
    }
    
    beta0.abu ~ dunif(1, 3.4)       #abundance intercept
    beta1.abu ~ dnorm(0, 0.1)       #abundance slope for distance
    mean.beta0.abu <- exp(beta0.abu)  #mean abundance at value 0 for covariates (mean distance)
    
    beta0.pm <- logit(mean.pm)        #intercept for proportion melanic
    mean.pm ~ dunif(0, 0.5)             #mean proportion melanic at value = 0 for covariates (mean distance)
    beta1.pm ~ dnorm(0, 0.1)        #proportion melanic slope for distance


    #Likelihood
    
    #Ecological process model for abundance
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])                                 #total squirrel abundance
    log(lambda[i]) <- beta0.abu + beta1.abu*distance[i]     #expected abundance as a function of distance
    Nm[i] ~ dbinom(pm[i], N[i])                             #coat color process / abundance of melanic morph
    logit(pm[i]) <- beta0.pm + beta1.pm*distance[i]         #expected proportion melanic as a function of distance
    Ng[i] <- N[i] - Nm[i]                                   #abundance of gray morph       
    }
    
    #Observation model for detection probability - occupancy data
    for(i in 1:nsites) {
    for(j in 1:nsurveys.occ) {
    ym.occ[i,j] ~ dbern(pstar.m.occ[i,j])                         #observed detections for melanic morph              
    pstar.m.occ[i,j] <- 1-(1-pdet.m.occ[i,j])^Nm[i]               #Pstar = P(detect melanic morph), pdet = ind. detection prob.
    logit(pdet.m.occ[i,j]) <- alpha0[1] +                         #ind. detection prob. as a function of survey temperature
                              alpha1[1]*tmp.occ[i,j] +
                              alpha2[1]*pow(tmp.occ[i,j], 2)
                                
    yg.occ[i,j] ~ dbern(pstar.g.occ[i,j])                         #observed detections for gray morph         
    pstar.g.occ[i,j] <- 1-(1-pdet.g.occ[i,j])^Ng[i]               #Pstar = P(detect gray morph), pdet = ind. detection prob.     
    logit(pdet.g.occ[i,j]) <- alpha0[2] +                         #ind. detection prob. as a function of survey temperature
                              alpha1[2]*tmp.occ[i,j] +
                              alpha2[2]*pow(tmp.occ[i,j], 2)
    }
    }
    
    #Observation model for detection probability - count data
    for(i in 1:nsites) {
    for(j in 1:nsurveys.ct) {
    ym.ct[i,j] ~ dbinom(pdet.m.ct[i,j], Nm[i])                   #observed melanic count; pdet = ind detection prob.
    logit(pdet.m.ct[i,j]) <- alpha0[1] +                         #ind. detection prob. as a function of survey temperature
                             alpha1[1]*tmp.ct[i,j] +
                             alpha2[1]*pow(tmp.ct[i,j], 2)
    
    yg.ct[i,j] ~ dbinom(pdet.g.ct[i,j], Ng[i])                   #observed gray count; pdet = ind detection prob.
    logit(pdet.g.ct[i,j]) <- alpha0[2] +                         #ind. detection prob. as a function of survey temperature
                             alpha1[2]*tmp.ct[i,j] +
                             alpha2[2]*pow(tmp.ct[i,j], 2)
    }
    }
    
    #Derived quantitites
    
    #predicted abundances and cline in melanism
    for(k in 1:100) {
    
      #total abundance as a function of distance
      lam.pred[k] <- exp(beta0.abu + beta1.abu*Xdist[k])
      
      #proportion melanic as a function of distance
      cline.pred[k] <- (exp(beta0.pm + beta1.pm*Xdist[k])) / (1 + (exp(beta0.pm + beta1.pm*Xdist[k])))
      
      #abundance of gray morph as a function of distance
      lam.g.pred[k] <- lam.pred[k] - (lam.pred[k]*cline.pred[k])
      
      #abudnance of melanic morph as a function of distance
      lam.m.pred[k] <- lam.pred[k]*cline.pred[k]
      
    }
    
    #predicted ind. detection probability as a function of survey temperature
    for(k in 1:100) {
    
      #melanic morph
      logit(pdet.m.pred[k]) <- alpha0[1] + alpha1[1]*Xtemp[k] + alpha2[1]*pow(Xtemp[k], 2)
      
      #gray morph
      logit(pdet.g.pred[k]) <- alpha0[2] + alpha1[2]*Xtemp[k] + alpha2[2]*pow(Xtemp[k], 2)
    }
    
    #proportion melanic at each sampled site
    for(i in 1:nsites) {
    pm.site[i] <- Nm[i]/ifelse(N[i]==0, 1e6, N[i])
    }
    
    }    
    ", fill=TRUE)
sink()

#Initial values - warnings for -Inf when no data; warn = F sets values to 0
y.g.max.ct <- apply(y.ct.g, 1, max, na.rm=TRUE, warn=F)
y.m.max.ct <- apply(y.ct.m, 1, max, na.rm=TRUE, warn=F)
y.g.max.occ <- apply(y.occ.g, 1, max, na.rm=TRUE, warn=F)
y.m.max.occ <- apply(y.occ.m, 1, max, na.rm=TRUE, warn=F)

inits <- function(){list(Nm = apply(cbind(y.m.max.occ,y.m.max.ct), 1, max), 
                         N = apply(cbind(y.m.max.occ,y.m.max.ct), 1, max) + apply(cbind(y.g.max.occ,y.g.max.ct), 1, max))}

#Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "mean.p",
            "beta0.abu", "beta1.abu", "mean.beta0.ab",
            "beta0.pm", "beta1.pm", "mean.pm",
            "lam.pred", "lam.g.pred", "lam.m.pred", "cline.pred",
            "pdet.m.pred", "pdet.g.pred",
            "N", "Nm", "Ng", "pm.site")


##MCMC
ni <- 11000 ; nt <- 10 ;  nb <- 1000 ; nc <- 3


m.live <- jags(jags.data.live.live, inits, params, "cline_integratedModel.txt", 
               n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=20000,
               parallel=TRUE, n.cores=3)


##Model output----

#traceplots
traceplot(m.live, parameters=c("alpha0", "alpha1", "alpha2",
                               "beta0.abu", "beta1.abu",
                               "beta0.pm", "beta1.pm"))

#parameter estimates with Rhat
MCMCsummary(m.live, 
            params = c("alpha0", "alpha1", "alpha2",
                       "beta0.abu", "beta1.abu",
                       "beta0.pm", "beta1.pm"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

#prediction plots
X.dist <- jags.data.live.live$Xdist
X.dist.km <- jags.data.live.live$Xdist*sd(distance) + mean(distance)


par(mfrow=c(2,2),
    mai=c(0.7, 0.7, 0.2, 0), #c(bottom, left, top, right) 
    omi=c(0.1,0.1,0.1, 0.1),
    mgp=c(2.5, 0.8, 0))

plot(m.live$mean$lam.pred ~ X.dist.km, type="l", lwd=2,
     main="", 
     xlab="", ylab="Total abundance",
     xlim=c(0, 11), ylim=c(0,40), cex.main = 1.8, cex.lab=1.3, cex.axis=1.2)
polygon(x = c(X.dist.km, rev(X.dist.km)),
        y = c(m.live$q97.5$lam.pred, rev(m.live$q2.5$lam.pred)),
        col =  adjustcolor("#4393C3", alpha.f = 0.50), border = NA)
lines(m.live$mean$lam.pred ~ X.dist.km, lwd=2)
points(m.live$mean$N ~ sites$dist.cent.km)
mtext(quote(bold("(") * bolditalic(a) * bold(") Total abundance")),
      side=3, line=0.2, at=-0.15, adj=0, cex=1.22)


plot(m.live$mean$lam.m.pred ~ X.dist.km, type="l", lwd=2,
     main="", xlab="", ylab="Melanic morph abundance",
     xlim=c(0, 11), ylim=c(0,30), cex.main = 1.8, cex.lab=1.3, cex.axis=1.2)
polygon(x = c(X.dist.km, rev(X.dist.km)),
        y = c(m.live$q97.5$lam.m.pred, rev(m.live$q2.5$lam.m.pred)),
        col =  adjustcolor("#4393C3", alpha.f = 0.50), border = NA)
lines(m.live$mean$lam.m.pred ~ X.dist.km, lwd=2)
points(m.live$mean$Nm ~ sites$dist.cent.km)
mtext(quote(bold("(") * bolditalic(b) * bold(") Melanic abundance")),
      side=3, line=0.2, at=-0.15, adj=0, cex=1.22)

plot(m.live$mean$lam.g.pred ~ X.dist.km, type="l", lwd=2,
     main="", xlab="Distance from city center (km)", ylab="Gray morph abundance",
     xlim=c(0, 11), ylim=c(0,30), cex.main = 1.8, cex.lab=1.3, cex.axis=1.2)
polygon(x = c(X.dist.km, rev(X.dist.km)),
        y = c(m.live$q97.5$lam.g.pred, rev(m.live$q2.5$lam.g.pred)),
        col =  adjustcolor("#4393C3", alpha.f = 0.50), border = NA)
lines(m.live$mean$lam.g.pred ~ X.dist.km, lwd=2)
points(m.live$mean$Ng ~ sites$dist.cent.km)
mtext(quote(bold("(") * bolditalic(c) * bold(") Gray abundance")),
      side=3, line=0.2, at=-0.15, adj=0, cex=1.22)


plot(m.live$mean$cline.pred ~ X.dist.km, type="l", lwd=2,
     cex.main = 1.8, cex.lab=1.3, cex.axis=1.2,
     main="", xlab="Distance from city center (km)", ylab="Proportion melanic",
     xlim=c(0, 11), ylim=c(0,1))
polygon(x = c(X.dist.km, rev(X.dist.km)),
        y = c(m.live$q97.5$cline.pred, rev(m.live$q2.5$cline.pred)),
        col =  adjustcolor("#4393C3", alpha.f = 0.40), border = NA)
lines(m.live$mean$cline.pred ~ X.dist.km, lwd=2)
points(m.live$mean$pm.site ~ sites$dist.cent.km, cex=0.8)
mtext(quote(bold("(") * bolditalic(d) * bold(") Proportion melanic")),
      side=3, line=0.2, at=-0.15, adj=0, cex=1.22)

#m.dor$mean$cline.pred.rd ~ X.dist.km

plot(m.live$mean$cline.pred ~ X.dist.km, type="l", lwd=2,
     cex.main = 1.8, cex.lab=1.3, cex.axis=1.2,
     main="", xlab="Distance from city center (km)", ylab="Proportion melanic",
     xlim=c(0, 11), ylim=c(0,1))
polygon(x = c(X.dist.km, rev(X.dist.km)),
        y = c(m.live$q97.5$cline.pred, rev(m.live$q2.5$cline.pred)),
        col =  adjustcolor("black", alpha.f = 0.40), border = NA)
lines(m.live$mean$cline.pred ~ X.dist.km, lwd=2)
lines(m.dor$mean$cline.pred.rd ~ X.dist.km, lwd=2, lty=2)
polygon(x = c(X.dist.km, rev(X.dist.km)),
        y = c(m.dor$q97.5$cline.pred.rd, rev(m.dor$q2.5$cline.pred.rd)),
        col =  adjustcolor("grey", alpha.f = 0.30), border = NA)
points(m.live$mean$pm.site ~ sites$dist.cent.km, cex=0.8)

#detection
X.temp <- jags.data.live.live$Xtemp*sd(temp.avg) + mean(temp.avg)

par(mfrow=c(1,1))
plot(m.live$mean$pdet.g.pred ~ X.temp, type="l", lwd=2,
     main="", xlab="Daily Mean Temperature (Celcius)", ylab="Individual detection probability",
     xlim=c(-20, 30), ylim=c(0,0.2))
polygon(x = c(X.temp, rev(X.temp)),
        y = c(m.live$q97.5$pdet.g.pred, rev(m.live$q2.5$pdet.g.pred)),
        col =  adjustcolor("skyblue", alpha.f = 0.30), border = NA)

lines(m.live$mean$pdet.m.pred ~ X.temp, lwd=2)
polygon(x = c(X.temp, rev(X.temp)),
        y = c(m.live$q97.5$pdet.m.pred, rev(m.live$q2.5$pdet.m.pred)),
        col =  adjustcolor("gray", alpha.f = 0.30), border = NA)

legend(-20, 0.2, c("Gray", "Melanic"), lwd=c(2,2), pch=c(15,15), col=c("skyblue", "gray"), bty="n")



#####Compare p(melanic) between living and DOR squirrels#####

library(boot)

pm.delta <- matrix(NA, nrow=100, ncol=3000)
for(i in 1:3000){
  pm.live <- inv.logit(m.live$sims.list$beta0.pm[i] + m.live$sims.list$beta1.pm[i]*X.dist)
  pm.dor <- inv.logit(m.dor$sims.list$beta0.rd[i] + m.dor$sims.list$beta1.rd[i]*X.dist)
  pm.delta[,i] <- pm.live - pm.dor
}

pm.delta.mean <- apply(pm.delta, 1, mean)
pm.delta.lcl <- apply(pm.delta, 1, function(x) quantile(x, prob=c(0.025)))
pm.delta.ucl <- apply(pm.delta, 1, function(x) quantile(x, prob=c(0.975)))

plot(pm.delta.mean ~ X.dist.km, type="l", lwd=2,
     main="", 
     xlab="Distance from city center (km)", ylab="Difference in proportion melanic (Live - DOR)",
     xlim=c(0, 11), ylim=c(-0.5,0.5), cex.main = 1.2, cex.lab=1.2, cex.axis=1.1)
polygon(x = c(X.dist.km, rev(X.dist.km)),
        y = c(pm.delta.ucl, rev(pm.delta.lcl)),
        col =  adjustcolor("#4393C3", alpha.f = 0.50), border = NA)
abline(0, 0, lty=2)

