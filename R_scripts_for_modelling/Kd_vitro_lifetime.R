## This script was used to determine the in vitro Kd for Tq-Ca-FLITS from lifetime data
# also the lifetime data of jGCaMP7c is visualized

## Part 1:contains
# calculation of the G and S coordinates of both sensors from the lifetime
# determines the fraction of Tq-Ca-FLITS in the calcium bound state for each concentration
# fitting of model to Tq-Ca-FLITS lifetimes and fraction to determine Kd

## Part 2: contains
# plotting the models with the measured data of Tq=Ca-FLITS
# plotting the measured data of both sensors on a polar plot

#for summarizing of data
require(dplyr)
#for plotting
require(ggplot2)

setwd("/Users/Franka/surfdrive/mTQ2 paper/to_Github")
FLITS <- read.csv("Kd_vitro_lifetime_Tq-Ca-FLITS.csv")
jGCaMP7c <- read.csv("Kd_vitro_lifetime_jGCaMP7c.csv")

###### PART 1: Calculate into G and S coordinates and fit Kd to Tq-Ca-FLITS data######
#####calculate lifetimes to Phi, M, G and S####
addPolarCoordinates <- function(data, MHz = 40) {
  tphi <- data$tauPhi
  tmod <- data$tauM
  omega = 2*pi*MHz/1000
  Phi = atan(omega*tphi)
  M=sqrt(1/ (1 + (omega*tmod)^2 ) )
  G = M * cos(Phi)
  S = M * sin(Phi)
  new <- data.frame(Phi=Phi, M=M, G=G, S=S)
  data <- cbind(data, new)
}

FLITS <- addPolarCoordinates(FLITS)
jGCaMP7c <- addPolarCoordinates(jGCaMP7c)

##summarize the data
sum_FLITS <- FLITS %>%
  group_by(free_Ca) %>%
  summarize(n=n(),
            mean_tauPhi=mean(tauPhi),
            mean_tauM=mean(tauM),
            mean_G=mean(G),
            mean_S=mean(S))

sum_jGCaMP7c <- jGCaMP7c %>%
  group_by(free_Ca) %>%
  summarize(n=n(),
            mean_tauPhi=mean(tauPhi),
            mean_tauM=mean(tauM),
            mean_G=mean(G),
            mean_S=mean(S))

#### Fit the tphi and tmod, to find the calcium dependency of the lifetime
#define a modeling function
fit_variable_to_L <- function(data, colL, colVariable, Ka.init, n.init) {
  L <- as.matrix(data[colL])
  Fb <- as.matrix(data[colVariable])
  max.init <- max(Fb)
  min.init <- min(Fb)
  #start values for fitting sensor, free min and max lifetime
  startvalues <- c(Ka=Ka.init, n=n.init, max=max.init, min=min.init)
  nls(Fb ~ ((max-min)/((Ka/L)^n +1) + min), start=startvalues,
      model=TRUE, trace=FALSE, algorithm = "port")
}

#make a smooth fit. import coef(model) as parameters
make_smooth_fit <- function(parameters){
  conc_range <- c(0, 10^seq(-2,1.7,0.1))
  fit <- data.frame(free_Ca = conc_range)
  Ka <- as.numeric(parameters["Ka"])
  n <- as.numeric(parameters["n"])
  max <- as.numeric(parameters["max"])
  min <- as.numeric(parameters["min"])
  fit$pred <- ((max-min) / ((Ka/fit$free_Ca)^n + 1)+min)
  return(fit)
}

# do the modelling, find variables and make a smooth fit dataframe
model_tauPhi <- fit_variable_to_L(FLITS, "free_Ca", "tauPhi", 0.23, 1.8)
para_tauPhi <- coef(model_tauPhi)
ci_tauPhi <- confint(model_tauPhi)
fit_tauPhi <- make_smooth_fit(para_tauPhi)

model_tauM <- fit_variable_to_L(FLITS, "free_Ca", "tauM", 0.23, 1.8)
para_tauM <- coef(model_tauM)
ci_tauM <- confint(model_tauM)
fit_tauM <- make_smooth_fit(para_tauM)

##Use the top and bottom tauPhi and tauM from measurements to calculate the fraction, via Polar coordinates.
#separate the extreme states (fully bound or unbound)
extremes <- filter(sum_FLITS, free_Ca==0 | free_Ca==39)

## calculate dS and dG, using the mean_G and mean_S from the low state
FLITS$dG <- FLITS$G - extremes[1,]$mean_G
FLITS$dS <- FLITS$S - extremes[1,]$mean_S
extremes$dG <- extremes$mean_G - extremes[1,]$mean_G
extremes$dS <- extremes$mean_S - extremes[1,]$mean_S

##calculate inproduct (projection on the line between low and high state).
FLITS$IN <- FLITS$dG * extremes[2,]$dG + FLITS$dS * extremes[2,]$dS
extremes$IN <- extremes$dG * extremes[2,]$dG + extremes$dS * extremes[2,]$dS

## calculate a line fraction (inproduct divided by total length of line {= dG max state ^ 2 + dS max state ^ 2} = IN max state)
FLITS$lfraction <- FLITS$IN/extremes[2,]$IN
extremes$lfraction <- extremes$IN/extremes[2,]$IN

##convert to true fraction
Ratio <- 3.51 ##from pKa fitting, pH7.
FLITS$Fraction <- FLITS$lfraction/(Ratio*(1-FLITS$lfraction)+FLITS$lfraction)
extremes$Fraction <- extremes$lfraction/(Ratio*(1-extremes$lfraction)+extremes$lfraction)

# add dG, dS, IN, lfraction to summary (sum_FLITS)
extra <- FLITS %>%
  group_by(free_Ca) %>%
  summarize(n=n(),
            mean_dG=mean(dG),
            mean_dS=mean(dS),
            mean_IN=mean(IN),
            mean_lfraction=mean(lfraction),
            mean_Fraction=mean(Fraction))
sum_FLITS <- merge(sum_FLITS, extra)
rm(extra)


#Fitting to line fraction
model_lfraction <- fit_variable_to_L(FLITS, "free_Ca", "lfraction", 0.3, 1)
para_lfr <- coef(model_lfraction)
ci_lfr <- confint(model_lfraction)
fit_lfraction <- make_smooth_fit(para_lfr)

#Fitting to real Fraction
model_Fraction <- fit_variable_to_L(FLITS, "free_Ca", "Fraction", 0.3, 1)
para_Fr <- coef(model_Fraction)
ci_Fr <- confint(model_Fraction)
fit_Fraction <- make_smooth_fit(para_Fr)

#collect all parameters
parameters <- list(Fraction_fit = para_Fr, ci=ci_Fr,
               lfraction_fit = para_lfr, ci_lfraction = ci_lfr,
               tphi_fit = para_tauPhi, ci_tphi = ci_tauPhi,
               tmod_fit = para_tauM, ci_tmod = ci_tauM)



###### PART 2: plotting of the models with the measurement data.
## show plot with the fitted models and the measurements
#model of modulation lifetime
ggplot(FLITS, aes(x=free_Ca, y=tauM)) + geom_point(shape=1) + 
  geom_path(data=fit_tauM, aes(x=free_Ca, y=pred)) + scale_x_log10() +
  ggtitle(expression(paste("Modulation Lifetime vs calcium ", italic("in vitro")))) +
  xlab("Free Ca2+ (uM)") + ylab("Modulation Lifetime (ns)") +
  coord_cartesian(xlim=c(0.001,30))+
  theme_light()

#model of phase lifetime
ggplot(FLITS, aes(x=free_Ca, y=tauPhi)) + geom_point(shape=1) + 
  geom_path(data=fit_tauPhi, aes(x=free_Ca, y=pred)) + scale_x_log10() +
  ggtitle(expression(paste("Phase Lifetime vs calcium ", italic("in vitro")))) +
  xlab("Free Ca2+ (uM)") + ylab("Phase Lifetime (ns)") +
  coord_cartesian(xlim=c(0.001,30))+
  theme_light()

#model of the fraction in the calcium bound state
ggplot(FLITS, aes(x=free_Ca, y=Fraction)) + geom_point(shape=1) + 
  geom_path(data=fit_Fraction, aes(x=free_Ca, y=pred)) + scale_x_log10() +
  ggtitle(expression(paste("Fraction vs. calcium ", italic("in vitro")))) + 
  xlab("Free Ca2+ (uM)") + ylab("Fraction") +
  coord_cartesian(xlim=c(0.001,30))+
  theme_light()


### plot the lifetime measurements in a polar plot
#define a function for a circle
make_half_circle <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
#call the function, with proper centre
polar <- make_half_circle(c(0.5,0))
#add 0,0 to the dataframe, order the data according to x
polar <- rbind(polar,c(0,0))
#plot an empty polar plot
empty_polar <- ggplot(polar,aes(x,y)) + geom_path() + 
  coord_fixed(ratio=1, xlim=c(0,1), ylim=c(0,0.5)) + xlab("G") +ylab("S") + theme_light()

#set calciumconcentration to factor for plotting
FLITS$free_Ca <- as.factor(FLITS$free_Ca)
sum_FLITS$free_Ca <- as.factor(sum_FLITS$free_Ca)
jGCaMP7c$free_Ca <- as.factor(jGCaMP7c$free_Ca)
sum_jGCaMP7c$free_Ca <- as.factor(sum_jGCaMP7c$free_Ca)

#plot triplo and mean for of Tq-Ca-FLITS
empty_polar + geom_point(data=FLITS, aes(x=G, y= S, color=free_Ca), alpha=.3, shape=16, size=1.5) +
  geom_point(data=sum_FLITS, aes(x=mean_G, y= mean_S, color=free_Ca), shape=3, size=2) +
  geom_line(data=extremes, aes(x=mean_G, y= mean_S), color="black") +
  geom_point(data=extremes, aes(x=mean_G, y= mean_S), color="black", size=0.3) +
  ggtitle("Tq-Ca-FLITS in vitro calibration")

#plot triplo and mean of jGCaMP-data
empty_polar + 
  geom_point(data=jGCaMP7c, aes(x=G, y= S, color=free_Ca), alpha=.3, shape=16, size=1.5) +
  #geom_point(data=sum_jGCaMP7c, aes(x=mean_G, y= mean_S, color=free_Ca), shape=3, size=2) +
  ggtitle("jGCaMP7c in vitro lifetimes for changing calcium")
