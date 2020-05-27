## This script was used to determine the in vivo Kd for Tq-Ca-FLITS from lifetime data

## Part 1:contains
# calculation of the G and S coordinates from the lifetime
# determines the fraction of Tq-Ca-FLITS in the calcium bound state for each concentration
# fitting of model to Tq-Ca-FLITS lifetimes and fraction to determine Kd

## Part 2: contains
# plotting the models with the measured data of Tq-Ca-FLITS
# plotting the measured data of both sensors on a polar plot


# cleaning and tidying data
library(dplyr)
#plotting
library(ggplot2)

link <- "https://raw.githubusercontent.com/Franka-van-der-Linden/Quantitative-Calcium-Imaging/master/R_scripts_for_modelling/Kd_vivo_lifetime_Tq-Ca-FLITS.csv"
FLITS_vivo <- read.csv(link)

###### PART 1: calculations and fitting
##calculate lifetimes to Phi, M, G and S####
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

FLITS_vivo <- addPolarCoordinates(FLITS_vivo)

#summarize the data of my sensor, before continuing to fitting to tphi and tmod
sum_vivo <- FLITS_vivo %>%
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

#make a smooth fit. Import coef(model) as parameters
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

model_tauPhi <- fit_variable_to_L(FLITS_vivo, "free_Ca", "tauPhi", 0.3, 1)
para_tauPhi <- coef(model_tauPhi)
ci_tauPhi <- confint(model_tauPhi)
fit_tauPhi <- make_smooth_fit(para_tauPhi)

model_tauM <- fit_variable_to_L(FLITS_vivo, "free_Ca", "tauM", 0.1, 1)
para_tauM <- coef(model_tauM)
ci_tauM <- confint(model_tauM)
fit_tauM <- make_smooth_fit(para_tauM)


##Use the min and max tauPhi and tauM from the fitting to calculate the fraction, via Polar coordinates.
#separate the extreme states (fully bound or unbound)
extremes <- data.frame(tauPhi=para_tauPhi,tauM=para_tauM)[c(3,4),]
extremes$free_Ca <- c(39,0)

#calculated G and S coordinate of these lifetimes
extremes <- addPolarCoordinates(extremes)

## calculate dS and dG in FLITS_vivo, using the mean_G and mean_S from the low state
FLITS_vivo$dG <- FLITS_vivo$G - extremes["min",]$G
FLITS_vivo$dS <- FLITS_vivo$S - extremes["min",]$S
extremes$dG <- extremes$G - extremes["min",]$G
extremes$dS <- extremes$S - extremes["min",]$S

##calculate inproduct (projection on the line between low and high state).
FLITS_vivo$IN <- FLITS_vivo$dG * extremes["max",]$dG + FLITS_vivo$dS * extremes["max",]$dS
extremes$IN <- extremes$dG * extremes["max",]$dG + extremes$dS * extremes["max",]$dS

## calculate a line fraction (inproduct divided by total length of line {= dG max state ^ 2 + dS max state ^ 2} = IN max state)
FLITS_vivo$lfraction <- FLITS_vivo$IN/extremes[1,]$IN
extremes$lfraction <- extremes$IN/extremes[1,]$IN

##convert to true fraction
Ratio <- 3.02 ##from HeLa's stimulated with ionomycin: Fmax/F0, bg corrected
FLITS_vivo$Fraction <- FLITS_vivo$lfraction/(Ratio*(1-FLITS_vivo$lfraction)+FLITS_vivo$lfraction)
extremes$Fraction <- extremes$lfraction/(Ratio*(1-extremes$lfraction)+extremes$lfraction)

# add dG, dS, IN, lfraction, Fraction  sum_vivo
extra <- FLITS_vivo %>%
  group_by(free_Ca) %>%
  summarize(n=n(),
            mean_dG=mean(dG),
            mean_dS=mean(dS),
            mean_IN=mean(IN),
            mean_lfraction=mean(lfraction),
            mean_Fraction=mean(Fraction))
sum_vivo <- merge(sum_vivo, extra)
rm(extra)


##Now I want to fit a model to the line fraction or the real Fraction
## Fix the min and max: I already fitted these with the lifetimes free.
## So no more changing the min and max lifetimes.

fit_variable_to_L2 <- function(data, colL, colVariable, Ka.init, n.init) {
  L <- as.matrix(data[colL])
  Fb <- as.matrix(data[colVariable])
  max = 1
  min = 0
  #start values for fitting sensor, fixed min and max lifetime
  startvalues <- c(Ka=Ka.init, n=n.init)
  nls(Fb ~ ((max-min)/((Ka/L)^n +1) + min), start=startvalues,
      model=TRUE, trace=FALSE, algorithm = "port")
}

#make a smooth fit. Import coef(model) as parameters
make_smooth_fit2 <- function(parameters){
  conc_range <- c(0, 10^seq(-2,1.7,0.1))
  fit <- data.frame(free_Ca = conc_range)
  Ka <- as.numeric(parameters["Ka"])
  n <- as.numeric(parameters["n"])
  max <- 1
  min <- 0
  fit$pred <- ((max-min) / ((Ka/fit$free_Ca)^n + 1)+min)
  return(fit)
}

#Fitting to line fraction. Use new function with fixed min and max
model_lfraction <- fit_variable_to_L2(FLITS_vivo, "free_Ca", "lfraction", 0.1, 1)
para_lfraction <- coef(model_lfraction)
ci_lfraction <- confint(model_lfraction)
fit_lfraction <- make_smooth_fit2(para_lfraction)

#Fitting to real Fraction. Use new function with fixed min and max
model_Fraction <- fit_variable_to_L2(FLITS_vivo, "free_Ca", "Fraction", 0.1, 1)
para_Fraction <- coef(model_Fraction)
ci_Fraction <- confint(model_Fraction)
fit_Fraction <- make_smooth_fit2(para_Fraction)


##store all fitting parameters. Don't forget to add the min and max values with NA as ci.
para_lfraction <- c(para_lfraction, "max"=1, "min"=0)
ci_lfraction <- rbind(ci_lfraction,ci_lfraction)
row.names(ci_lfraction) <- c("Ka","n","max", "min")
ci_lfraction[c("max","min"),] <- NA

para_Fraction <- c(para_Fraction, "max"=1, "min"=0)
ci_Fraction <- rbind(ci_Fraction,ci_Fraction)
row.names(ci_Fraction) <- c("Ka","n","max", "min")
ci_Fraction[c("max","min"),] <- NA

params <- list(lfraction_fit = para_lfraction, ci_lfraction = ci_lfraction,
               Fraction_fit = para_Fraction, ci_Fraction = ci_Fraction,
               tphi_fit = para_tauPhi, ci_tphi = ci_tauPhi,
               tmod_fit = para_tauM, ci_tmod = ci_tauM)



###### PART 2: plotting data
#### plot on normal plotarea, all models and the measured data

# Fraction in high state
ggplot(FLITS_vivo, aes(x=free_Ca, y=Fraction)) + geom_point(shape=16, size=1.8) + 
  geom_path(data=fit_Fraction, aes(x=free_Ca, y=pred)) + scale_x_log10() +
  ggtitle(expression(paste("Fraction vs calcium ", italic("in vivo")))) + 
  xlab("Free Ca2+ (uM)") + ylab("Fraction") + 
  coord_cartesian(xlim=c(0.001,30))+
  theme_light()

# line fraction
ggplot(FLITS_vivo, aes(x=free_Ca, y=lfraction)) + geom_point(shape=16, size=1.8) + 
  geom_path(data=fit_lfraction, aes(x=free_Ca, y=pred)) + scale_x_log10() +
  ggtitle(expression(paste("Line fraction vs calcium ", italic("in vivo")))) + 
  xlab("Free Ca2+ (uM)") + ylab("Line fraction") + 
  coord_cartesian(xlim=c(0.001,30))+
  theme_light()

# phase lifetime
ggplot(FLITS_vivo, aes(x=free_Ca, y=tauPhi)) + geom_point(shape=16, size=1.8) + 
  geom_path(data=fit_tauPhi, aes(x=free_Ca, y=pred)) + scale_x_log10() +
  ggtitle(expression(paste("Phase lifetime vs calcium ", italic("in vivo")))) + 
  xlab("Free Ca2+ (uM)") + ylab("tauPhi (ns)") + 
  coord_cartesian(xlim=c(0.001,30))+
  theme_light()

# modulation lifetime
ggplot(FLITS_vivo, aes(x=free_Ca, y=tauM)) + geom_point(shape=16, size=1.8) + 
  geom_path(data=fit_tauM, aes(x=free_Ca, y=pred)) + scale_x_log10() +
  ggtitle(expression(paste("Modulation lifetime vs calcium ", italic("in vivo")))) + 
  xlab("Free Ca2+ (uM)") + ylab("tauM (ns)") + 
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

#### now add the real data
#set calciumconcentration to factor for plotting
FLITS_vivo$free_Ca <- as.factor(round(FLITS_vivo$free_Ca, digits=3))
sum_vivo$free_Ca <- as.factor(round(sum_vivo$free_Ca, digits=3))

#plot single measurements and their mean. 
#add line from the min to max as calculated from the tphi and tmod fit, this is already stored in the extremes parameter
empty_polar + geom_point(data=FLITS_vivo, aes(x=G, y= S, color=free_Ca), alpha=.3, shape=16, size=1.5) +
  geom_point(data=sum_vivo, aes(x=mean_G, y= mean_S, color=free_Ca), shape=3, size=2) +
  geom_line(data=extremes, aes(x=G, y= S), color="black") +
  geom_point(data=extremes, aes(x=G, y= S), color="black", size=0.3) +
  ylab("S") + xlab("G") +
  scale_x_continuous(breaks=seq(0.6,1,0.2), minor_breaks = seq(0.6,1,0.1)) +
  scale_y_continuous(breaks=seq(0.2,0.5,0.2), minor_breaks = seq(0.2,0.5,0.1)) +
  theme_light() + 
  coord_fixed(ratio=1, xlim=c(0.6,1), ylim=c(0.2,0.5))

