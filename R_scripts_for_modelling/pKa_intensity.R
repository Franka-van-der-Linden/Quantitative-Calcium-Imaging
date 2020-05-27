## This script was used to determine the pKa for Tq-Ca-FLITS and jGCaMP7c from intensity data

## Part 1:contains
# Fitting of model to the measurement data
# Calculation of the dynamic range vs pH

## Part 2: contains
# plotting the models with the measured data


# add predictions
require(modelr)
#plotting data
require(ggplot2)

link1 <- "https://raw.githubusercontent.com/Franka-van-der-Linden/Quantitative-Calcium-Imaging/master/R_scripts_for_modelling/pKa_FLITS_Ca.csv"
link2 <- "https://raw.githubusercontent.com/Franka-van-der-Linden/Quantitative-Calcium-Imaging/master/R_scripts_for_modelling/pKa_FLITS_EGTA.csv"
link3 <- "https://raw.githubusercontent.com/Franka-van-der-Linden/Quantitative-Calcium-Imaging/master/R_scripts_for_modelling/pKa_jGCaMP7c_Ca.csv"
link4 <- "https://raw.githubusercontent.com/Franka-van-der-Linden/Quantitative-Calcium-Imaging/master/R_scripts_for_modelling/pKa_jGCaMP7c_EGTA.csv"
link5 <- "https://raw.githubusercontent.com/Franka-van-der-Linden/Quantitative-Calcium-Imaging/master/R_scripts_for_modelling/pKa_mTurquoise2.csv"
FLITS_Ca <- read.csv(link1)
FLITS_EGTA <- read.csv(link2)
jGCaMP7c_Ca <- read.csv(link3)
jGCaMP7c_EGTA <- read.csv(link4)
mTurquoise2 <- read.csv(link5)


###### PART 1: fitting of the pKa curves.

## function for fitting one pKa value
model_pH <- function(pH, intensity, pKa.init, n.init) {
  max.init=max(intensity)
  min.init=0
  startvals = c(pKa = pKa.init, max=max.init, min=min.init, n=n.init)
  model <- nls(intensity ~ min + (max-min) / ((10^(pKa-pH))^n + 1), 
               start = startvals, lower = list(min = 0), algorithm = "port", model = TRUE,trace=FALSE)
  return(model)
}

FLITS_EGTA_model <- model_pH(FLITS_EGTA$pH, FLITS_EGTA$Intensity, pKa.init = 4, n.init = 1)
jGCaMP7c_Ca_model <- model_pH(jGCaMP7c_Ca$pH, jGCaMP7c_Ca$Intensity, pKa.init = 4, n.init = 1)
jGCaMP7c_EGTA_model <- model_pH(jGCaMP7c_EGTA$pH, jGCaMP7c_EGTA$Intensity, pKa.init = 4, n.init = 1)
mTurquoise2_model <- model_pH(mTurquoise2$pH, mTurquoise2$Intensity, pKa.init=4, n.init=1)

## function for fitting two pKa value's
model_pH2 <- function(pH, intensity, pKa1.init, pKa2.init, n1.init, n2.init) {
  max.init=max(intensity)
  med.init=max.init/10
  min.init=0
  startvals = c(pKa = pKa1.init,pKa2 = pKa2.init, max=max.init, med=med.init, min=min.init, n=n1.init, n2=n2.init)
  model <- nls(intensity ~ min + (med-min) / ((10^(pKa-pH))^n + 1) + (max-med) / ((10^(pKa2-pH))^n2 + 1), 
               start = startvals, lower = list(min = 0), algorithm = "port", model = TRUE, trace=FALSE)
  return(model)
}

FLITS_Ca_model <- model_pH2(FLITS_Ca$pH, FLITS_Ca$Intensity, pKa1.init = 4, pKa2.init=6, n1.init = 1, n2.init=1)

#store the parameters to view
parameters <- data.frame(FLITS_Ca=coef(FLITS_Ca_model), FLITS_EGTA=NA, 
                         jGCaMP7c_Ca=NA, jGCaMP7c_EGTA=NA, 
                         mTurquoise2=NA)
parameters[c("pKa", "max", "min", "n"),2:5] <- c(coef(FLITS_EGTA_model),
                                                 coef(jGCaMP7c_Ca_model),
                                                 coef(jGCaMP7c_EGTA_model),
                                                 coef(mTurquoise2_model))

#make dataframe for a smooth fit and calculate the dynamic range
pH_range <- seq(2.7,10.2,0.1)

FLITS_fit <- add_predictions(data.frame(pH = pH_range), FLITS_Ca_model, var="Ca")
FLITS_fit <- add_predictions(FLITS_fit, FLITS_EGTA_model, var= "EGTA")
FLITS_fit$dynamic_range <- FLITS_fit$Ca / FLITS_fit$EGTA

jGCaMP7c_fit <- add_predictions(data.frame(pH = pH_range), jGCaMP7c_Ca_model, var="Ca")
jGCaMP7c_fit <- add_predictions(jGCaMP7c_fit, jGCaMP7c_EGTA_model, var="EGTA")
jGCaMP7c_fit$dynamic_range <- jGCaMP7c_fit$Ca / jGCaMP7c_fit$EGTA

mTurquoise2_fit <- add_predictions(data.frame(pH = pH_range), mTurquoise2_model, var="Intensity")


###### PART 2: plotting of the models
#plot the models with the data
ggplot(FLITS_fit, aes(x=pH)) +
  geom_point(data=FLITS_EGTA, aes(y=Intensity), color="turquoise4", alpha=0.5, shape=16) +
  geom_point(data=FLITS_Ca, aes(y=Intensity), color="turquoise1", alpha=0.5, shape=16) +
  geom_line(aes(y=EGTA), color="turquoise4") +
  geom_line(aes(y=Ca), color="turquoise1") +
  ylab("Intensity (a.u.)") +
  ggtitle("Tq-Ca-FLITS") +
  theme_light()

ggplot(jGCaMP7c_fit, aes(x=pH)) +
  geom_point(data=jGCaMP7c_Ca, aes(y=Intensity), color="green1", alpha=0.5, shape=16) +
  geom_point(data=jGCaMP7c_EGTA, aes(y=Intensity), color="green4", alpha=0.5, shape=16) +
  geom_line(aes(y=EGTA), color="green4") +
  geom_line(aes(y=Ca), color="green1") +
  ylab("Intensity (a.u.)") +
  ggtitle("jGCaMP7c") +
  theme_light()

ggplot(mTurquoise2_fit, aes(x=pH, y= Intensity)) +
  geom_point(data=mTurquoise2, color="cyan", alpha=0.5, shape=16) +
  geom_line(color="cyan") +
  ylab("Intensity (a.u.)") +
  ggtitle("mTurquoise2") +
  theme_light()

#plot the dynamic range
ggplot(FLITS_fit, aes(x=pH, y=dynamic_range)) + 
  geom_line(aes(y=dynamic_range), color="gray") +
  ylab("Dynamic range") +
  ggtitle("Tq-Ca-FLITS") +
  theme_light()

ggplot(jGCaMP7c_fit, aes(x=pH, y=dynamic_range)) + 
  geom_line(aes(y=dynamic_range), color="gray") +
  ylab("Dynamic range") +
  ggtitle("jGCaMP7c") +
  theme_light()

