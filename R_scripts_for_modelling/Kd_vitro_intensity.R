## This script was used to determine the in vitro Kd for jGCaMP7c and Tq-Ca-FLITS from intensity data
# the script fits a curve to the data to determine the Kd.
# also some plots are generated, showing the model and the measurements

#for plotting
require(ggplot2)
#add predictions
require(modelr)

# #import data
link1 <- "https://raw.githubusercontent.com/Franka-van-der-Linden/Quantitative-Calcium-Imaging/master/R_scripts_for_modelling/Kd_vitro_intensity_Tq-Ca-FLITS.csv"
FLITS <- read.csv(link1)
link2 <- "https://raw.githubusercontent.com/Franka-van-der-Linden/Quantitative-Calcium-Imaging/master/R_scripts_for_modelling/Kd_vitro_intensity_jGCaMP7c.csv"
jGCaMP7c <- read.csv(link2)

## I want to fit a model for the Kd, without normalization of the data

#startvalues for fitting of Tq-Ca-FLITS data
Ka.init <- 0.3
n.init <- 1
L <- FLITS$concentration
Fb <- FLITS$Intensity
max.init <-5000
min.init <-1
# fitting equation
model_FLITS <- nls(Fb ~ ((max-min) / ((Ka/L)^n + 1)+min), start = c(Ka = Ka.init,n = n.init, max=max.init, min=min.init),model = TRUE)

#startvalues for fitting of jGcaMP7c data
Ka.init <- 0.3
n.init <- 1
L <- jGCaMP7c$concentration
Fb <- jGCaMP7c$Intensity
max.init <-5000
min.init <-1
# fitting equation
model_jGCaMP7c <- nls(Fb ~ ((max-min) / ((Ka/L)^n + 1)+min), start = c(Ka = Ka.init,n = n.init, max=max.init, min=min.init),model = TRUE)


#store parameters
(param_FLITS <- coef(model_FLITS))
(param_jGCaMP7c <- coef(model_jGCaMP7c))
#calculate confidence interval (95%) for the models
(ci_FLITS <- confint(model_FLITS))
(ci_jGCaMP7c <- confint(model_jGCaMP7c))

## make dataframe for a smooth fit
# Tq-Ca-FLITS
conc_range <- c(0, 10^seq(-2,1.7,0.1))
fit_FLITS<- data.frame(L = conc_range) 
fit_FLITS <- add_predictions(fit_FLITS, model_FLITS)
names(fit_FLITS)[1] <- "concentration"

# jGCaMP7c
conc_range <- c(0, 10^seq(-2,1.7,0.1))
fit_jGCaMP7c <- data.frame(L = conc_range) 
fit_jGCaMP7c <- add_predictions(fit_jGCaMP7c, model_jGCaMP7c)
names(fit_jGCaMP7c)[1] <- "concentration"

## plot individual measurements and model
# Tq-Ca-FLITS
ggplot(FLITS, aes(x=concentration, y=Intensity)) + 
  geom_point(alpha=0.5) + 
  geom_line(data=fit_FLITS, aes(y=pred)) + 
  scale_x_log10() +
  ggtitle("Tq-Ca-FLITS") + xlab("Free calcium (uM)") + ylab("Fluorescence intensity (a.u.)")

# jGCaMP7c
ggplot(jGCaMP7c, aes(x=concentration, y=Intensity)) + 
  geom_point(alpha=0.5) + 
  geom_line(data=fit_jGCaMP7c, aes(y=pred)) + 
  scale_x_log10() +
  ggtitle("jGCaMP7c") + xlab("Free calcium (uM)") + ylab("Fluorescence intensity (a.u.)")
