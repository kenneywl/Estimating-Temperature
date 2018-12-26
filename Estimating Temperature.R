###################################################
#Part 0: Description.
#The goal of the following project is to forcast
#global temperature just a few years. 2006-2011.
#After preprocessing we will endup with
#two data frames. Data frame "warm" is composed of
#year,temp (temperature),co2 (carbon dioxide),ch4 (methane),
#n20 (dinitrogen monoxide),cfc (chlorofluorocarbon).
#from 1880 to 2005.
#Data frame "carbon" is composed of
#year,temp, and co2 from 2006 to 2011.
#In this analysis we try to predict the temperatures
#just in data frame "carbon" (6 data points).
#####################################################
#Our data comes from GISTEMP Team, 
#2016: GISS Surface Temperature Analysis 
#(GISTEMP). NASA Goddard Institute for Space Studies
#Temperature is from the Global
#Land Ocean Index. The GLOI is the deviation from the 
#thirty year average global temperature between 1951 and 1980. 
#That temperature is estimated to be 14C.
#The factor "temp" is global deviation from this value.
#i.e temp=0 => temp=14c.
#Global temperatures were computed by GISS.
#http://data.giss.nasa.gov/gistemp/graphs_v3/
#Carbon and other green house gasses were made by artic 
#ice core samples being taken and analyized.
#http://data.giss.nasa.gov/modelforce/ghgases/
#This data was interpolated and "smoothed" by GISS.
###########################################################
#This is an excersise in linear regression with
#principle compenent analysis. I also fit a polynomial on time,
#just to see what I could get to work. Overall, youll see that
#simple linear regression is probabily best.

#####################################################
#Part 1: DATA INPUT
#If you follow along with this in your own terminal,
#you need to change the working directory to wherever
#the files are kept.
setwd("/home/kenneywl/Documents/Global Warming/Estimating Temperature")
#I do very little commenting in this section because
#it is mostly just all preprocessing and getting the
#data in a usable form. (No analysis.)
######################################################
#We're off!
#global temp average 1880 to 2011
temp <- read.fwf("temp1880.txt",widths=c(5,11,11))
temp <- temp[-1,]
colnames(temp) <- c("Year","Annual.Mean","Five.Year.Mean")
temp <- data.frame(temp)

temp <- sapply(temp,FUN=function(x) trimws(x,"both"))
temp[temp=="*"] <- NA
temp <- apply(temp,MARGIN=2,FUN=function(x) as.numeric(x))

temp <- temp[-133:-137,]
rownames(temp) <- 1:132

#co2 data 1880 to 2011
co2 <- read.fwf("icecore.txt",widths=c(7,6,6,14,6,9,6,7,6,6,9))
co2 <- sapply(co2,FUN=function(x) trimws(x,"both"))
co2 <- co2[,-c(1,6,9)]
co2 <- apply(co2,MARGIN=2,FUN=function(x) as.numeric(x))
co2 <- data.frame(co2)
colnames(co2) <- c("year","co2.ppm")

co2[51:100,1:2] <- co2[1:50,3:4]
co2[101:150,1:2] <- co2[1:50,5:6]
co2[151:200,1:2] <- co2[1:50,7:8]
co2 <- co2[,-3:-8]
co2 <- co2[-(163:200),]
co2 <- co2[-1:-30,]

#ch4 1850 to 2005
ch4 <- read.fwf("ch4.txt",widths=c(5,8))
colnames(ch4) <- c("year","ch4.ppm")

#cfc's 1850 to 2005

cfc <- read.fwf("cfc.txt",widths=c(5,7))
cfc <- data.frame(cfc)
colnames(cfc) <- c("year","cfc.ppm")

#n2o 1850 to 2005
n2o <- read.fwf("n2o.txt",widths=c(5,8))
n2o <- data.frame(n2o)
colnames(n2o) <- c("year","n2o.ppm")

#temp by carbon only 1880 to 2011 in df "carbon"
carbon <- data.frame(temp[127:132,c(1,2)],co2[127:132,2])
colnames(carbon) <- c("year","temp","co2")

#combine temp+co2ppm+ch4+n2o+cfc's 1880 to 2005
warm <- data.frame(year=co2[1:126,1],temp=temp[1:126,2],
                   co2=co2[1:126,2],ch4=ch4[31:156,2],n2o=n2o[31:156,2],cfc=cfc[11:136,2])
rownames(warm) <- 1:nrow(warm)

#Available for use "warm" and "carbon".




#############################################################################
#Part 2: MLR
#We have all the green house data and temperature from 1880 to 2005 in "warm"
#and the carbon data and temperature from 2006 to 2011 in "carbon".
#We aim to use "warm" to build a model to predict the temperature in "carbon"
##############################################################################
#First lets take a look at our data:

pairs(warm[,3:6],main="Linear Correlation of Factors")

#Everything is linearly correlated!

#Our aim is to predict 2006 to 2011 temperature.
#We have the carbon data for those years, but not
#the other green house gasses.
#We want to use "warm" to predict the the temp in "carbon".
#First we'll fit a multiple linear regression (MLR).
#Problem with this data set is that all the predictors
#are correleted. Well use PCA to overcome this:

pc <- princomp(warm[,3:6])

#Nifty way to transform the points into orthogonal factors by the eigenvectors from PCA.
pc <- apply(pc$loadings[1:4,1:4],MARGIN=2,FUN=function(x)
  as.matrix(warm[,c(3,4,5,6)]) %*% as.matrix(x))
pc <- data.frame(pc)
warm_pc <- data.frame(year=warm$year,temp=warm$temp,comp1=pc$Comp.1,
                      comp2=pc$Comp.2,comp3=pc$Comp.3,comp4=pc$Comp.4)


#Fit a model with all components
#note these factors are orthoganal by PCA
#so there are no interaction terms.

l <- lm(temp~comp1+comp2+comp3+comp4,data=warm_pc)
summary(l)
anova(l)

#Comp 2 and 3 are not at all significant so we remove them and refit.

l <- lm(temp~comp1+comp4,data=warm_pc)
summary(l)
anova(l)

#Our adjusted r-squared goes up slightly, but what is important is we have a simpler model.
#Also note HOW significant comp1 and comp4 are. VERY. A pvalue of .05 is a 
#1/.05=20, 1:20 chance that the data (or more extreme) would be observed given 
#no relationship between predictors and response. A 2*10^-16 is a 1/(2*10^-16)=5*10^15, 
#a 1:5*10^15 chance we see this or more extreme.

#Model adequecy checking:
r <- resid(l)

#Normality of residuals and for outliers
#Residuals are normal enough. We may have one technical outlier.
#But it should be kept in as it doesn't affect
#the results significantly.
par(mfrow=c(1,2))
hist(r,main="Histogram of Residuals")
qqnorm(r);qqline(r)

#Just to have on the record. These are possible outliers.
warm_pc[which(cooks.distance(l)>4/nrow(warm_pc)),1]

#constant variance and independence of errors.
par(mfrow=c(2,2))
plot(r~l$fitted.values,main="Residual vs Fitted")
plot(r~comp1,data=warm_pc,main="Residual vs Comp1")
plot(r~comp4,data=warm_pc,main="Residual vs Comp4")
plot(r~year,data=warm_pc,main="Residual vs Year")

#in about 1940 there are some correlation.
#not enough to invalidate, (I'll say)
#through a little research I found that it is
#not that there was some one time event causing rapid warming,
#but that a high concentration of sulphate aerosols in the atmosphere
#countered the warming trend causing a cooling of what would have been exponential
#warming. Anyway, the deviation is not so pronounced
#as to invalidate the model. What should be done (but isn't
#in this analysis) is obtain data on sulphate aerosols
#and make them a factor in the model.


#And lets take a look at our model, plotted against time.
pred <- l$fitted.values
actual <- warm_pc$temp

par(mfrow=c(1,1))
plot(actual~warm_pc$year,col="red",main="MLR with PCA",ylab="Temperature Differential",xlab="Year")
points(pred~warm_pc$year,col="blue",type="l")
legend(1880,.6,legend=c("Data Points","Regression Line"),lwd=c(2.5,2.5),col=c("red","blue"))




#####################################################################
#Part 3: SLR 
#We will build a SLR using just co2
#####################################################################
#lets look at correlation again.
pairs(warm[,c(2,3)],main="Linear Correlation of Factors")

#These are very linearly correleted.
#lets make a fit

l2 <- lm(temp~co2,data=warm)
summary(l2)
anova(l2)

#Notice adjusted r-squared is slightly less then the MLR PCA version
#.81 compared to .83. MSE is .0137 compared to .0122.
#It appears that this model is just a hair worse then the MLR version.
#It is imporant to note that the majority of predicting power of the MLR
#version is contained just in the carbon data.

#Model adequecy checking.
r2 <- resid(l2)

#normality of residuals #check! There may be an outlier
#all the way to the right. But the Quantile
#and Shapiro Wilks test both look great.
par(mfrow=c(1,2))
qqnorm(r2);qqline(r2)
hist(r2,main="Histogram of Residuals")
shapiro.test(r2)


#Homoscedasticity. There appears to be a slight decrease
#in variance. This is proababily due to increased sensitivity
#of equipment. Independence of errors. OK.
par(mfrow=c(1,2))
plot(r2~l2$fitted.values,main="Residual vs Fitted")
plot(r2~warm$temp,main="Residual vs Temp")

#lets take a look at the model plotted alongside the MLR version.
pred <- l$fitted.values
pred2 <- l2$fitted.values
actual <- warm_pc$temp

par(mfrow=c(1,1))
plot(actual~warm_pc$year,col="red",main="MLR vs. SLR",ylab="Temperature Differential",xlab="Year")
points(pred~warm_pc$year,col="blue",type="l")
points(pred2~warm_pc$year,col="black",type="l")
legend(1880,.6,legend=c("MLR Regression","SLR  Regression"),lwd=c(2.5,2.5),col=c("blue","black"))

#It is clear to see why MLR does slightly better,
#it appears to be slightly more sensitive.


########################################################################################
#Part 4: Imputation
#Next we take a look at "carbon." This database has year,temp, and carbon from 2006 to 2011.
#Our goal is to predict those 6 data points.
#We can use the SLR directly, because we have the carbon data for those years.
#But we don't have any of the other greenhouse gasses.
#We'll try to impute them here.
#######################################################################################

#Lets get a feel for the data.
pairs(warm[,3:6])

#lets look at just the last portion of our data set:
last13 <- 114:126
pairs(warm[last13,3:6])

#Since we need to forcast only 6 points into the future,
#and the last 13 points are fairly linear with carbon,
#well do SLR against each independantly.
#the hope is here that the extra data from correlated with carbon
#gives us more information than carbon alone. Imputation only goes so far.

lch4 <- lm(ch4~co2,data=warm[last13,])
anova(lch4)

#Model adequecy checking.
r3 <- resid(lch4)
par(mfrow=c(2,2))
qqnorm(r3);qqline(r3)
hist(r3,main="Histogram of Residuals")
plot(r3~lch4$fitted.values,main="Residual vs Fitted")
plot(r3~warm[last13,3],main="Residual vs Co2")

#Everything is fine except residuals vs fitted. This will decrease our accuracy.

ln2o <- lm(n2o~co2,data=warm[last13,])
anova(ln2o)

#Model adequecy checking.
par(mfrow=c(2,2))
r4 <- resid(ln2o)
qqnorm(r4);qqline(r4)
hist(r4,main="Histogram of Residuals")
plot(r4~ln2o$fitted.values,main="Residual vs Fitted")
plot(r4~warm[last13,3],main="Residual vs Co2")

#We press on despite some devitations. (specifically residual vs. fitted and
#residual vs Co2 look bad.)

lcfc <- lm(cfc~co2,data=warm[last13,])
anova(lcfc)

#Model adequecy checking.
par(mfrow=c(2,2))
r5 <- resid(lcfc)
qqnorm(r5);qqline(r5)
hist(r5,main="Histogram of Residuals")
plot(r5~ln2o$fitted.values,main="Residual vs Fitted")
plot(r5~warm[last13,3],main="Residual vs Co2")

#We use our models to predict.

carbon$ch4 <- predict(lch4,carbon)
carbon$n2o <- predict(ln2o,carbon)
carbon$cfc <- predict(lcfc,carbon)

#We take the eignvectors from our PCA from "warm" and use them in "carbon"
#to get our corresponding compenents.

pc <- princomp(warm[,c(3,4,5,6)])
pc <- apply(pc$loadings[1:4,1:4],MARGIN=2,FUN=function(x) 
  as.matrix(carbon[,c(3,4,5,6)]) %*% as.matrix(x))
pc <- data.frame(pc)
carbon_pc <- data.frame(year=carbon$year,temp=carbon$temp,comp1=pc$Comp.1,
                        comp2=pc$Comp.2,comp3=pc$Comp.3,comp4=pc$Comp.4)

#and we have imputed our principle compenents.

#####################################################################
#Part 5: Linear Regression Testing
#Now we predict 2006 to 2011 temperature and compare.
####################################################################

mlr_predict <- predict(l,carbon_pc)
slr_predict <- predict(l2,carbon)

#MSE for both  
mlr_mse <- mean((mlr_predict-carbon$temp)^2)
slr_mse <- mean((slr_predict-carbon$temp)^2)
mlr_mse
slr_mse
#not terrible. slr is slightly better by .00115
mlr_mse-slr_mse

#lets look at the plots:
par(mfrow=c(1,1))
plot(temp~year,data=carbon,ylim=c(.5,.75),main="MLR and SLR Regression Lines")
lines(slr_predict~carbon$year,col="blue")
lines(mlr_predict~carbon$year,col="red")

legend(2006,.55,legend=c("MLR","SLR"),lwd=c(2,2),col=c("red","blue"))

#and lets look at the whole thing together ploted against year again.

all_data <- rbind(warm,carbon)
all_slr <- c(l2$fitted.values,slr_predict)
all_mlr <- c(l$fitted.values,mlr_predict)

plot(temp~year,data=all_data,main="MLR vs SLR")
points(all_slr~all_data$year,col="blue",type="l")
points(all_mlr~all_data$year,col="red",type="l")
abline(v=2005)

legend(1880,.6,legend=c("MLR","SLR"),lwd=c(2.5,2.5),col=c("red","blue"))

#The vertical line separates what was used for training and what was used to test.
#We can see that both are pretty good, SLR might be just a hair better,
#that is because we actually have the data for SLR to make the prediction.
#For MLR, we had to find a way to impute the values for the principle compentents.

#The final result is usually better the less one has to impute, that is because
#imputation is just using the data we already have to "stand in" for missing data points.
#As could be expected, actual data is more informative then no data.

#As a final note, temperature forcasting is increadably complicated, with many
#factors that come into play. These factors go way beyond this analysis.
#But it is amazing what we CAN do with simple models
#Provided we just pick out what is most important.

############################################################################
#Part 6: Polynomial regression. (Bonus section)
###########################################################################
#An interesting thing about polynomial regression is that
#you can get a better and better fit just by increasing the
#order of the polynomial,but at some point you are overfitting.
#It is also a common mantra to not extrapolate into the future
#when you do polynomial regression. Unfortunately, this is exactly
#what forecasting is. There is no way around it. The future is nebulus,
#but hopefully the past can tell us something about it.
#For the following analysis I only use time as a predictor.
############################################################################
#Lets fit a polynomial model.

p <- lm(temp~poly(year,2),data=warm)
summary(p)

#We stop at the 2nd order polynomial because the 3rd is not significant.
#We have an adjsuted r squared of .80, which is on par with the other models,
#although just a little worse. 
#Here is a plot, I put in the the confidence intervals too.
temp_year <- predict(p,all_data,interval="confidence",level=.95)
temp_year <- cbind(data.frame(year=1880:2011),temp_year)
plot(temp~year,data=all_data,main="2nd Order Polynomial Regression")
lines(fit~year,data=temp_year,col="red")
lines(lwr~year,data=temp_year,col="blue")
lines(upr~year,data=temp_year,col="blue")
abline(v=2005)

#It appears to be a little too conservative. If we eye the derivative of post 1980
#temperature change, we get a sharper incline. The high end of our confidence interval
#looks much closer to what we want. The MSE for 2005 to 2011 is:

mean((temp_year[127:132,2]-carbon[,2])^2)

#Surprisingly, MSE for these points is about equal to the SLR case with carbon only.
#Take a look at the zoomed in graph of SLR and Poly.

par(mfrow=c(1,1))
plot(temp~year,data=carbon,ylim=c(.5,.75),main="SLR and Polynomial Regression Lines")
lines(slr_predict~carbon$year,col="blue")
lines(fit~year,data=temp_year[127:132,],col="green")
legend(2006,.55,legend=c("SLR","Poly"),lwd=c(2,2),col=c("blue","green"))

#Model adequecy (quick version):
par(mfrow=c(1,2))
hist(resid(p),main="Histogram of Residuals")
plot(resid(p)~fitted(p),main="Residuals vs Fitted")

#Just like linear regressinon, it doesn't do a good job at passing Residual vs Fitted.

#############################################################################################
#All told, the real lesson here is that a simple analysis goes far. We might be able to
#estimate global cooling with an indicator variable pre 1940 and post 1940, assuming the actual
#numerical values of sulphate aerosols are not procurable.

