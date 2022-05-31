
#################################################
#Code by Nasim Khayyer
#Code by Nasim Khayyer - M1 DASEE - Spatial point pattern analysis
#################################################

library(spatstat)
library(dbmss)
load("./CSBIGS.Rdata")

# graph data
plot(Pop, cols = "red", use.marks = TRUE, legend = TRUE, markscale = 0.4, axes=TRUE, main = "Full data - CSBIGS")

#About data
cat("Number of points in dataset is:" , npoints(Pop))
#marks(Pop)
#coords(Pop)
as.owin(Pop)
#df <- as.data.frame(Pop)

#Reduce Scale
Pop_rescale <- affine(Pop, mat=diag(c(1/1000,1/1000)))
plot(Pop_rescale, cols = "red", use.marks = TRUE, legend = TRUE, markscale = 0.0003, axes=TRUE)

#Subset data to a study subregion for fitting models
subregion<- owin(c(515000, 535000), c(1835000, 1855000))
region <- owin(c(513273.8, 540696), c(1825747.7, 1858266.6))
plot(Pop[subregion],cols = "red", use.marks = TRUE, legend = TRUE, markscale = 0.15, axes=TRUE, main = "study subregion of CSBIGS")

# datasummary
Pop_summary <- summary(Pop)
Pop_summary

############################### Density plot Analysis #########################
#Density plots
plot(density(Pop,500), main = "density of observations")
plot(density(Pop,1100), main = "density of observations")
plot(density(Pop[subregion],1100), main = "density of observations - subregion")

#Contour plot
contour(density(Pop,400), axes = F, main = "contour plot of CSBIGS")
contour(density(Pop,500), axes = F, main = "contour plot of CSBIGS")
contour(density(Pop[subregion],400), axes = F, main = "contour plot of CSBIGS - SUBREGION")

# Kernel density analysis
den_400 <- density(Pop,kernel = "epanechnikov",sigma=400)  
den_500 <- density(Pop,kernel = "epanechnikov",sigma=500)  
plot(den_400)
plot(den_500)
plot(Pop, add = TRUE, cols = "green", markscale = 0.3)

#3D representation of data
persp(den_500, col = "red", main = "3D representation of kernel density")

#Density in the selected sub region of CSBIGS
den_500_subregion <- density(Pop[subregion],kernel = "epanechnikov",sigma=500)  
plot(den_500_subregion, main = "kernel density for subregion")
plot(Pop[subregion], add = TRUE, cols = "green", markscale = 0.3, )
persp(den_500_subregion, col = "red", main = "3D representation of intensity - subregion")


##########################  INTENSITY ANALYSIS ########################
lamb <- summary(Pop)$intensity
cat( "Point intensity is equal to:", lamb )

#Density analysis of data
den_500 <- density(Pop,kernel = "epanechnikov",sigma=500)  # The resulting object is a pixel image. This class has methods for print, summary, plot, contour and persp
den_800 <- density(Pop,kernel = "epanechnikov",sigma=800)  # The resulting object is a pixel image. This class has methods for print, summary, plot, contour and persp
plot(den_500)
plot(den_800)
plot(Pop, add = TRUE, cols = "white", markscale = 0.3)

boxplot(den_500$v)

#Quadrat counting map  of CSBIGS and its subregion
quadrat_Pop <- quadratcount(Pop, nx = 7, ny = 7)
plot(quadrat_Pop, main = "quadrat counting the data")
quadrat_Pop <- quadratcount(Pop[subregion], nx = 7, ny = 7)
plot(quadrat_Pop, main = "quadrat counting the data - Subregion")

#Quadrat test
q_test10 <- quadrat.test(Pop[subregion], nx = 10, ny = 10)
q_test10
q_test15 <- quadrat.test(Pop[subregion], nx = 15, ny = 15)
q_test15
q_test20 <- quadrat.test(Pop[subregion], nx = 20, ny = 20)
q_test20

#iNTENSITY ANALYSIS BASED ON DISTANCE OF FIRES
plot(distmap(Pop))
plot(rhohat(ppp(Pop$x, Pop$y, subregion), distmap(Pop)), main ="intensity of events based on distance")

#kest
kest_data <- Kest(Pop)
kest_data
plot(kest_data)

#envelope
E <- envelope(Pop[subregion],Kest,nsim=39)
plot(E)

#Quadrat test for intensity analysis
q_test10 <- quadrat.test(Pop[subregion], nx = 10, ny = 10)
q_test10
q_test15 <- quadrat.test(Pop[subregion], nx = 15, ny = 15)
q_test15
q_test20 <- quadrat.test(Pop[subregion], nx = 20, ny = 20)
q_test20
plot(subregion)
plot(q_test10, add= TRUE)

############## Intensity analysis using Poisson models and MAXIMUM LIKELIHOOD ESITMATIONS

# MODEL SPECIFCATION FOR POISSON PROCESS
fit_1 <- ppm(ppp(Pop$x, Pop$y, subregion),~1)  
fit_1
fit_2 <- ppm(ppp(Pop$x, Pop$y, subregion),~x+y)  
fit_2
lam_2 <- predict(fit_2)
plot(lam_2, main = "fitted model of the Possion MLE - first order")
plot(Pop[subregion],  pch = ".", cols="white", add = TRUE)
#plot(effectfun(fit_2), main = "effect function of")
fit_3 <-ppm(ppp(Pop$x, Pop$y, subregion),~polynom(x,y,2))
fit_3
lam_3 <- predict(fit_3)
plot(lam_3, main = "fitted model of the Possion MLE - polynomial power 2")
plot(Pop[subregion], pch = ".", cols="white", add = TRUE)

#polynomial model power 4
fit_4 <-ppm(ppp(Pop$x, Pop$y, subregion),~polynom(x,y,4))
fit_4
lam_4 <- predict(fit_4)
plot(lam_4,  main = "fitted model of the Possion MLE - polynomial power 4")
plot(Pop[subregion],  pch = ".", cols="white", add = TRUE)
plot(effectfun(fit_2,"x", y=0.5))
plot(effectfun(fit_2,"y", x=0.5))

# MODEL WITH DISTANCE AS COVARIATE with power 2
z_distance <- distmap(Pop[subregion])
fit_5 <- ppm(ppp(Pop$x, Pop$y, subregion),~polynom(z_distance,2))
fit_5
lam_5 <- predict(fit_5)
plot(lam_5,  main = "fitted model of the Possion MLE - distance model")
plot(Pop[subregion],  pch = "*", cols="white", add = TRUE)
plot(effectfun(fit_5), main = "intensity as function of distance")

#ANOVA TEST
anova(fit_1, fit_4, test = "Chi")

#########################   Goodness of Fit ##############################

#quadrat MODEL FOR SUB REGION - LINEAR
fit_2 <- ppm(ppp(Pop$x, Pop$y, subregion),~x+y)  
fit_2
quadrat_fit_2 <- quadrat.test(fit_2, nx = 7, ny = 7)
quadrat_fit_2
plot(Pop[subregion], pch = "*", main = "linear model")
plot(quadrat_fit_2, add = TRUE, col = "red")
#quadrat model for sub region - polynomial power 4
fit_3 <-ppm(ppp(Pop$x, Pop$y, subregion),~polynom(x,y,2))
fit_3
quadrat_fit_3 <- quadrat.test(fit_3, nx = 7, ny = 7)
quadrat_fit_3
plot(Pop[subregion], pch = "*", main = "polynomial power 4 model")
plot(quadrat_fit_3, add = TRUE, col = "red")
#quadrat model for sub region - polynomial power 4
fit_4 <-ppm(ppp(Pop$x, Pop$y, subregion),~polynom(x,y,4))
fit_4
quadrat_fit_4 <- quadrat.test(fit_4, nx = 7, ny = 7)
quadrat_fit_4
plot(Pop[subregion], pch = ".", main = "polynomial power 4 model")
plot(quadrat_fit_4, add = TRUE, col = "red")

#Residual diagnostics of the linear model and polynomial model 
diagnose.ppm(fit_2)
diagnose.ppm(fit_4)
diagnose.ppm(fit_2, which = "smooth", main = "Smoother residuals for first order model")
diagnose.ppm(fit_4, which = "smooth", main = "Smoother residuals for power 4 polynomial model")

z_distance <- distmap(Pop[subregion])
lurking(fit_4, z_distance, type = "raw", main = "lurking variable plot for distance - before adding distance")

fit_5 <-ppm(ppp(Pop$x, Pop$y, subregion),~polynom(x,y,4) + polynom(z_distance,2))
fit_5
diagnose.ppm(fit_5, which = "smooth", main = "Smoother residuals for power 4 polynomial model + distance")
lurking(fit_5, z_distance, type = "raw", main = "lurking variable plot for distance + after adding distance")

quadrat_fit_5 <- quadrat.test(fit_5, nx = 7, ny = 7)
quadrat_fit_5


################################# POINT PATTERN DEPENDENCE ANALYSIS ################
#MORISHITA PLOT
miplot(ppp(Pop$x, Pop$y, subregion), main = "Morishita plot of point dependence in CLM subregion")
# FRYPLOT
fryplot(ppp(Pop$x, Pop$y, subregion), axes = TRUE, main = "fryplot for CSBIGS - SUBREGION")

#Pair distances
d_1 <- nndist(ppp(Pop$x, Pop$y, subregion))
d_2 <- pairdist(ppp(Pop$x, Pop$y, subregion))
d_3 <- pairdist(ppp(Pop$x, Pop$y, subregion), periodic=TRUE)
d_4 <- pairdist(ppp(Pop$x, Pop$y, subregion), squared=TRUE)

#POINT DEPENDENCE ANALYSIS USING K FUNCTION
#Kest
Gc <- Kest(ppp(Pop$x, Pop$y, subregion))
Gc
par(pty = "s")
plot(Gc)
#L (R)
L <- Lest(ppp(Pop$x, Pop$y, subregion))
plot(L, main = "L function")
# ENVELOPE CRITICAL THRESHOLDS
E <- envelope(ppp(Pop$x, Pop$y, subregion), Kest, nsim = 39, rank = 1)
E
plot(E, main = "pointwise envelope")
E <- envelope(ppp(Pop$x, Pop$y, subregion), Lest, nsim = 39, rank = 1, global = TRUE)  # This stabilizes the variance
E
plot(E, main = "global envelopes of L(r)")
v <- varblock(ppp(Pop$x, Pop$y, subregion), Kest, nx=4, ny=2) # CI for the K function
plot(v)

#### NON POISSON MODELS
fit_10 <- kppm(ppp(Pop$x, Pop$y, subregion),~1,"MatClust")
fit_10
plot(fit_10)
E <- envelope(fit_10,Kest,nsim= 59, global = TRUE, correction = "border")
plot(E)
######################### mark analysis ########################
summary(Pop)
hist(Pop$marks, breaks = 100, main = "histogram of CSBIGS Marks")

#model with marks as covariates
fit_mark_1 <- ppm(ppp(Pop$x, Pop$y, subregion, marks = Pop$marks),~marks + polynom(x,y,2))
fit_mark_1
lam_mark_1 <- predict(fit_mark_1)
plot(lam_mark_1, main = "fitted model with marks_1")

#model with marks as covariates
fit_mark_1 <- ppm(ppp(Pop$x, Pop$y, subregion, marks = Pop$marks),~marks + polynom(x,y,4))
fit_mark_1
lam_mark_1 <- predict(fit_mark_1)
plot(lam_mark_1, main = "fitted model with marks_1")



