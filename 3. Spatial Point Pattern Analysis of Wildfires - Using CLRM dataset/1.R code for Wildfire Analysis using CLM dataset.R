

#################################################
#Code by Nasim Khayyer
#Code by Nasim Khayyer - M1 DASEE - Spatial point pattern analysis
#################################################
library(spatstat)
data("clmfires")
data <- clmfires

#Plotting data by cause of fire
plot(data, which.marks="cause", cols=2:5, cex=1,
     main="Castilla-La Mancha forest fires - by cause", axes = TRUE)

cat("Number of points in dataset is:" , npoints(data))
#marks(data)
#coords(data)
as.owin(data)
df <- as.data.frame(data)

#select a region in CLS
selected_data <- owin(c(200, 300), c(200, 300))
plot(data[selected_data], which.marks="cause", cols=2:5, cex=1,
     main="Castilla-La Mancha forest fires for selected coordinates - by cause")

clm_region <- owin(c(8.247, 385.345), c(24.220, 377.177))
clm_subregion <- owin(c(170, 350), c(70, 250))

#Summerize data
data_summary <- summary(data)
data_summary

#Histogram plots
plot(clmfires$marks$cause)
hist(clmfires$marks$burnt.area, breaks = 100)
hist(log(clmfires$marks$burnt.area), breaks = 100)

# Plotting density of fires in CLM
plot(density(data,20), main = "Density plot of fires in CLM region - high focus")
plot(density(data,50), main = "Density plot of fires in CLM region - low focus")

#contour plot 
contour(density(data,8), axes = F, main = "Contour plot of CLM fires")

#Quadrant counting 7 by 7
quadrat_data <- quadratcount(data, nx = 7, ny = 7)
quadrat_data
plot(quadrat_data, main = "quadrat counting the data")

#kest
kest_data <- Kest(data)
kest_data
plot(kest_data)

#envelope
E <- envelope(data[selected_data],Kest,nsim=39)
plot(E)

#Strauss -> not working
#clm_region <- owin(c(8.247, 385.345), c(24.220, 377.177))
#clm_subregion <- owin(c(170, 350), c(70, 250))
#fit <- ppp(clmfires$x, clmfires$y, clm_subregion, ~1, Strauss(10))   # Fit a Strauss process to the data
#fit
#plot(simulate(fit))  # Plot the simulated data with the process
#plot(envelope(fit,Kest,nsim=39))


#############################    INTENSITY ANALYSIS     #########################

lamb <- summary(data)$intensity
cat( "Point intensity is equal to:", lamb )

#Kernel Plots of section 2
den_20 <- density(data,kernel = "epanechnikov",sigma=20)  # The resulting object is a pixel image. This class has methods for print, summary, plot, contour and persp
den_50 <- density(data,kernel = "epanechnikov",sigma=50)  # The resulting object is a pixel image. This class has methods for print, summary, plot, contour and persp
plot(den_20)
plot(den_50)
plot(data, add = TRUE, cols = "green", markscale = 0.3)
summary(den_20)
contour(den_20)
# 3D kernel representation of section 2
persp(den_20, col = "red", main = "3D representation of fire intensity")
boxplot(den_20$v)

# MAPS FOR SUBREGION
den_subregion_20 <- density(clmfires[clm_subregion],kernel = "epanechnikov",sigma=20)  
den_subregion_50 <- density(clmfires[clm_subregion],kernel = "epanechnikov",sigma=50)  
plot(den_subregion_20,main = "CLM subregion kernel density")
plot(clmfires[clm_subregion], add = TRUE, cols = "white", markscale = 0.3, pch = "*")

#3D representation of data
persp(den_subregion_20, col = "red", main = "3D representation of kernel density")

quadrat_data_subregion <- quadratcount(data[clm_subregion], nx = 7, ny = 7)
plot(quadrat_data_subregion, main = "quadrat count _ subregion")

# Map according to elavation
Z<- clmfires.extra$clmcov100$elevation
plot(Z[clm_subregion], main ="CLM region elavation")
b <- quantile(Z[clm_subregion],probs=(0:3)/3)
Zcut <- cut(Z[clm_subregion],breaks=b,labels=1:3)
V <- tess(image=Zcut)
plot(V[clm_subregion], main = "CLM region elavation - 3 categories")
plot(clmfires[clm_subregion], which.marks = "cause",add=TRUE,pch="+", cols = "green")   
qb <- quadratcount(clmfires[clm_subregion],tess=V) 
plot(qb, main = "Number of fire observation in each elavation category - subregion")

#Selection of the sub region for CLM to have a full fit of data
clm_region <- owin(c(8.247, 385.345), c(24.220, 377.177))
clm_subregion <- owin(c(170, 350), c(70, 250))
plot(clmfires[clm_subregion], pch = "+", which.marks = "cause", axes= TRUE)

#Intensity and covariates analysis for CLM 
plot(rhohat(ppp(clmfires$x, clmfires$y, clm_region), clmfires.extra$clmcov100$elevation), main ="intensity of fires based on elevation")
plot(rhohat(ppp(clmfires$x, clmfires$y, clm_region), clmfires.extra$clmcov100$slope), main ="intensity of fires based on slope")
plot(rhohat(ppp(clmfires$x, clmfires$y, clm_region), clmfires.extra$clmcov100$orientation), main ="intensity of fires based on orientation")
#Intensity and covariates analysis for CLM  - SUBREGION
plot(rhohat(ppp(clmfires$x, clmfires$y, clm_subregion), clmfires.extra$clmcov100$elevation), main ="intensity of fires based on elevation - subregion")
plot(rhohat(ppp(clmfires$x, clmfires$y, clm_subregion), clmfires.extra$clmcov100$slope), main ="intensity of fires based on slope - subregion")
plot(rhohat(ppp(clmfires$x, clmfires$y, clm_subregion), clmfires.extra$clmcov100$orientation), main ="intensity of fires based on orientation - subregion")
#iNTENSITY ANALYSIS BASED ON DISTANCE OF FIRES
plot(distmap(clmfires))
plot(rhohat(ppp(clmfires$x, clmfires$y, clm_region), distmap(clmfires)), main ="intensity of fires based on distance")

q_test7 <- quadrat.test(clmfires[clm_subregion], nx = 7, ny = 7)
q_test7
q_test30 <- quadrat.test(clmfires[clm_subregion], nx = 30, ny = 30)
q_test30
q_test15 <- quadrat.test(clmfires[clm_subregion], nx = 15, ny = 15)
q_test15

plot(clm_subregion)
plot(q_test7, add= TRUE)

## Poisson MLE TEST FOR Complete spatial randomness in CLM subregion
#Poisson model fitting
fit_1 <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~1)  
fit_1
fit_2 <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~x+y)  
fit_2
lam_2 <- predict(fit_2)
plot(lam_2, main = "fitted model of the Possion MLE - linear")
plot(clmfires[clm_subregion], which.marks = "cause", pch = "+", cols="white", add = TRUE)
#plot(effectfun(fit_2), main = "effect function of")
fit_3 <-ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~polynom(x,y,2))
fit_3
lam_3 <- predict(fit_3)
plot(lam_3, main = "fitted model of the Possion MLE - polynomial power 2")
plot(clmfires[clm_subregion], which.marks = "cause", pch = "+", cols="white", add = TRUE)

#polynomial model power 4
fit_4 <-ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~polynom(x,y,4))
fit_4
lam_4 <- predict(fit_4)
plot(lam_4, main = "fitted model of the Possion MLE - polynomial power 4")
plot(clmfires[clm_subregion], which.marks = "cause", pch = "+", cols="white", add = TRUE)
plot(effectfun(fit_4,"x", y=0.5))
plot(effectfun(fit_4,"y", x=0.5))

# MODEL WITH DISTANCE AS COVARIATE
z_distance <- distmap(clmfires[clm_subregion])
fit_5 <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~ polynom(z_distance,2))
fit_5
lam_5 <- predict(fit_5)
plot(lam_5, main = "fitted values of the Possion MLE - Distance power 2 model")
plot(effectfun(fit_5), main = "intensity as function of distance")

#ANOVA TEST
anova(fit_1, fit_3, test = "Chi")

# model with intensity function of elevation
elevation <- clmfires.extra$clmcov100$elevation
fit_elevation <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion), ~elevation, covariates = list(elevation = elevation)) 
fit_elevation
lam_6 <- predict(fit_elevation)
plot(lam_6)

#########################   Goodness of Fit ##############################
#quadrat MODEL FOR SUB REGION - LINEAR
fit_2 <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~x+y)  
fit_2
quadrat_fit_2 <- quadrat.test(fit_2, nx = 7, ny = 7)
quadrat_fit_2
plot(clmfires[clm_subregion], which.marks = "cause", pch = "*", main = "first order")
plot(quadrat_fit_2, add = TRUE, col = "red")
#quadrat model for sub region - polynomial power 4
fit_4 <-ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~polynom(x,y,4))
fit_4
quadrat_fit_4 <- quadrat.test(fit_4, nx = 7, ny = 7)
quadrat_fit_4
plot(clmfires[clm_subregion], which.marks = "cause", pch = "*", main = "polynomial power 4 model")
plot(quadrat_fit_4, add = TRUE, col = "red")

#MODEL POLYNOMIAL + DISTANCE + ELEVATION
fit_5 <-ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~polynom(x,y,4) + polynom(z_distance,2)+elevation)
fit_5
quadrat_fit_5 <- quadrat.test(fit_5, nx = 7, ny = 7)
quadrat_fit_5
plot(clmfires[clm_subregion], which.marks = "cause", pch = "*", main = "polynomial power 4 + distance + elevation model")
plot(quadrat_fit_5, add = TRUE, col = "red")


#Residual diagnostics of the linear model and polynomial model 
diagnose.ppm(fit_subregion_2)
diagnose.ppm(fit_subregion_4)
diagnose.ppm(fit_2, which = "smooth", main = "Smoother residuals for first order model")
diagnose.ppm(fit_4, which = "smooth", main = "Smoother residuals for power 4 polynomial model")
diagnose.ppm(fit_5, which = "smooth", main = "Smoother residuals for power 4 polynomial model + distance + elevation")

elevation <- clmfires.extra$clmcov100$elevation
orientation <- clmfires.extra$clmcov100$orientation
slope <- clmfires.extra$clmcov100$slope
z_distance <- distmap(clmfires[clm_subregion])
lurking(fit_4, elevation, type = "raw", main = "lurking plot elevation - only polynomial")
lurking(fit_4, slope, type = "raw", )
lurking(fit_4, orientation, type = "raw")
lurking(fit_4, z_distance, type = "raw", main = "lurking plot distance - only polynomial")

lurking(fit_5, elevation, type = "raw", main = "lurking plot elevation - with covariates")
lurking(fit_5, slope, type = "raw", )
lurking(fit_5, orientation, type = "raw")
lurking(fit_5, z_distance, type = "raw", main = "lurking plot distance - with covariates")


fit_5 <-ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~polynom(x,y,4) + z_distance)
diagnose.ppm(fit_5, which = "smooth", main = "Smoother residuals for power 4 polynomial model + distance")
lurking(fit_4, z_distance, type = "raw", main = "lurking variable plot for distance for polynomial distance model")

################################# POINT PATTERN DEPENDENCE ANALYSIS ################
#MORISHITA PLOT
miplot(ppp(clmfires$x, clmfires$y, clm_subregion), main = "Morishita plot of point dependence in CLM subregion")
# FRYPLOT
clm_subregion1 <- owin(c(200, 300), c(150, 250))
fryplot(ppp(clmfires$x, clmfires$y, clm_subregion1), axes = TRUE, main = "fryplot for subregion 1 - CLM")

#Pair distances
d_1 <- nndist(ppp(clmfires$x, clmfires$y, clm_subregion))

d_2 <- pairdist(ppp(clmfires$x, clmfires$y, clm_subregion))
d_3 <- pairdist(ppp(clmfires$x, clmfires$y, clm_subregion), periodic=TRUE)
d_4 <- pairdist(ppp(clmfires$x, clmfires$y, clm_subregion), squared=TRUE)

#POINT DEPENDENCE ANALYSIS USING K FUNCTION
#Kest
Gc <- Kest(ppp(clmfires$x, clmfires$y, clm_subregion))
Gc
par(pty = "s")
plot(Gc)
#L (R)
L <- Lest(ppp(clmfires$x, clmfires$y, clm_subregion))
plot(L, main = "L function")
# ENVELOPE CRITICAL THRESHOLDS
E <- envelope(ppp(clmfires$x, clmfires$y, clm_subregion), Kest, nsim = 39, rank = 1)
E
plot(E, main = "pointwise envelope")
E <- envelope(ppp(clmfires$x, clmfires$y, clm_subregion), Lest, nsim = 39, rank = 1, global = TRUE)  # This stabilizes the variance
E
plot(E, main = "global envelopes of L(r)")
v <- varblock(ppp(clmfires$x, clmfires$y, clm_subregion), Kest, nx=4, ny=2) # CI for the K function
plot(v)

#### NON POISSON MODELS
fit_10 <- kppm(ppp(clmfires$x, clmfires$y, clm_subregion),~1,"MatClust")
fit_10
plot(fit_10)
E <- envelope(fit_10,Kest,nsim= 59, global = TRUE, correction = "border")
plot(E)
#Covariates analysis
fit_11 <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion),~elevation, covariates = clmfires.extra)
lam_11 <- predict(fit_11,locations=clmfires[clm_subregion])
Ki <- Kinhom(clmfires[clm_subregion],lam_11)
plot(Ki,main = "Inhomogeneous K function")

######################### mark analysis ########################
summary(clmfires)

plot(data, which.marks="cause", cols=2:5, cex=1,
     main="Castilla-La Mancha forest fires - by cause", axes = TRUE)

plot(density(split(clmfires)), main = "density by cause categorical mark")

plot(clmfires$marks$cause)
hist(clmfires$marks$burnt.area, breaks = 100)
hist(log(clmfires$marks$burnt.area), breaks = 100)

#model with marks as covariates
fit_mark_1 <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion, marks = clmfires$marks$cause),~marks + polynom(x,y,2))
fit_mark_1
lam_mark_1 <- predict(fit_mark_1)
plot(lam_mark_1, main = "fitted model with marks_1")

# model with marks + polynomial power 2
fit_mark_2 <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion, marks = clmfires$marks$cause),~marks * polynom(x,y,2))
fit_mark_2
lam_mark_2 <- predict(fit_mark_2)
plot(lam_mark_2, main = "fitted model_Polynomial power 2")

# model with marks + polynomial power 4
fit_mark_3 <- ppm(ppp(clmfires$x, clmfires$y, clm_subregion, marks = clmfires$marks$cause),~marks * polynom(x,y,4))
fit_mark_3
lam_mark_3 <- predict(fit_mark_3)
plot(lam_mark_3, main = "fitted model_Polynomial power 4")

#ANOVA for fit
anova(fit_mark_2,fit_mark_3,test="Chi")


