# The following is the script for the main analysis described in the manuscript:
# 

# Dylan Jon Keel dkeel@res.us
# Last updated 11-27-2024
# 


#library required packages

# library(doSNOW)
# library(doParallel)
# library(doMPI)
library(tidyverse)
library(vegan)
library(tictoc)
library(fields)
library(ggplot2)
library(sp)
library(sf)
library(lattice)
library(INLA)
library(plotly)
library(geometry)
library(viridis)
#library(devtools) #might need devtools to install some packaages
library(tictoc)
library(kableExtra)

library(rgdal)

#sometimes this happens...
#install.packages("gstat",force=TRUE)
library(gstat)
library(remotes)
#library(INLAOutputs)


#library(inlamesh3d)
library(inlatools)
library(INLAutils)
library(ggregplot)

#getwd()

###############

# source useful scripts
source("HighstatLibV13.R")
source("INLA_plotting_functions.R")

flow.data=read.csv("flow.cor.data.klamath.meta.2023.csv")

meta_table<-read.csv("klamath.meta.covariates.2023.csv",check.names=FALSE)

meta_table <- meta_table[meta_table$Year == 2023, ]


# Step 2: Calculate Shannon's Diversity Index and Species Richness
diversity_data.flow <- flow.data%>%
  rowwise() %>%
  mutate(
    ShannonsDI = diversity(c_across(Chinook.Salmon:Speckled.Dace.Dace.spp.), index = "shannon"),  # Shannon's Diversity Index
    Richness = sum(c_across(Chinook.Salmon:Speckled.Dace.Dace.spp.) > 0)  # Species Richness (non-zero species count)
  ) %>%
  select(Site.Name, ShannonsDI, Richness)  # Keep only necessary columns

# Rename empty columns in meta_table
names(meta_table)[names(meta_table) == ""] <- paste0("Unnamed_", seq_along(names(meta_table))[names(meta_table) == ""])
meta_table$Unnamed_1=NULL

flow.data <- dplyr::left_join(diversity_data.flow, meta_table, by = "Site.Name")


par(mar=c(2,2,2,2))
# hist(total.data$ShannonsDI,breaks=50)
# hist(total.data$Richness,breaks=10)
# 
# hist(vol.data$ShannonsDI,breaks=50)
# hist(vol.data$Richness,breaks=10)

hist(flow.data$ShannonsDI,breaks=50)
hist(flow.data$Richness,breaks=10)

## check for multicolinearity
colnames(flow.data)

temp_data <- data.frame(lapply(flow.data, function(x) {
  if (is.character(x) || is.factor(x)) {
    as.numeric(factor(x))  # Convert character or factor to numeric based on unique values
  } else {
    x  # Leave numeric columns unchanged
  }
}))
colnames(temp_data)
MyVar.temp<-c("Air.Temperature"   ,     "Water.Temperature"    ,  "Dissolved.Oxygen"   ,   
              "Specific.Conductance" ,  "pH"    ,       # "Q.CMS",              
              "Habitat"        ,        "Size"       ,            "Reference"     ,         "Control"   )


MyVar.cont <- c("Air.Temperature"  ,    "Water.Temperature" ,   "Dissolved.Oxygen"  ,   "Specific.Conductance", "pH"  # , "Total_Volume"    ,           
                #"Q.CMS"   
                )

Mypairs(temp_data[,MyVar.temp])

# procede when finding none


#scale the covariates

flow.data.scale <- dplyr::select(flow.data, all_of(MyVar.cont)) 
names(flow.data.scale) <- paste(MyVar.cont, '.z', sep='')

scale.par = do.call(data.frame, 
                    list(mean = apply(flow.data.scale, 2, mean),
                         sd = apply(flow.data.scale, 2, sd),
                         min = apply(flow.data.scale, 2, min),
                         max = apply(flow.data.scale, 2, max)))

flow.data.scale <- as.data.frame(scale(flow.data.scale))
flow.data <- cbind(flow.data, flow.data.scale)
#total.data <-cbind(total.data,flow.data.scale)
#flow.data<-cbind(flow.data,flow.data.scale)

#now check all the covariates for multicolinearity
colnames(flow.data)
#
MyVar <- c("Air.Temperature.z"    ,  # "Total_Volume.z", #"Q.CMS.z",
           "Water.Temperature.z"    ,  "Dissolved.Oxygen.z"  ,     "Specific.Conductance.z"  ,
           "pH.z"             ,           
           "Habitat"    ,            "Size"            ,       "Reference"       ,      
           "Control"    )


#Mypairs(flow.data[,MyVar])
# plot them against the response variable


# MyMultipanel.ggp2(Z = vol.data, 
#                   varx = MyVar, 
#                   vary = "ShannonsDI", 
#                   ylab = "ShannonsDI",
#                   addSmoother = T,
#                   addRegressionLine = F,
#                   addHorizontalLine = FALSE)
# 
# 
# MyMultipanel.ggp2(Z = vol.data, 
#                   varx = MyVar, 
#                   vary = "Richness", 
#                   ylab = "Richness",
#                   addSmoother = T,
#                   addRegressionLine = F,
#                   addHorizontalLine = FALSE)
# 
# 
# MyMultipanel.ggp2(Z = total.data, 
#                   varx = MyVar, 
#                   vary = "ShannonsDI", 
#                   ylab = "ShannonsDI",
#                   addSmoother = T,
#                   addRegressionLine = F,
#                   addHorizontalLine = FALSE)
# 
# 
# MyMultipanel.ggp2(Z = total.data, 
#                   varx = MyVar, 
#                   vary = "Richness", 
#                   ylab = "Richness",
#                   addSmoother = T,
#                   addRegressionLine = F,
#                   addHorizontalLine = FALSE)


MyMultipanel.ggp2(Z = flow.data, 
                  varx = MyVar, 
                  vary = "ShannonsDI", 
                  ylab = "ShannonsDI",
                  addSmoother = T,
                  addRegressionLine = F,
                  addHorizontalLine = FALSE)


MyMultipanel.ggp2(Z = flow.data, 
                  varx = MyVar, 
                  vary = "Richness", 
                  ylab = "Richness",
                  addSmoother = T,
                  addRegressionLine = F,
                  addHorizontalLine = FALSE)



nonzero_values <- flow.data$ShannonsDI[flow.data$ShannonsDI > 0]

# Plot histogram
hist(flow.data$ShannonsDI, breaks = 10, col = "red", border = "black", 
     main = "Histogram with Zero Values", 
     xlab = "Value", ylab = "Frequency")
hist(nonzero_values, breaks = 10, col = "blue", border = "black", 
     main = "Histogram of Non-Zero Values", 
     xlab = "Value", ylab = "Frequency")




#make sure data is in the right format
####FIX BEFORE RUNNING OTHER SITES####
flow.data$Habitat=as.factor(flow.data$Habitat)

#flow.data$Habitat<-recode_factor(All_Data$YearR, "2021"=2, "2022"=3,"2023"=4,"2024"=5)
#YearR=as.numeric(All_Data$YearR)
#All_Data=subset(All_Data,select=-YearR)

flow.data$Size=as.factor(flow.data$Size)
flow.data$Reference=as.factor(flow.data$Reference)
flow.data$Control=as.factor(flow.data$Control)
flow.data$Stream=gsub("[0-9]", "", flow.data$Site.Name)
flow.data$Stream=gsub("Klamath_Ranch", "Klamath", flow.data$Stream)

# Assign the new specific stream names
flow.data$Stream[flow.data$Site.Name %in% c( "Klamath34", "Klamath36")] <- "JCBoyle"
flow.data$Stream[flow.data$Site.Name %in% c("Klamath12", "Klamath10", "Beaver1")] <- "Copco"
flow.data$Stream[flow.data$Site.Name %in% c("Klamath6", "Klamath4", "Jenny1", "Scotch1", "Klamath2")] <- "IronGate"

flow.data$Stream=as.factor(flow.data$Stream)


flow.data$ShannonsDI <- flow.data$ShannonsDI + 1e-6

ShannonsDI=flow.data$ShannonsDI

# 
# #No quads yet


plot(flow.data$ShannonsDI~flow.data$Habitat)
# NumZones.z=as.numeric(All_Data$NumZones.z)
# Number.of.Days.Inundated.z=as.numeric(All_Data$Number.of.Days.Inundated.z)
# MaxD.cm.z=as.numeric(All_Data$MaxD.cm.z)
# Area.m.z=as.numeric(All_Data$Area.m.z)
# NeighborCount.z=as.numeric(All_Data$NeighborCount.z)
# PerimeterM.z=as.numeric(All_Data$PerimeterM.z)
# PerimAreaRatio.z=as.numeric(All_Data$PerimAreaRatio.z)
# MaxD.cm.z2=as.numeric(All_Data$MaxD.cm.z)*as.numeric(All_Data$MaxD.cm.z)
# Number.of.Days.Inundated.z2=as.numeric(All_Data$Number.of.Days.Inundated.z)*as.numeric(All_Data$Number.of.Days.Inundated.z)


#make base model forms
base.form.pcp <- ShannonsDI ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z+ f(Stream,model="iid")

#pcprior <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))

control.fixed <- list(
  mean = 0,  # Mean for all fixed effects
  prec = 1e-6  # Very low precision (high variance), effectively flat
)

# Set uninformative prior for the precision of the Gamma distribution
# Inverse-gamma prior for precision
control.family <- list(
  hyper = list(
    prec = list(
      prior = "loggamma",  # Use a log-gamma distribution for precision
      param = c(1, 0.00001),  # Uninformative prior: low shape (1) and high variance (small rate)
      initial = log(1),  # Set initial value for the precision
      fixed = FALSE  # Let INLA estimate the precision
    )
  )
)



######here we fit some models without spatial effects
####normal non-spatial with distance from shore
# Prepare data: Add a column to categorize based on ShannonsDI

# Update Habitat labels

flow.data <- flow.data %>%
  mutate(Habitat_Label = ifelse(as.character(Habitat) == "Reservoir", "Reservoir", "Stream"),
         Category = ifelse(ShannonsDI > 0.05, "Above 0.05", "0.05 or Below"))
#names(flow.data)

# Ensure Stream factor order aligns with Site.Name order
flow.data$Stream <- factor(flow.data$Stream, levels = unique(flow.data$Stream[order(flow.data$Site.Name)]))

# Generate the plot with Site.Name as the x-axis labels
ggplot(flow.data, aes(x = reorder(Site.Name, ShannonsDI), y = ShannonsDI, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Habitat_Label, scales = "free_x") +  # Facet by Habitat_Label (Reservoir/Stream)
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels for readability
    strip.text = element_text(size = 12),  # Enhance facet label size
    panel.grid.major.x = element_blank()  # Remove x-axis major grid lines for a cleaner look
  ) +
  labs(
    title = "Shannon's Diversity Index (Native Fishes) by Site and Habitat Type",
    x = "Site Name",
    y = "Shannon's DI"
  ) +
  scale_fill_manual(values = c("Above 0.05" = "#00577b", "0.05 or Below" = "#97a641"))


set.seed(1993)



flow.data <- flow.data %>%
  mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))

# Step 1: Fit the zero-inflation model
zero_inflated <- flow.data  # replace flow.data with your actual data frame
zero_inflated$has_non_zero <- as.numeric(flow.data$ShannonsDI > 0.05)

# Logistic regression for zero-inflation
binary_formula <- has_non_zero ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z+  f(Stream,model="iid")

new_data=zero_inflated

binary_fit <- inla(binary_formula, 
                   family = "binomial", 
                   data = zero_inflated, 
                   control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE,return.marginals.predictor=TRUE),
                   control.predictor = list(link = 1, compute = TRUE))

summary(binary_fit)


# Extract fitted values (mean posterior predictions)
fitted_values <-binary_fit$summary.fitted.values$mean

# Combine fitted values with original data
zero_inflated$fitted <- fitted_values

library(pROC)

# Fit your INLA model (assuming you already have a model)
# model <- inla(has_non_zero ~ Habitat, family = "binomial", data = zero_inflated)

# Extract the fitted probabilities from your INLA model
fitted_probs <- binary_fit$summary.fitted.values$mean

# Calculate the ROC curve and AUC
roc_curve <- roc(zero_inflated$has_non_zero, fitted_probs)

# Display AUC value
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))

# Optionally plot the ROC curve
par(mar=c(3,3,3,3))
par(mfrow = c(1,1))
plot(roc_curve, col = "black", lwd = 2, main = "ROC Curve for INLA Binomial Model")
abline(a = 0, b = 1, col = "grey", lty = 2)  # Add a diagonal reference line



# Step 1: Generate posterior predictive replications
n_replications <- 10000  # Number of replicated datasets
predictions <- inla.posterior.sample(n_replications, binary_fit)

# Step 2: Compute discrepancy measure (e.g., mean residuals for simplicity)
obs_discrepancy <- abs(zero_inflated$has_non_zero - binary_fit$summary.fitted.values$mean)

rep_discrepancy <- sapply(predictions, function(rep) {
  replicated_data <- sapply(rep$latent, mean)  # Average over replicated values
  abs(zero_inflated$has_non_zero - replicated_data)  # Calculate discrepancy
})

# Step 3: Calculate Bayesian p-value
bayesian_p_value <- mean(rep_discrepancy > obs_discrepancy)

bayesian_p_value


summary(binary_fit)


exp(7.585)
exp(31.511)


1- exp(-7.198)
1 -exp(-.873)

plot(flow.data$ShannonsDI~flow.data$Size)
######

# After accounting for the random effect of waterbody...

# Streams are at least 1968 times (1968 - 4.84 x 10^13) more likely to have native 
# fish diversity that is not extremely low (Shannon's Diversity Index < 0.05).

# The Scott River reference sites have between 0 and 42% less native fish diversity than
# the Upper Klamath River and tributaries.

# Our model fit the data well, with an AUC of 0.928 indicating a 92.8% true positive rate. 
# Additionally, our bayesian p-value was between 0.1 and 0.9 (0.8), indicating 
# goodness of fit.

# Because our posterior probability distributions spanned zero for Control/Impact,
#  Dissolved Oxygen Concentration, and Specific Conductance, there
# is no evidence that these variables described a meaningful amount of variance in
#native fish species richness.


#### we assessed the diversity model for spatial autocorrelation following the same steps as the richness model above
#### we did not find evidence of residual spatial autocorrelation




#average distance between sites


# check out the model
par(mar=c(3,3,3,3))

observed<-zero_inflated$has_non_zero

D<-INLAutils::plot_inla_residuals(binary_fit, observed)



# Let's make a simplified world to view our data in 2D to check for
# spatial-autocorrelation. Let's assume that the river is a straight
# line and that sinuosity is not important.


#average distance between sites

xy <- with(flow.data, data.frame(Site.Name, X.Coordinate.x, Y.Coordinate.x))

# Convert the lat/long dataframe into a spatial object
latlong_sf <- st_as_sf(xy, coords = c("X.Coordinate.x", "Y.Coordinate.x"), crs = 4326)  # EPSG:4326 is WGS84

# Transform the lat/long coordinates to UTM
# sf automatically selects the appropriate UTM zone based on your data location
utm_sf <- st_transform(latlong_sf, crs = 32610)  # Example: UTM zone 10N, adjust accordingly

# View the results
#print(utm_sf)

# Extract the UTM coordinates
utm_coords <- st_coordinates(utm_sf)

# Replace the lat/long in the xy object with UTM coordinates
xy$UTM_Easting <- utm_coords[,1]  # UTM Easting (X coordinate)
xy$UTM_Northing <- utm_coords[,2] # UTM Northing (Y coordinate)

# Remove old Latitude and Longitude columns if you want
xy$Latitude <- NULL
xy$Longitude <- NULL

# View the updated dataframe
#print(xy)


#nonspat w/out velocity adj
Pi1 <- binary_fit$summary.fitted.values[, "mean"]



D1  <- (observed - Pi1) / sqrt(Pi1)

summary(D1)

#D1
# e <- Cs.Copies-Pi1

# dev.res <- sign(e)*sqrt(-2*(Cs.Copies*log(Pi1) + (1 - Cs.Copies)*log(1 - Pi1)))
# #D1
#Pi1
### shrink x 2 orders of magnitude
MyData <- data.frame(D1  = as.numeric(D1),
                     X = as.numeric(xy$UTM_Easting),
                     Y = as.numeric(xy$UTM_Northing))



###### here we look for spatial patterns in the residuals
MyData$MySize <- 2 * abs(MyData$D1) / max(MyData$D1)
MyData$MyCol <- ifelse(MyData$D1> 0, 1, 2)
lattice::xyplot(Y ~ X,
                data = MyData,
                cex = MyData$MySize,
                col = MyData$MyCol,
                pch = 1)

## we can see that sites 3 and 6 have clear spatial patterns
# in residuals
#View(MyData)


V1a=NULL

## now we make a semivariogram

coordinates(MyData) <- c("X","Y")
V1a <- gstat::variogram(D1 ~ X + Y, 
                        data = MyData, 
                        cutoff = 6000,
                        cressie = TRUE)

hist(V1a$dist)

p <- ggplot()
p <- p + geom_point(data = V1a,
                    aes(x = dist,
                        y = gamma))
p <- p + geom_smooth(data = V1a,
                     span = 0.9,
                     se = FALSE,
                     aes(x = dist,
                         y = gamma))
p <- p + xlab("Distance(m)") + ylab("Sample variogram")
p1 <- p + theme(text = element_text(size=15))+ylim(0,.25)
p1 

#with a sill of ~0.625 and a nugget of ~.25 we find:

Sill.Nug.R=.075/.125
1-Sill.Nug.R

# Your plot
p <- ggplot() +
  geom_point(data = V1a, aes(x = dist, y = gamma)) +
  geom_smooth(data = V1a, span = 0.9, se = FALSE, aes(x = dist, y = gamma)) +
  xlab("Distance (m)") +
  ylab("Sample variogram") +
  theme(text = element_text(size = 15)) +
  ylim(0, .25) +
  # Add horizontal line for the Sill
  geom_hline(yintercept = 0.125, color = "green", linetype = "dashed") +
  annotate("text", x = max(V1a$dist, na.rm = TRUE) * 0.9, y = 0.125, label = "Sill", color = "green", hjust = 0) +
  # Add horizontal line for the Nugget
  geom_hline(yintercept = 0.075, color = "orange", linetype = "dashed") +
  annotate("text", x = max(V1a$dist, na.rm = TRUE) * 0.9, y = 0.075, label = "Nugget", color = "orange", hjust = 0) +
  # Add vertical line for the Range
  geom_vline(xintercept = 2500, color = "red", linetype = "dashed") +
  annotate("text", x = 2500, y = max(V1a$gamma, na.rm = TRUE) * 0.95, label = "Range", color = "red", angle = 90, hjust = -0.5)

# Display the plot
p


#we see that spatial autocorrelation increases from 0 to 2500 meters
# in our data

Loc<-NULL
Loc <- cbind(xy$UTM_Easting,xy$UTM_Northing)
#Loc
#what are the distances between the points?
D <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between samples (m)",
     ylab = "Frequency")
#View(Loc)


RangeGuess <-2500


# There is still a spatial correlation in the residuals of spat.inla.  The first thing that you should try is to
# use a mesh with more triangles. Use a smaller MaxEdge and also a smaller
# cutoff. By doing that we allow for smaller-scale correlation.
# Right now we have:

require(splancs)



sum(is.na(Loc))
Loc <- na.omit(Loc)

Hull <- inla.nonconvex.hull(Loc, convex = -0.2)
MaxEdge  <- RangeGuess/10
mesh2d     <- inla.mesh.2d(boundary = Hull,
                           max.edge = c(1,40) * MaxEdge,
                           cutoff = MaxEdge/20 ,
                           offset = c(5, 100) ,
                           max.n=20000)  #control max.n to be <3000
# 
mesh2d$n
# 

# #back to 2D
par(mfrow = c(1,1), mar=c(0,0,0,0))
plot(mesh2d, asp=1, main = "")
points(Loc, col = 2, pch = 1, cex = 1.5)



####fix na values in covariates for spatial model
zero_inflated$Specific.Conductance.z[is.na(zero_inflated$Specific.Conductance.z)] <- mean(zero_inflated$Specific.Conductance.z, na.rm = TRUE)



# Define projector matrices for the mesh.
A <- inla.spde.make.A(mesh2d, loc = Loc)


spde <- inla.spde2.pcmatern(mesh2d, 
                            prior.range = c(RangeGuess, 0.5), 
                            prior.sigma = c(1, .05))

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)
#w.index


stackform<-as.formula(binary_formula)


modrando<-c("Stream")

terms <- base.form.pcp[3]
terms
terms <- gsub(" ", "", as.character(terms))
terms <- unlist(strsplit(terms, '\\+'))
terms <- terms[-grep('iid', terms)]
if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
terms <- c(terms, modrando)
terms

## make the formula for the stack
stackform <- formula(paste('~', paste(terms, collapse='+')))
stackform

Xm <- model.matrix(stackform, data = zero_inflated,na.action=na.pass)


Xm <- as.data.frame(Xm)
N <- nrow(zero_inflated)
#w.index


head(zero_inflated)


has_non_zero=zero_inflated$has_non_zero

StackFit <- inla.stack(
  remove.unused = T,
  tag = "Fit",
  data = list(
    has_non_zero = has_non_zero
  ),
  A = list(1, 1, A),
  effects = list(
    Intercept = rep(1, N),
    Xm = Xm[,-1],  # Use covariates without the intercept (includes dummy variables)
    w = w.index    # Ensure w.index is correctly aligned with random effects
  )
)
# 
# StackFit <- inla.stack(
#   remove.unused = TRUE,
#   tag = "Fit",
#   data = list(Shannon = vol.data$ShannonsDI),  # Ensure this variable exists in vol.data
#   A = list(1,1, A),  # Double-check that A is correctly sized for vol.data
#   effects = list(
#     Intercept = rep(1, nrow(vol.data)),
#     Xm = model.matrix(~ Habitat + Size + Control + Reference + Dissolved.Oxygen.z + Specific.Conductance.z - 1, data = vol.data), 
#     w = w.index
#   )
# )
# 
# # Use only relevant terms from Xm and provide the response separately
# StackFit <- inla.stack(
#  # remove.unused = TRUE,
#   tag = "Fit",
#   data = list(Shannon = ShannonDI),  # Response variable as separate vector
#   A = list(1, 1, A),  # A matrix for projections (make sure this aligns with w.index)
#   effects = list(
#     Intercept = rep(1, nrow(Xm)),  # Intercept term
#     Xm = Xm[,-1],  # Use the covariates from Xm without the intercept
#     w = w.index    # Ensure w.index aligns correctly
#   )
# )

#inla.stack.data(StackFit)

#Fit the global model
# 
# base.form.pcp <- ShannonsDI ~ Habitat + Size + Control+ Reference +
#   Dissolved.Oxygen.z + Specific.Conductance.z+  f(Stream,model="iid")

base.form.pcp <- has_non_zero ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z +#+  f(Stream,model="iid")
  f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)



tic()
spat.inla.diversity<- inla(formula = update(base.form.pcp, . ~ . +
                                            f(Stream,model="iid")),
                         family = "binomial",
                         data = inla.stack.data(StackFit),
                         control.compute = list(waic = TRUE, config=TRUE,cpo=TRUE, po=TRUE,return.marginals.predictor=TRUE),
                         # control.family = control.family,
                         # control.inla = list(strategy = "gaussian"), 
                         control.predictor = list(link = 1, compute = TRUE,
                                                  A = inla.stack.A(StackFit)))

toc()


summary(spat.inla.diversity)



# Extract fitted values (mean posterior predictions)
fitted_values <-spat.inla.diversity$summary.fitted.values$mean[1:42]

# Combine fitted values with original data
zero_inflated$fitted <- fitted_values

library(pROC)

# Fit your INLA model (assuming you already have a model)
# model <- inla(has_non_zero ~ Habitat, family = "binomial", data = zero_inflated)

# Extract the fitted probabilities from your INLA model
fitted_probs <- spat.inla.diversity$summary.fitted.values$mean[1:42]

# Calculate the ROC curve and AUC
roc_curve <- roc(zero_inflated$has_non_zero, fitted_probs)

# Display AUC value
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))

# Optionally plot the ROC curve
par(mar=c(3,3,3,3))
par(mfrow = c(1,1))
plot(roc_curve, col = "black", lwd = 2, main = "ROC Curve for INLA Binomial Model")
abline(a = 0, b = 1, col = "grey", lty = 2)  # Add a diagonal reference line



# Step 1: Generate posterior predictive replications
n_replications <- 1000  # Number of replicated datasets
predictions <- inla.posterior.sample(n_replications, spat.inla.diversity)

# Step 2: Compute discrepancy measure (e.g., mean residuals for simplicity)
obs_discrepancy <- abs(zero_inflated$has_non_zero - spat.inla.diversity$summary.fitted.values$mean)

rep_discrepancy <- sapply(predictions, function(rep) {
  replicated_data <- sapply(rep$latent, mean)  # Average over replicated values
  abs(zero_inflated$has_non_zero - replicated_data)  # Calculate discrepancy
})

# Step 3: Calculate Bayesian p-value
bayesian_p_value <- mean(rep_discrepancy > obs_discrepancy)

bayesian_p_value


summary(spat.inla.diversity)

#habitat
exp(12.6)
exp(53.6)

#size
1- exp(-41.1)
1 -exp(-7.4)

#control
1- exp(-19.4)
1 -exp(-0.18)

#reference
1- exp(-22.6)
1 -exp(-1.3)


# After accounting for the random effect of waterbody and spatial position...

# The log of the odds that streams have non-low native fish diversity is 
# ~13-54x greater than the log of the odds that reservoirs have non-low
# fish diversity that is extremely low (Shannon's Diversity Index < 0.05).

# The odds that mainstem sites have non-low native fish diversity is a near certainty
# to be lower than the odds that tributaries have non-low
# fish diversity that is extremely low (Shannon's Diversity Index < 0.05).

# Although the model finds that the odds that control and reference sites have 
# non-low native fish diversity that is a near certainty
# to be lower than the odds that impact sites have non-low
# fish diversity that is extremely low (Shannon's Diversity Index < 0.05), 
# low sample size of both makes these parameter estimates difficult to interpret.

# Our model fit the data well, with an AUC of 0.96 indicating a 96% true positive rate. 
# Additionally, our bayesian p-value was between 0.1 and 0.9 (0.64), indicating 
# goodness of fit.

# Because our posterior probability distributions spanned zero for Control/Impact,
#  Dissolved Oxygen Concentration, and Specific Conductance, there
# is no evidence that these variables described a meaningful amount of variance in
#native fish species richness.















#### Now we fit the richness models in the same way
# Logistic regression for zero-inflation
base.form.rich <- Richness ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z+  f(Stream,model="iid")


nonspat.rich <- inla(base.form.rich,
                     
                     control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE,return.marginals.predictor=TRUE),
                     family = "poisson", 
                     
                     control.predictor = list(link = 1, compute = TRUE),
                     data = flow.data)
summary(nonspat.rich)



nonspat.d2.nbn <- inla(base.form.rich,
                       
                       control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                       family = "nbinomial", 
                       
                       control.predictor = list(link = 1, compute = TRUE),
                       data = flow.data)
summary(nonspat.d2.nbn)

#we'll check to see if zero inflation is warranted

nonspat.d2.zip <- inla(base.form.rich,
                       
                       control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                       family = "zeroinflatedpoisson1", 
                       
                       control.predictor = list(link = 1, compute = TRUE),
                       data = flow.data)
summary(nonspat.d2.zip)



nonspat.d2.zinbn <- inla(base.form.rich,
                         
                         control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                         family = "zeroinflatednbinomial1", 
                         
                         control.predictor = list(link = 1, compute = TRUE),
                         data = flow.data)
summary(nonspat.d2.zinbn)


###########################Check for Overdispersion

#source("Modified Distribution Check inlatools.R")
#source("Modified Dispersion Check inlatools.R")



nonspat.d2_dc <- dispersion_check(nonspat.rich)
nonspat.d2.nbn_dc <- dispersion_check(nonspat.d2.nbn)
nonspat.d2.zip_dc <- dispersion_check(nonspat.d2.zip)
nonspat.d2.zinbn_dc <- dispersion_check(nonspat.d2.zinbn)



# Here's what we want to test our models for overdisperssion
plot(nonspat.d2_dc)

mean(nonspat.d2_dc$data)
mean(nonspat.d2_dc$model)
#the poisson model is clearly underdispersed we want the 
#  0.1 < P(D|data>D|model) < 0.9

mean(nonspat.d2.nbn_dc$model)

# that's low, but underdispersion can arise from an autocorrelation
# structure. Let's not rule out the NB distribution until we check
# for spatial autocorrelation


# check out the model
par(mar=c(3,3,3,3))

observed<-flow.data$Richness

D<-INLAutils::plot_inla_residuals(nonspat.rich, observed)


#P<-autoplot(nonspat.d2.zinbn,which=(c(1:5)))
#P

#summary(nonspat.d2.zinbn)
summary(nonspat.rich)


# Let's make a simplified world to view our data in 2D to check for
# spatial-autocorrelation. Let's assume that the river is a straight
# line and that sinuosity is not important.


#average distance between sites

xy <- with(flow.data, data.frame(Site.Name, X.Coordinate.x, Y.Coordinate.x))

# Convert the lat/long dataframe into a spatial object
latlong_sf <- st_as_sf(xy, coords = c("X.Coordinate.x", "Y.Coordinate.x"), crs = 4326)  # EPSG:4326 is WGS84

# Transform the lat/long coordinates to UTM
# sf automatically selects the appropriate UTM zone based on your data location
utm_sf <- st_transform(latlong_sf, crs = 32610)  # Example: UTM zone 10N, adjust accordingly

# View the results
#print(utm_sf)

# Extract the UTM coordinates
utm_coords <- st_coordinates(utm_sf)

# Replace the lat/long in the xy object with UTM coordinates
xy$UTM_Easting <- utm_coords[,1]  # UTM Easting (X coordinate)
xy$UTM_Northing <- utm_coords[,2] # UTM Northing (Y coordinate)

# Remove old Latitude and Longitude columns if you want
xy$Latitude <- NULL
xy$Longitude <- NULL

# View the updated dataframe
#print(xy)


#nonspat w/out velocity adj
Pi1 <- nonspat.rich$summary.fitted.values[, "mean"]



D1  <- (observed - Pi1) / sqrt(Pi1)

summary(D1)

#D1
# e <- Cs.Copies-Pi1

# dev.res <- sign(e)*sqrt(-2*(Cs.Copies*log(Pi1) + (1 - Cs.Copies)*log(1 - Pi1)))
# #D1
#Pi1
### shrink x 2 orders of magnitude
MyData <- data.frame(D1  = as.numeric(D1),
                     X = as.numeric(xy$UTM_Easting),
                     Y = as.numeric(xy$UTM_Northing))



###### here we look for spatial patterns in the residuals
MyData$MySize <- 2 * abs(MyData$D1) / max(MyData$D1)
MyData$MyCol <- ifelse(MyData$D1> 0, 1, 2)
lattice::xyplot(Y ~ X,
                data = MyData,
                cex = MyData$MySize,
                col = MyData$MyCol,
                pch = 1)

## we can see that sites 3 and 6 have clear spatial patterns
# in residuals
#View(MyData)


V1a=NULL

## now we make a semivariogram

coordinates(MyData) <- c("X","Y")
V1a <- gstat::variogram(D1 ~ X + Y, 
                        data = MyData, 
                        cutoff = 6000,
                        cressie = TRUE)

hist(V1a$dist)

p <- ggplot()
p <- p + geom_point(data = V1a,
                    aes(x = dist,
                        y = gamma))
p <- p + geom_smooth(data = V1a,
                     span = 0.9,
                     se = FALSE,
                     aes(x = dist,
                         y = gamma))
p <- p + xlab("Distance(m)") + ylab("Sample variogram")
p1 <- p + theme(text = element_text(size=15))+ylim(0,1.2)
p1 

#with a sill of ~0.625 and a nugget of ~.25 we find:

Sill.Nug.R=.25/.625
1-Sill.Nug.R

# Your plot
p <- ggplot() +
  geom_point(data = V1a, aes(x = dist, y = gamma)) +
  geom_smooth(data = V1a, span = 0.9, se = FALSE, aes(x = dist, y = gamma)) +
  xlab("Distance (m)") +
  ylab("Sample variogram") +
  theme(text = element_text(size = 15)) +
  ylim(0, 1.2) +
  # Add horizontal line for the Sill
  geom_hline(yintercept = 0.625, color = "green", linetype = "dashed") +
  annotate("text", x = max(V1a$dist, na.rm = TRUE) * 0.9, y = 0.625, label = "Sill", color = "green", hjust = 0) +
  # Add horizontal line for the Nugget
  geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed") +
  annotate("text", x = max(V1a$dist, na.rm = TRUE) * 0.9, y = 0.25, label = "Nugget", color = "orange", hjust = 0) +
  # Add vertical line for the Range
  geom_vline(xintercept = 2500, color = "red", linetype = "dashed") +
  annotate("text", x = 2500, y = max(V1a$gamma, na.rm = TRUE) * 0.95, label = "Range", color = "red", angle = 90, hjust = -0.5)

# Display the plot
p


#we see that spatial autocorrelation increases from 0 to 2500 meters
# in our data

Loc<-NULL
Loc <- cbind(xy$UTM_Easting,xy$UTM_Northing)
#Loc
#what are the distances between the points?
D <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between samples (m)",
     ylab = "Frequency")
#View(Loc)


RangeGuess <-2500


# There is still a spatial correlation in the residuals of spat.inla.  The first thing that you should try is to
# use a mesh with more triangles. Use a smaller MaxEdge and also a smaller
# cutoff. By doing that we allow for smaller-scale correlation.
# Right now we have:

require(splancs)



sum(is.na(Loc))
Loc <- na.omit(Loc)

Hull <- inla.nonconvex.hull(Loc, convex = -0.2)
MaxEdge  <- RangeGuess/10
mesh2d     <- inla.mesh.2d(boundary = Hull,
                           max.edge = c(1,40) * MaxEdge,
                           cutoff = MaxEdge/20 ,
                           offset = c(5, 100) ,
                           max.n=20000)  #control max.n to be <3000
# 
mesh2d$n
# 

# #back to 2D
par(mfrow = c(1,1), mar=c(0,0,0,0))
plot(mesh2d, asp=1, main = "")
points(Loc, col = 2, pch = 1, cex = 1.5)



####fix na values in covariates for spatial model
flow.data$Specific.Conductance.z[is.na(flow.data$Specific.Conductance.z)] <- mean(flow.data$Specific.Conductance.z, na.rm = TRUE)



# Define projector matrices for the mesh.
A <- inla.spde.make.A(mesh2d, loc = Loc)


spde <- inla.spde2.pcmatern(mesh2d, 
                            prior.range = c(RangeGuess, 0.5), 
                            prior.sigma = c(1, .05))

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)
#w.index

base.form.pcp <- Richness ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z+ f(Stream,model="iid")

stackform<-as.formula(base.form.pcp)


modrando<-c("Stream")

terms <- base.form.pcp[3]
terms
terms <- gsub(" ", "", as.character(terms))
terms <- unlist(strsplit(terms, '\\+'))
terms <- terms[-grep('iid', terms)]
if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
terms <- c(terms, modrando)
terms

## make the formula for the stack
stackform <- formula(paste('~', paste(terms, collapse='+')))
stackform

Xm <- model.matrix(stackform, data = flow.data,na.action=na.pass)


Xm <- as.data.frame(Xm)
N <- nrow(flow.data)
#w.index

Richness=flow.data$Richness

StackFit <- inla.stack(
  remove.unused = T,
  tag = "Fit",
  data = list(
    Richness = Richness
  ),
  A = list(1, 1, A),
  effects = list(
    Intercept = rep(1, N),
    Xm = Xm[,-1],  # Use covariates without the intercept (includes dummy variables)
    w = w.index    # Ensure w.index is correctly aligned with random effects
  )
)
# 
# StackFit <- inla.stack(
#   remove.unused = TRUE,
#   tag = "Fit",
#   data = list(Shannon = vol.data$ShannonsDI),  # Ensure this variable exists in vol.data
#   A = list(1,1, A),  # Double-check that A is correctly sized for vol.data
#   effects = list(
#     Intercept = rep(1, nrow(vol.data)),
#     Xm = model.matrix(~ Habitat + Size + Control + Reference + Dissolved.Oxygen.z + Specific.Conductance.z - 1, data = vol.data), 
#     w = w.index
#   )
# )
# 
# # Use only relevant terms from Xm and provide the response separately
# StackFit <- inla.stack(
#  # remove.unused = TRUE,
#   tag = "Fit",
#   data = list(Shannon = ShannonDI),  # Response variable as separate vector
#   A = list(1, 1, A),  # A matrix for projections (make sure this aligns with w.index)
#   effects = list(
#     Intercept = rep(1, nrow(Xm)),  # Intercept term
#     Xm = Xm[,-1],  # Use the covariates from Xm without the intercept
#     w = w.index    # Ensure w.index aligns correctly
#   )
# )

#inla.stack.data(StackFit)

#Fit the global model
# 
# base.form.pcp <- ShannonsDI ~ Habitat + Size + Control+ Reference +
#   Dissolved.Oxygen.z + Specific.Conductance.z+  f(Stream,model="iid")

base.form.pcp <- Richness ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z +#+  f(Stream,model="iid")
  f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)



tic()
spat.inla.shannon<- inla(formula = update(base.form.pcp, . ~ . +
                                            f(Stream,model="iid")),
                         family = "poisson",
                         data = inla.stack.data(StackFit),
                         control.compute = list(waic = TRUE, config=TRUE,cpo=TRUE, po=TRUE,return.marginals.predictor=TRUE),
                         # control.family = control.family,
                         # control.inla = list(strategy = "gaussian"), 
                         control.predictor = list(link = 1, compute = TRUE,
                                                  A = inla.stack.A(StackFit)))

toc()

spat.inla.rich=spat.inla.shannon
summary(spat.inla.rich)



autoplot(spat.inla.rich)

#INLAutils::plot_inla_residuals(spat.inla.d2.d.m, observed)

par(mar=c(4,4,4,4))
D<-INLAutils::plot_inla_residuals(spat.inla.rich, Richness)



########Seems we still have trouble with the lowest values

SpatField.w <- inla.spde2.result(inla = spat.inla.rich,
                                 name = "w",
                                 spde = spde,
                                 do.transfer = TRUE)
Kappa <- inla.emarginal(function(x) x,
                        SpatField.w$marginals.kappa[[1]] )

Sigma_u <- inla.emarginal(function(x) sqrt(x),
                          SpatField.w$marginals.variance.nominal[[1]] )

Range <- inla.emarginal(function(x) x,
                        SpatField.w$marginals.range.nominal[[1]] )


result = inla.spde2.result(spat.inla.rich, "w", spde)

Range

par(mar=c(4,4,4,4))
plot(result[["marginals.range.nominal"]][[1]], type = "l",
     main = "Nominal range, posterior density")


LocMesh <- mesh2d$loc[,1:2]
D<-dist(LocMesh)
#D <- as.matrix(D)
#D
# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, #max(D)
             144392, length = 10000)
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1)
Cor.M[1] <- 1

Sigma_u
Range
max(D)

cor.plot.data=data.frame(d.vec,Cor.M)
cor.plot.data

cor.plot=ggplot()+geom_line(data=cor.plot.data,aes(x=d.vec,y=Cor.M))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlim(0,5000)+
  ylim(0,1)+
  geom_abline(aes(intercept = 0.05, slope = 0),linetype=3)+
  xlab("Distance (m)")+
  ylab("Matern Correlation Values")+
  theme(axis.title = element_text(face="bold"))+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))+
  geom_text(x = 1000, y = 0.1, aes(label = "5% autocorrelation"))

cor.plot
# 
# par(mfrow=c(1,1), mar = c(5,5,2,2))
# plot(x = d.vec,
#      y = Cor.M,
#      pch = 16,
#      type = "l",
#      cex.lab = 1.5,
#      xlab = "Distance",
#      ylab = "Correlation",
#      xlim = c(0, 10000))
# abline(h = 0.1, lty = 2)


cor.plot.data
Kappa

fitIndex <- inla.stack.index(StackFit, tag='Fit')$data


fitted.b<- spat.inla.rich$summary.fitted.values[fitIndex,]

#fitted.nonspat<-Chin.inla.null$summary.fitted.values[fitIndex,]



Pred   <- fitted.b[, "mean"]
VarY <- Pred
resid   <- (Richness - Pred) / sqrt(VarY)

MyData3 <- data.frame(resid = resid,
                      Xkm = utm_coords[,1],
                      Ykm = utm_coords[,2])


INLAutils::plot_inla_residuals(spat.inla.rich, observed)

# let's look at the residuals again
resid.plot=ggplot_inla_residuals(spat.inla.rich,Richness,CI = TRUE,binwidth = 0.1)

ggplot_inla_residuals2(spat.inla.rich,Richness, CI=TRUE,method = NA)


MyData3$MySize <- 2 * abs(MyData3$resid) / max(MyData3$resid)
MyData3$MyCol <- ifelse(MyData3$resid> 0, 1, 2)

#View(MyData3)
lattice::xyplot(MyData3$Ykm ~ MyData3$Xkm,
                data = MyData3,
                cex = MyData3$MySize,
                col = MyData3$MyCol,
                pch = 1)
par(mfrow=c(1,1))

hist(MyData3$resid, breaks = 20)

summary(spat.inla.rich)
#After accounting for the random effects of Site and Space 
#for no tested effects were included in our best model.


#fitted.d=nospat.inla.t.m$summary.fitted.values[fitIndex,]
rmse <- sqrt(mean((Richness - fitted.b$mean)^2))
rmse

#



fitIndex <- inla.stack.index(StackFit, tag='Fit')$data
fitted <- spat.inla.rich$summary.fitted.values[fitIndex,]

Pi1=fitted$mean

# Pi1
# 
# observed
# 

D1  <- (Richness - Pi1) / sqrt(Pi1)

summary(D1)



fitted.d<- spat.inla.rich$summary.fitted.values
RRtot=sum((observed-mean(observed))^2)
RRes=sum((observed-fitted.b$mean)^2)
pseudo_r2_val.d=1-RRes/RRtot
pseudo_r2_val.d
txt="R^{2} == 0.706"


CSD2=data.frame(flow.data,fitted.b)
# 49%% of the variance is accounted for by random effects
max(CSD2$mean)

plot=ggplot(CSD2,aes(x=mean,y=Richness))+
  geom_point()+
  theme_bw()+
  xlim(c(0,8))+
  ylim(c(0,8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  
  geom_abline(aes(intercept = 0, slope = 1))+
  xlab("Fitted native fish species richness")+
  ylab("Observed native fish species richness")+
  theme(axis.title = element_text(face="bold"))

plot=plot+geom_text(x =2, y = 6, aes(label = txt),parse=T)+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
plot



# Generate posterior predictive replications
n_replications <- 1000  # Number of replicated datasets
predictions <- inla.posterior.sample(n_replications, spat.inla.rich)

# Calculate observed discrepancies
obs_discrepancy <- abs(flow.data$Richness - spat.inla.rich$summary.fitted.values$mean)

# Calculate discrepancies for replicated data
rep_discrepancies <- sapply(predictions, function(rep) {
  replicated_data <- rep$latent[flow.data$Richness]  # Use appropriate index for observed data points
  abs(flow.data$Richness - replicated_data)
})

# Calculate Bayesian p-value
bayesian_p_value <- mean(rowMeans(rep_discrepancies) > obs_discrepancy)
print(bayesian_p_value)

# 
# # Create the plot (for a subset of replications, if necessary)
# discrepancy_data <- data.frame(
#   Observed = obs_discrepancy[1:42],
#   Replicated = rep_discrepancies[, sample(1:n_replications, length(obs_discrepancy))]
# )
# 
# ggplot(discrepancy_data, aes(x = Observed, y = Replicated.205)) +
#   geom_point(color = 'blue') +
#   geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed") +
#   labs(
#     title = "Bayesian p-value Comparison",
#     x = "Observed Discrepancy",
#     y = "Replicated Discrepancy"
#   ) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 


summary(spat.inla.rich)



exp(0.842)

exp(2.064)


1-exp(-0.917)
1-exp(-0.167)



# Generate the plot with Site.Name as the x-axis labels
ggplot(flow.data, aes(x = reorder(Site.Name, Richness), y = Richness, fill = Habitat_Label)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Habitat_Label, scales = "free_x") +  # Facet by Habitat_Label (Reservoir/Stream)
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels for readability
    strip.text = element_text(size = 12),  # Enhance facet label size
    panel.grid.major.x = element_blank()  # Remove x-axis major grid lines for a cleaner look
  ) +
  labs(
    title = "Species Richness by Site and Habitat Type",
    x = "Site Name",
    y = "Richness"
  ) +
  scale_fill_manual(values = c("#00577b", "#97a641"))  # Adjust colors as needed





# After accounting for the random effects of waterbody and spatial position....

# Streams are expected to have 2.3 - 7.9 times greater native fish richness than 
#reservoirs.

# Tributaries are expected to have between 15 and 60% less native fish richness 
#than mainstem river sites.

# Because our posterior probability distributions spanned zero for Control/Impact,
# Reference/Impact, Dissolved Oxygen Concentration, and Specific Conductance, there
# is no evidence that these variables described a meaningful amount of variance in
#native fish species richness.


# Our model described approximately 64% of the variance in native fish species
# richness, had normally distributed residuals, and had a bayesian p-value between
# 0.1 and 0.9 (0.45), indicating that the model fit the data well.





library(ggplot2)
library(patchwork)
library(cowplot)
# AUC Plot B (ROC Curve using ggplot2 with step function, flipped axes)
roc_coords <- coords(roc_curve, ret = c("specificity", "sensitivity"), transpose = FALSE)

roc_df <- data.frame(
  specificity = roc_coords$specificity,
  sensitivity = roc_coords$sensitivity
)

# Extract fixed effect summary statistics
fixed_effects <- spat.inla.diversity$summary.fixed

# Extract the posterior mean, 2.5% and 97.5% credible intervals
forest_data <- data.frame(
  Variable = rownames(fixed_effects),
  OR = exp(fixed_effects$mean),  # Convert log-odds to odds ratio
  CI_Lower = exp(fixed_effects$`0.025quant`),  # 2.5% credible interval
  CI_Upper = exp(fixed_effects$`0.975quant`)  # 97.5% credible interval
)

print(forest_data)  # View extracted data

# Remove the intercept and reformat the variable names without reentering numeric values
forest_data <- forest_data[-1, ]  # Remove the intercept row

# Reverse the order of Variable if not already reversed
forest_data$Variable <- factor(forest_data$Variable, levels = rev(forest_data$Variable))

# Reorder and rename the Poisson model variables as per your requirements
forest_data$Variable <- recode(forest_data$Variable,
                                       "Habitat" = "Habitat (reservoir vs. stream)",
                                       "Size" = "Stream size (mainstem vs. tributary)",
                                       "Reference" = "Reference sites vs. impact sites",
                                       "Control" = "Control sites vs. impact sites",
                                       "Dissolved.Oxygen.z" = "Scaled dissolved oxygen (mg/L)",
                                       "Specific.Conductance.z" = "Scaled specific conductance (\u03BCS/cm)")

# Add asterisks to rows where the credible interval does not span 1 (Poisson model)
forest_data$Label <- ifelse(forest_data$CI_Lower > 1 | forest_data$CI_Upper < 1, "*", "")

# Extract fixed effect summary statistics for the Poisson model
fixed_effects_poisson <- spat.inla.rich$summary.fixed

# Extract the posterior mean, 2.5% and 97.5% credible intervals for the Poisson model
forest_data_poisson <- data.frame(
  Variable = rownames(fixed_effects_poisson),
  OR = exp(fixed_effects_poisson$mean),  # Convert log-odds to odds ratio (Poisson case)
  CI_Lower = exp(fixed_effects_poisson$`0.025quant`),  # 2.5% credible interval
  CI_Upper = exp(fixed_effects_poisson$`0.975quant`)  # 97.5% credible interval
)

# Remove the intercept and reformat the variable names without reentering numeric values
forest_data_poisson <- forest_data_poisson[-1, ]  # Remove the intercept row

# Reverse the order of Variable if not already reversed
forest_data_poisson$Variable <- factor(forest_data_poisson$Variable, levels = rev(forest_data_poisson$Variable))

# Reorder and rename the Poisson model variables as per your requirements
forest_data_poisson$Variable <- recode(forest_data_poisson$Variable,
                                       "Habitat" = "Habitat (reservoir vs. stream)",
                                       "Size" = "Stream size (mainstem vs. tributary)",
                                       "Reference" = "Reference sites vs. impact sites",
                                       "Control" = "Control sites vs. impact sites",
                                       "Dissolved.Oxygen.z" = "Scaled dissolved oxygen (mg/L)",
                                       "Specific.Conductance.z" = "Scaled specific conductance (\u03BCS/cm)")

# Add asterisks to rows where the credible interval does not span 1 (Poisson model)
forest_data_poisson$Label <- ifelse(forest_data_poisson$CI_Lower > 1 | forest_data_poisson$CI_Upper < 1, "*", "")

# Add a column to distinguish between models (Binomial and Poisson)
forest_data$Model <- "Non-Extreme Low Diversity (Binomial)"
forest_data_poisson$Model <- "Species Richness (Poisson)"

# Combine both datasets
combined_forest_data <- rbind(forest_data, forest_data_poisson)

# Reorder Variable to maintain the current order
combined_forest_data$Variable <- factor(combined_forest_data$Variable, levels = rev(unique(combined_forest_data$Variable)))

# Reorder the variables in both forest plots (A and C)
forest_data$Variable <- factor(forest_data$Variable, levels = c("Habitat (reservoir vs. stream)",
                                                               "Stream size (mainstem vs. tributary)",
                                                                "Reference sites vs. impact sites",
                                                                "Control sites vs. impact sites",
                                                                "Scaled dissolved oxygen (mg/L)",
                                                                "Scaled specific conductance (\u03BCS/cm)"))  # Specify the correct levels
forest_data_poisson$Variable <- factor(forest_data_poisson$Variable, levels = c("Habitat (reservoir vs. stream)",
                                                                                "Stream size (mainstem vs. tributary)",
                                                                                "Reference sites vs. impact sites",
                                                                                "Control sites vs. impact sites",
                                                                                "Scaled dissolved oxygen (mg/L)",
                                                                                "Scaled specific conductance (\u03BCS/cm)"))  # Match levels

# Update the AUC Plot B (ROC Curve with consistent font and size for text)
auc_plot <- ggplot(roc_df, aes(x = specificity, y = sensitivity)) +
  geom_step(color = "black", size = 1) +
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "grey") +
  annotate("text", x = 0.00, y = 0.05, label = "b", size = 8, fontface = "bold", hjust = 0) +
  annotate("text", x = 0.6, y = 0.2, label = "AUC = 0.96", size = 6, fontface = "plain") +  # Standardized font and size
  theme_classic() +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)
  )

# Update Forest Plot A (Non-Extreme Low Diversity)
plot_binomial <- ggplot(forest_data, aes(x = Variable, y = OR)) +
  geom_pointrange(aes(ymin = CI_Lower, ymax = CI_Upper), size = 1, color = "black") +
  geom_text(aes(label = Label), nudge_x = 0.2, color = "black", size = 8) +
  coord_flip() +
  theme_classic() +
  labs(y = "Odds ratio (non-extreme low diversity)",
       x = NULL) +
  annotate("text", x = 1, y = 1e30, label = "a", size = 8, fontface = "bold", hjust = 1) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  scale_y_log10() +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

# Update Forest Plot C (Species Richness)
plot_poisson <- ggplot(forest_data_poisson, aes(x = Variable, y = OR)) +
  geom_pointrange(aes(ymin = CI_Lower, ymax = CI_Upper), size = 1, color = "black") +
  geom_text(aes(label = Label), nudge_x = 0.2, color = "black", size = 8) +
  coord_flip() +
  theme_classic() +
  labs(x = "Variable", y = "Rate ratio (species richness)") +
  annotate("text", x = 6, y = 8, label = "c", size = 8, fontface = "bold", hjust = 1) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(hjust = 1.2),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  scale_y_log10() +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

# Update Richness Plot D (Observed vs Fitted Richness)
plot_richness <- ggplot(CSD2, aes(x = mean, y = Richness)) +
  geom_point() +
  theme_bw() +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  xlab("Fitted native fish species richness") +
  ylab("Observed native fish species richness") +
  annotate("text", x = 0.1, y = 7.6, label = "d", size = 8, fontface = "bold", hjust = 0) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  geom_text(x = 2, y = 6, aes(label = txt), parse = TRUE) +
  theme(
    strip.text.x = element_text(size = 20)
  )

# Combine the plots using patchwork with consistent spacing
combined_plot <- (plot_binomial + auc_plot) / (plot_poisson + plot_richness) +
  plot_layout(guides = "collect") & theme(legend.position = 'none')

# Display the final plot
print(combined_plot)




