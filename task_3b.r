# TASK3 Part 1 ---> GLM version ( Generalized linear models)

#_______________________________________________________________________________________
# Predict the dose that produces 80% gastric retention at 30 minutes.
#_______________________________________________________________________________________

#install.packages(("msm"))

# Load the data
data <- read.csv("Data_T3.csv", nrows = 30)

#Showing the data
head(data)

# Convert retention to proportions (0-1 scale) since 
# glm() with binomial family expects values between 0 and 1

data$retention_prop <- data$GE_content / 100

# Fit GLM with logistic link 
# This is exactly like: glm(y ~ x1, family = "binomial", data = binom_dat)

# EXPLANATION OF FUNCTION
# retention_prop ~ Dose --> GE content is dependent on dose
# family binomial       --> logistic regression - maps probabilities 0-1 to the real line

model <- glm(retention_prop ~ Dose, family = binomial(link = "logit"), data = data) # coefficients β0 and β1 estimated by maximum likelihood

# Showing the model
print(summary(model))   # Gives coefficients, errors, p-values etc 

# Define inverse logit function (from lecture)
# invlogit transforms back from logit scale to probability scale since the
# GLM coefficients are on the logit scale
invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}

# Find dose for 80% retention
# We need to solve: 0.80 = invlogit(β0 + β1*Dose)
# Rearranging: logit(0.80) = β0 + β1*Dose
# So: Dose = (logit(0.80) - β0) / β1

logit_80 <- log(0.80 / (1 - 0.80))  # logit(0.80) = 1.386
intercept <- coef(model)[1]         # we extract β0 and slope β1
slope <- coef(model)[2]

dose_for_80 <- (logit_80 - intercept) / slope    # calculate the dose for 80%

# Answer
print(paste("Dose for 80% gastric retention:", round(dose_for_80, 2), "mg"))




# Creating prediction data for plotting
#________________________________________________________________________________________

# From dose 0-150 with 200 points for smooth plotting
# pred_retention predicts probability at each dose using the inverse logit

dose_seq <- seq(min(data$Dose), (dose_for_80*1.5), length.out = 400)       
pred_retention <- invlogit(coef(model)[1] + coef(model)[2] * dose_seq)


# Divide-by-4 rule: 
# in logistic regression - maximum slope in porbability scale is β1/4
# purpose: gives a quick sense of how fast the probability changes per dose
max_diff <- slope / 4
print(paste("Divide-by-4 rule: Maximum probability change per unit dose =",round(max_diff, 3)))

# Plot on percentage scale 
plot(data$Dose, data$GE_content, 
     xlab = "Dose (mg)", 
     ylab = "Gastric Retention (%)",
     main = "Dose-Response Curve",
     pch = 19, col = "blue", ylim = c(0,100), xlim = c(0,dose_for_80))

lines(dose_seq, pred_retention * 100, col = "blue", lwd = 2) # pred_retention*100 for % again

abline(h = 80, col = "red", lty = 2, lwd = 2)
abline(v = dose_for_80, col = "red", lty = 2, lwd = 2)

text(dose_for_80, 90, 40, labels = paste("Dose =", round(dose_for_80, 2), "mg"), pos = 2, col = "red")

# Add legend
legend("bottomright", 
       legend = c("Observed data", "GLM fit", "80% target"),
       col = c("blue", "blue", "red"),
       pch = c(19, NA, NA),
       lty = c(NA, 1, 2),
       lwd = c(NA, 2, 2),
       bty = "n")


# Added code for computing the 95% CI :
#_________________________________________________________________________________________
# Extract coefficients and covariance matrix
coefs <- coef(model)
vcov_mat <- vcov(model)

# Target on logit scale
logit_80 <- log(0.80 / (1 - 0.80))

# Delta method for SE of the dose
library(msm)  # for deltamethod
dose_se <- deltamethod(~ (logit_80 - x1) / x2, coefs, vcov_mat)

# Point estimate (needed to center CI)
dose_for_80 <- (logit_80 - coefs[1]) / coefs[2]

# 95% CI
dose_CI_lower <- dose_for_80 - 1.96 * dose_se
dose_CI_upper <- dose_for_80 + 1.96 * dose_se

# Print CI only
print(paste("95% CI for dose at 80% retention: [",
            round(dose_CI_lower, 2), "-", round(dose_CI_upper, 2), "] mg"))

# ANSWER: 95% CI for dose at 80% retention: [ 59.59 - 148.35 ] mg




# Compare null model vs full model (Chi-squared test)
#______________________________________________________________________________________
null_model <- glm(retention_prop ~ 1, family = binomial(link = "logit"), data = data)
anova(null_model, model, test = "Chisq")

# Deviance residuals
resid_dev <- residuals(model, type = "deviance")

# Plot residuals vs fitted values
plot(fitted(model), resid_dev, 
     xlab = "Fitted values (predicted retention)", 
     ylab = "Deviance residuals",
     main = "Deviance Residuals vs Fitted")
abline(h = 0, lty = 2, col = "red")



#######################################################################################################################
# Task 3B part2
#______________________________________________________________________________________
# Comparing the predictions of model made on dosage 5-20mg with real data
#______________________________________________________________________________________


# Model made on lower dosage + new data
#______________________________________________________________________________________
#--------------------------------------------------------------------------------------

# Read data
data <- read.csv("Data_T3.csv", nrows = 30)
newData <- read.csv("Data_T3B_Adjusted_format.csv")

# Convert retention to proportions (0-1 scale) since 
data$retention_prop <- data$GE_content / 100
newData$retention_prop <- newData$GE_content / 100

# Create model 
model <- glm(retention_prop ~ Dose, family = binomial(link = "logit"), data = data)

# Showing the model
print(summary(model)) 

# Define inverse logit function
invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}

# Find dose for 80% retention
logit_80 <- log(0.80 / (1 - 0.80))  
intercept <- coef(model)[1]         
slope <- coef(model)[2]

dose_for_80 <- (logit_80 - intercept) / slope

# Answer
print(paste("Dose for 80% gastric retention:", round(dose_for_80, 2), "mg"))

# Creating prediction data for plotting
dose_seq <- seq(min(data$Dose), (160), length.out = 400)       
pred_retention <- invlogit(coef(model)[1] + coef(model)[2] * dose_seq)



# Plotting model curve and raw data
#_________________________________________________________________________________________

# Plot model predictions
plot(dose_seq, pred_retention,
     type = "l", lwd = 2, col = "blue",
     xlab = "Dose (mg)", ylab = "Gastric Retention (proportion)",
     ylim = c(0,1),
     main = "Model Prediction vs Raw Data")

# Add raw original data (training range 5–20 mg)
points(data$Dose, data$retention_prop,
       col = "black", pch = 16)

# Add raw new dataset (higher dose validation)
points(newData$Dose, newData$retention_prop,
       col = "red", pch = 17)

# Add a legend
legend("bottomright",
       legend = c("Model Prediction", "Training Data (5–20 mg)", "New Data (100-150 mg)"),
       col = c("blue", "black", "red", "orange"),
       lwd = c(2, NA, NA),
       pch = c(NA, 16, 17))

# Add horizontal line at 80% retention
abline(h = 0.8, col = "darkgreen", lty = 2, lwd = 2)

# Add vertical line at predicted dose for 80%
abline(v = dose_for_80, col = "purple", lty = 2, lwd = 2)





# Raw dose that actually gives 80% retention
#_________________________________________________________________________________________
threshold <- 0.8
cross_idx <- which(newData$retention_prop >= threshold)[1]

if (!is.na(cross_idx)) {
  observed_dose_for_80 <- newData$Dose[cross_idx]
  
  # Add vertical line at observed dose
  abline(v = observed_dose_for_80, col = "orange", lty = 2, lwd = 2)
}

print(paste("Dose for 80% gastric retention from raw data:", observed_dose_for_80, "mg"))


#####################################################################################################################
# Task3B Refined model (2 ways)
#________________________________________________________________________________________________

#install.packages("numDeriv")

library(dplyr)
library(tidyr)
library(ggplot2)
library(drc)

# Read data 
oldData <- read.csv("Data_T3.csv")
newData <- read.csv("Data_T3B_Adjusted_format.csv")

names(oldData)

oldData <- oldData %>%
  dplyr::select(Dose, GE_content) %>%          # To skip row "X"
  mutate(retention_prop = GE_content / 100)

# Match the columns
newData <- newData %>%
  dplyr::select(Dose, GE_content, retention_prop)

# Combine data
allData <- bind_rows(oldData, newData)

# Remove NAN rows
allData <- allData %>% filter(!is.na(Dose), !is.na(retention_prop))

# See the new data set
print(head(allData))
print(table(allData$Dose))


# Creating the model (with all data this time) (GLM)
#_________________________________________________________________________________________
model <- glm(retention_prop ~ Dose, data = allData, family = binomial)

# Print model 
print(summary(model))

# Predict up to 160 mg
dose_seq <- seq(min(allData$Dose), 160, length.out = 400)
pred_glm <- predict(model, newdata = data.frame(Dose = dose_seq), type = "response")

# Dose at 80% retention
logit_80 <- log(0.8 / (1 - 0.8))
intercept <- coef(model)[1]
slope <- coef(model)[2]
dose_for_80_glm <- (logit_80 - intercept) / slope

print(paste("Dose for 80% gastric retention:", round(dose_for_80_glm, 2), "mg"))


# Plot 
#______________________________________________________________________________________________
plot(allData$Dose, allData$retention_prop,
     col = "red", pch = 16,
     xlab = "Dose (mg)", ylab = "Gastric Retention (proportion)",
     main = paste("GLM: 80% retention at ~", round(dose_for_80_glm,1), "mg"),
     ylim = c(0,1))
lines(dose_seq, pred_glm, col = "blue", lwd = 2)
abline(h = 0.8, col = "darkgreen", lty = 2)
abline(v = dose_for_80_glm, col = "purple", lty = 2)





#------------------------------------------------------------------------------
# 4-Parameter Logistic (4PL) Model (Using DRM package)
#------------------------------------------------------------------------------
model_4pl <- drm(retention_prop ~ Dose, data = allData, fct = LL.4())

print(summary(model_4pl))

# Predict with 4PL
pred_4pl <- data.frame(Dose = dose_seq,
                       retention = predict(model_4pl, newdata = data.frame(Dose = dose_seq)))

# Dose at 80% retention from 4PL
ED_80 <- ED(model_4pl, 0.8, type = "absolute", interval = "delta")
print(paste("4PL predicted dose for 80% retention:", ED_80[1], "mg\n"))


# Plot 
#______________________________________________________________________________________________

plot(allData$Dose, allData$retention_prop,
     col = "red", pch = 16,
     xlab = "Dose (mg)", ylab = "Gastric Retention (proportion)",
     main = paste("4PL: 80% retention at ~", round(ED_80[1],1), "mg"),
     ylim = c(0,1))
lines(dose_seq, pred_4pl$retention, col = "blue", lwd = 2)
abline(h = 0.8, col = "darkgreen", lty = 2)
abline(v = ED_80[1], col = "purple", lty = 2)

print(table(allData$Dose))




##################################################################################################

# Estimating recommended dose based on derivative of the glm curve

library(numDeriv)

beta1 <- coef(model)[2]

dose_seq <- seq(min(allData$Dose), 160, length.out = 400)
pred_glm <- predict(model, newdata = data.frame(Dose = dose_seq), type = "response")

# analytisk derivata
deriv_glm <- beta1 * pred_glm * (1 - pred_glm)

# dos där lutningen är störst
dose_max_slope <- dose_seq[which.max(deriv_glm)]
max_slope_value <- max(deriv_glm)

print(paste("Dose with max slope (mg):", dose_max_slope))
print(paste("Max slope of GLM curve (proportion per mg):", max_slope_value))

# Answer on dosage at max slope: 41.1278195488722 mg















