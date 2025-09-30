# TASK3 Part 2 ---> GLM version ( Generalized linear models)

#_______________________________________________________________________________________
# Predict the dose that produces 80% gastric retention at 30 minutes.
#_______________________________________________________________________________________

#install.packages(("msm"))

# Load the data
data <- read.csv("Data_T3.csv")

# Showing the data
head(data)

# Convert retention to proportions (0-1 scale) since 
# glm() with binomial family expects values between 0 and 1

data$retention_prop <- data$GE_content / 100

# Fit GLM with logistic link 
# This is exactly like: glm(y ~ x1, family = "binomial", data = binom_dat)

# EXPLANATION OF FUNCTION
# retention_prop ~ Dose --> GE content is dependent on dose
# family binomial       --> logistic regression - maps probabilities 0-1 to the real line

model <- glm(retention_prop ~ Dose, family = binomial(link = "logit"), data = data)

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

# THE ANSWER
print(paste("Dose for 80% gastric retention:", round(dose_for_80, 2), "mg"))

# Creating prediction data for plotting
# From dose 0-150 with 200 points for smooth plotting
# pred_retention predicts probability at each dose using the inverse logit

dose_seq <- seq(0, 150, length.out = 200)       
pred_retention <- invlogit(coef(model)[1] + coef(model)[2] * dose_seq)

# Plotting the results - scatter plot of observed data (Dose vs retention)
plot(data$Dose, data$retention_prop, 
     xlab = "Dose (mg)", 
     ylab = "Gastric Retention (proportion)",
     main = "Dose-Response Curve (GLM)",
     pch = 19, col = "blue", ylim = c(0, 1))

# Add fitted curve
lines(dose_seq, pred_retention, col = "blue", lwd = 2)

# Add reference lines
abline(h = 0.80, col = "red", lty = 2, lwd = 2)
abline(v = dose_for_80, col = "red", lty = 2, lwd = 2)

# Add legend
legend("bottomright", 
       legend = c("Data", "GLM fit", "80% target"),
       col = c("blue", "blue", "red"),
       pch = c(19, NA, NA),
       lty = c(NA, 1, 2),
       lwd = c(NA, 2, 2))

# Divide-by-4 rule: 
# in logistic regression - maximum slope in porbability scale is β1/4
# purpose: gives a quick sense of how fast the probability changes per dose
max_diff <- slope / 4
print(paste("Divide-by-4 rule: Maximum probability change per unit dose =", 
            round(max_diff, 3)))

# Plot on percentage scale 
plot(data$Dose, data$GE_content, 
     xlab = "Dose (mg)", 
     ylab = "Gastric Retention (%)",
     main = "Dose-Response Curve",
     pch = 19, col = "blue")

lines(dose_seq, pred_retention * 100, col = "blue", lwd = 2) # pred_retention*100 for % again
abline(h = 80, col = "red", lty = 2, lwd = 2)
abline(v = dose_for_80, col = "red", lty = 2, lwd = 2)

text(dose_for_80, 20, 
     labels = paste("Dose =", round(dose_for_80, 2), "mg"),
     pos = 4, col = "red")


# Added code for computing the 95% CI :
#_________________________________________________________________________________________
# Extract coefficients and covariance matrix
coefs <- coef(model_glm)
vcov_mat <- vcov(model_glm)

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




# Summary of the approach:
#______________________________________________________________________________________ 
#  Convert % gastric retention to proportion.

#  Fit a logistic regression (GLM) of retention vs dose.

#  Use the inverse logit to convert logit-scale predictions back to probabilities.

#  Solve the logistic equation to find the dose that gives 80% retention.

#  Visualize the dose-response relationship and reference dose.

# Compute the 95% CI 


