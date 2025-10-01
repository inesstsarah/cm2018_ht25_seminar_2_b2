# Imports and others
library(tidyverse)
library(ggplot2)

install.packages("PKNCA")
library(PKNCA)
install.packages("Metrics")
library(Metrics)

install.packages("hydroGOF")
library(hydroGOF)

df = read.csv("data/Data_T2.csv")
selected_columns <- df[c("Time","GE_fasted")]

time_col = selected_columns["Time"]
fasted_col = selected_columns["GE_fasted"]

df
# Sample data
x <- df$Time
y <- df$GE_fasted

# Basic scatter plot
plot(x, y, main = "Gastric Content", xlab = "Time (minutes)", ylab = "Volume of Gastric Content (%)", pch = 19, col = "blue")


# Fit the exponential decay model
# Initial parameter estimates
start_vals <- list(C0 = max(y), lambda = 0.5)
start_vals

# Fit the model using nls()
fit <- nls(y ~ C0 * exp(-lambda * x), data = selected_columns, start = start_vals)
summary(fit)
fit

# Draw curve
plot(x, y, main = "Exponential Decay Fit for Gastric Content", pch = 16, xlab = "Time (minutes)", ylab = "Volume of Gastric Content (%)")
curve(predict(fit, newdata = data.frame(x = x)), add = TRUE, col = "red")

prediction = predict(fit, newdata = data.frame(x = x))
k <- coef(fit)["lambda"] # Extract the slope for time


# Calculate half-life, take the natural logarithm 
t_half <- log(2) / k
t_half
# Half life is calculated at 15.69482 minutes

# Plot residual error (RMSE -> root mean square error) and assess model performance 
prediction
rmse_regression <- sqrt(mean((y - prediction)^2))
rmse_regression
# Low RMSE of 4.263958 indicates that the predictions 
# are close to values
# Get NRMSE
nrmse_regression <- nrmse(prediction, y)
nrmse_regression
# The result of the NRMSE is 15.8%, which indicates
# a good fit

# Total Sum of Squares (TSS)
tss <- sum((y - mean(y))^2)

# Residual Sum of Squares (RSS)
rss <- sum((y - prediction)^2)

# R² formula
r_squared <- 1 - (rss / tss)
r_squared
# R-squared is 0.974671


# AUC can be used to describe the clearance of materials
# from the rat's gastric content
params <- coef(fit)
y0 <- params["C0"]
k <- params["lambda"]
fitted_fun <- function(t) y0 * exp(-k * t)

#


# AUC from t = 0 to t = 60
auc <- integrate(fitted_fun, lower = 0, upper = 60)$value
print(auc)

# Linear model with log transformation
log_y = log(y)
plot(x, log_y, main = "Gastric Content", xlab = "Time (minutes)", ylab = "Volume of Gastric Content (%)", pch = 19, col = "blue")
df
lm_model = lm(log(GE_fasted)~Time, data = df)
abline(lm_model, col = "red", lwd = 2) # col sets the line color, lwd sets the line width
lm_model

# Extract decay rate (negative slope)
k <- -coef(lm_model)["Time"]

# Calculate half-life
t_half <- log(2) / k

# Print result
print(t_half)
# The half life is 15.67384 minutes

predicted_values <- predict(lm_model)

# Calculate RMSE
rmse_lm <- sqrt(mean((log_y - predicted_values)^2))
rmse_lm
# The RMSE is 0.1000484

# Calculate NRMSE
nrmse_regression <- nrmse(predicted_values, log_y)
nrmse_regression
# The NRMSE is 12.8%

# R2
r_squared <- summary(lm_model)$r.squared
r_squared

summary(lm_model)
# The r_squared is closer to 1 (it is 0.9832208) which means it
# is a good fit


# Get AUC of model
fitted_fun <- function(t) y0 * exp(-k * t)

# AUC from t = 0 to t = 60
auc <- integrate(fitted_fun, lower = 0, upper = 60)$value
print(auc)

# 