# Imports and others
library(tidyverse)
library(ggplot2)



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


# Draw curve
plot(x, y, main = "Exponential Decay Fit for Gastric Content", pch = 16, xlab = "Time (minutes)", ylab = "Volume of Gastric Content (%)")
curve(predict(fit, newdata = data.frame(x = x)), add = TRUE, col = "red")

prediction = predict(fit, newdata = data.frame(x = x))


# Get the half life (still not sure how to do it)
prediction 

# Plot residual error (RMSE -> root mean square error) and assess model performance 


# Add other models and compare fits