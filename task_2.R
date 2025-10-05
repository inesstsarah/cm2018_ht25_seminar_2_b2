# Imports and others
# Install packages
install.packages("pracma")

# Libraries
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


# Compare the fasted and non-fasted plots
library(ggplot2)
library(pracma)
data <- read.csv("data/Data_T2.csv")

# Plot data
ggplot(data, aes(x = Time, y = GE_fed)) + 
  geom_point() + 
  labs(title = "Gastric content over time", x = "Time (minutes)", y = "Gastric content (volume)") +
  theme_minimal()

# ------------------ Fit models ------------------------------------------------
exp_fit <- nls(GE_fed ~ GE0 * exp(-k * Time),
               data = data,
               start = list(GE0 = data$GE_fed[1], k = 0.01))

exp_lag_fit <- nls(GE_fed ~ ifelse(Time <= tlag, GE0, GE0 * exp(-k * (Time - tlag))),
                   data = data,
                   start = list(GE0 = data$GE_fed[1], k = 0.01, tlag = 2))

weib_fit <- nls(GE_fed ~ GE0 * exp(- (Time/alpha)^beta),
                data = data,
                start = list(GE0 = data$GE_fed[1], alpha = 50, beta = 1.2))

# ---------------- Exponential model -------------------------------------------
# Get predicted GE_fed values from the exponential model
pred_exp <- predict(exp_fit, newdata = data)

# Plot data against predited values
ggplot(data, aes(x = Time, y = GE_fed)) +
  geom_point() +
  geom_line(aes(x = Time, y = pred_exp), color = "orchid", linewidth = 1) +
  labs(title = "Fed gastric emptying with fitted exponential model",
       x = "Time (minutes)",
       y = "Gastric content") +
  theme_minimal()

# ---------------- Lagged exponential  model -----------------------------------
# Get predicted GE_fed values from the lagged exponential model
pred_exp_lag <- predict(exp_lag_fit, newdata = data)

# Plot
ggplot(data, aes(x = Time, y = GE_fed)) +
  geom_point() +
  geom_line(aes(y = pred_exp_lag), color = "orchid", linewidth = 1) +
  labs(title = "Fed gastric emptying with fitted lagged exponential model",
       x = "Time (minutes)",
       y = "Gastric content") +
  theme_minimal()

# ---------------- Weibull model -----------------------------------------------
# Get predicted GE_fed values from the weibull model
pred_weib <- predict(weib_fit, newdata = data)

# Plot data against predicted values
ggplot(data, aes(x = Time, y = GE_fed)) +
  geom_point() +
  geom_line(aes(x = Time, y = pred_weib), color = "orchid", linewidth = 1) +
  labs(title = "Fed gastric emptying with fitted exponential model",
       x = "Time (minutes)",
       y = "Gastric content") +
  theme_minimal()

# -------------- Comparison between models -------------------------------------
# Compare fits
fits <- list(exp = exp_fit, lag_exp = exp_lag_fit, weibull = weib_fit)
tibble(
  model = names(fits),
  AIC = sapply(fits, AIC),
  BIC = sapply(fits, BIC),
  RMSE = sapply(fits, function(m) sqrt(mean(resid(m)^2)))
)

# --------------- Secondary parameters -----------------------------------------

# Extract secondary parameter for half-life (weibull)
coef_weib <-  coef(weib_fit)

# Calculate half-life (weibull)
t_half_weib <- coef_weib["alpha"] * (log(2)) ^ (1 / coef_weib["beta"])


# Tlag - first time point where gastric content drops below 99% of initial value
initial_value <- data$GE_fed[1]
Tlag <- min(data$Time[data$GE_fed < 0.99 * initial_value])

# Tmax - time of maximum gastric content
Tmax <-  data$Time[which.max(data$GE_fed)]


# AUC - area under curve
AUC <- trapz(data$Time, data$GE_fed)


# Times for different observed volumes, based on initial value (not max value)
t_at_retention <- function(z) {
  approx(x = data$GE_fed/100, y = data$Time, xout = z)$y
}

t90 <- t_at_retention(0.9)
t50 <- t_at_retention(0.5)
t10 <- t_at_retention(0.1)


