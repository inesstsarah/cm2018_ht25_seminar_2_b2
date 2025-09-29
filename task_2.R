# Install packages
install.packages("pracma")

# Libraries
library(tidyverse)
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


