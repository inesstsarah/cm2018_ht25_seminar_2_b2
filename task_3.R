# ==============================================================================

# CM2018 HT25-1 Statistics for Medical Engineering 7.5 credits
# Task 3 - Dose–response at 30 min


# PART A =======================================================================

# --- SETUP --------------------------------------------------------------------

# Install Packages (install once if needed)
#install.packages(c("readr","dplyr","ggplot2"))
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# --- FUNCTION -----------------------------------------------------------------

# Returns the dose at the specified GE_content value, given beta factors
# (intercept and slope, respectively).
getDoseAtGE <- function(target_GE, b0, b1) {
  (target_GE - b0) / b1
}

# --- LOAD DATA & QUICK SUMMARY ------------------------------------------------

# Load data (long format)
data_3A <- read_csv("Data_T3.csv", show_col_types = FALSE, n_max = 30)

# Rename columns
names(data_3A) <- c("ID", "dose", "GE_content")

# Group by dose and calculate mean of GE_content
data_3A_by_dose <- data_3A %>%
  group_by(dose) %>%
  summarise(mean_GE_content = mean(GE_content))

# --- LINEAR REGRESSION --------------------------------------------------------

# Fitting Linear Models (linear regression)
m_3A <- lm(GE_content ~ dose, data = data_3A)

summary(m_3A)

# Beta values
b0_3A <- coef(m_3A)[1]   # intercept
b1_3A <- coef(m_3A)[2]   # slope

# Solve analytically for the dose that gives 50% and 80% retention
dose_3A_GE50 <- getDoseAtGE(50, b0_3A, b1_3A)
dose_3A_GE80 <- getDoseAtGE(80, b0_3A, b1_3A)

# --- PLOTS --------------------------------------------------------------------

# Histogram of gastric emptying content
hist(data_3A$GE_content,
     main = "Histogram of GE content",
     xlab = "GE content",
     col = "lightblue",
     border = "black")
title("Part A", outer = TRUE, line = -1)

# Boxplot of GE content by dose
boxplot(GE_content ~ dose, data = data_3A,
        main = "Boxplot of GE content by dose [mg]",
        xlab = "Dose [mg]",
        ylab = "GE content",
        col = "lightgreen")
title("Part A", outer = TRUE, line = -1)

# Scatterplot of dose [mg] vs GE content
plot(data_3A$dose, data_3A$GE_content,
     main = "Scatterplot: Dose [mg] vs GE content",
     xlab = "Dose [mg]",
     ylab = "GE content",
     pch = 19, col = "blue")
title("Part A", outer = TRUE, line = -1)

# Line plot of mean GE content
plot(
  data_3A_by_dose$dose,
  data_3A_by_dose$mean_GE_content,
  xlab = "Dose",
  ylab = "Mean GE Content",
  main = "Mean GE Content by Dose",
  type = "b",
  pch = 19
)
title("Part A", outer = TRUE, line = -1)

# Find min/max for axes so lines fit
x_min_3A <- 0
x_max_3A <- max(c(data_3A$dose, dose_3A_GE80))   # extend to include 80% line
y_min_3A <- 0
y_max_3A <- 100

# Plot the dose response data points
plot(data_3A$dose, data_3A$GE_content,
     xlab = "Dose [mg]",
     ylab = "Gastric retention (%)",
     main = "Dose–response (linear fit)",
     xlim = c(x_min_3A, x_max_3A),
     ylim = c(y_min_3A, y_max_3A))
title("Part A", outer = TRUE, line = -1)

# Show the linear regression curve
abline(m_3A, col = "black", lwd = 2)

# Mark the 50% and 80% target lines
abline(h = 50, col = "red", lty = 2)
abline(h = 80, col = "blue", lty = 2)

# Mark intersection points
points(dose_3A_GE50, 50, pch = 19, col = "red")
points(dose_3A_GE80, 80, pch = 19, col = "blue")

# --- RESIDUAL ERRORS ----------------------------------------------------------

# 1 row, 3 columns
par(mfrow = c(1, 3))

# Residuals vs. Fitted values
# This shows if residuals are centered around 0 and whether variance changes 
# with fitted values (heteroscedasticity). This is the same as Residuals vs. 
# Dose because we only have one predictor (dose), not multiple (e.g., age, sex, 
# etc.)
plot(fitted(m_3A), resid(m_3A), 
     xlab = "Fitted", 
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "red")

# Normal Q–Q plot
qqnorm(resid(m_3A), 
       main = "Normal Q–Q")
qqline(resid(m_3A))

# Histogram of residuals
hist(resid(m_3A), 
     main = "Histogram of residuals", 
     xlab = "Residuals")

# Main title of the three plots
title("Part A Residuals", outer = TRUE, line = -1)

# 1 row, 1 column
par(mfrow = c(1, 1))


# PART B =======================================================================

# --- LOAD DATA & QUICK SUMMARY ------------------------------------------------

# Load data (wide format)
data_3B <- read_csv("Data_T3B.csv", show_col_types = FALSE)

# Convert from wide to long data
data_3B <- gather(data_3B, dose, GE_content, factor_key=TRUE)

# String of dose with removed "GE_content_" at start and "mg" at end
dose_string <- gsub("(^GE_content_|mg$)", "", data_3B$dose)

# Convert to numeric (will give NA for columns like "ID")
dose_numeric <- as.numeric(dose_string)

# Add the numeric dose values to the column
data_3B$dose <- dose_numeric

# Add ID column
data_3B <- data_3B %>%
  mutate(ID = 31:(30 + n()))

# Reorder columns
data_3B <- data_3B[c("ID", "dose", "GE_content")]

# Combine data_3A and data_3B
data_3B <- rbind(data_3A, data_3B)

# Quick summary
summary(data_3B)

# Group by dose and calculate mean of GE_content
data_3B_by_dose <- data_3B %>%
  group_by(dose) %>%
  summarise(mean_GE_content = mean(GE_content))

# --- LINEAR REGRESSION --------------------------------------------------------

# Fitting Linear Models (linear regression)
m_3B <- lm(GE_content ~ dose, data = data_3B)

# Beta values
b0_3B <- coef(m_3B)[1]   # intercept
b1_3B <- coef(m_3B)[2]   # slope

# Solve analytically for the dose that gives 50% and 80% retention
dose_3B_GE50 <- getDoseAtGE(50, b0_3B, b1_3B)
dose_3B_GE80 <- getDoseAtGE(80, b0_3B, b1_3B)

# --- PLOTS --------------------------------------------------------------------

# Histogram of gastric emptying content
hist(data_3B$GE_content,
     main = "Histogram of GE content",
     xlab = "GE content",
     col = "lightblue",
     border = "black")
title("Part B", outer = TRUE, line = -1)

# Boxplot of GE content by dose
boxplot(GE_content ~ dose, data = data_3B,
        main = "Boxplot of GE content by dose [mg]",
        xlab = "Dose [mg]",
        ylab = "GE content",
        col = "lightgreen")
title("Part B", outer = TRUE, line = -1)

# Scatterplot of dose [mg] vs GE content
plot(data_3B$dose, data_3B$GE_content,
     main = "Scatterplot: Dose [mg] vs GE content",
     xlab = "Dose [mg]",
     ylab = "GE content",
     pch = 19, col = "blue")
title("Part B", outer = TRUE, line = -1)

# Line plot of mean GE content
plot(
  data_3B_by_dose$dose,
  data_3B_by_dose$mean_GE_content,
  xlab = "Dose",
  ylab = "Mean GE Content",
  main = "Mean GE Content by Dose",
  type = "b",
  pch = 19
)
title("Part B", outer = TRUE, line = -1)

# Find min/max for axes so lines fit
x_min_3B <- 0
x_max_3B <- 200
y_min_3B <- 0
y_max_3B <- 100

# Plot the dose response data points
plot(data_3B$dose, data_3B$GE_content,
     xlab = "Dose [mg]",
     ylab = "Gastric retention (%)",
     main = "Dose–response (linear fit)",
     xlim = c(x_min_3B, x_max_3B),
     ylim = c(y_min_3B, y_max_3B))
title("Part B", outer = TRUE, line = -1)

# Show the linear regression curves
abline(m_3A, col = "black", lwd = 2)
abline(m_3B, col = "green", lwd = 2)

# Mark the 50% and 80% target lines
abline(h = 50, col = "red", lty = 2)
abline(h = 80, col = "blue", lty = 2)

# Mark intersection points
points(dose_3A_GE50, 50, pch = 19, col = "red")
points(dose_3A_GE80, 80, pch = 19, col = "blue")
points(dose_3B_GE50, 50, pch = 19, col = "red")
points(dose_3B_GE80, 80, pch = 19, col = "blue")

# --- RESIDUAL ERRORS ----------------------------------------------------------

# 1 row, 3 columns
par(mfrow = c(1, 3))

# Residuals vs. Fitted values
# This shows if residuals are centered around 0 and whether variance changes 
# with fitted values (heteroscedasticity). This is the same as Residuals vs. 
# Dose because we only have one predictor (dose), not multiple (e.g., age, sex, 
# etc.)
plot(fitted(m_3B), resid(m_3B), 
     xlab = "Fitted", 
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "red")

# Normal Q–Q plot
qqnorm(resid(m_3B), 
       main = "Normal Q–Q")
qqline(resid(m_3B))

# Histogram of residuals
hist(resid(m_3B), 
     main = "Histogram of residuals", 
     xlab = "Residuals")

# Main title of the three plots
title("Part B Residuals", outer = TRUE, line = -1)

# 1 row, 1 column
par(mfrow = c(1, 1))