# ==============================================================================

# CM2018 – Task 3
# Dose–response at 30 min


# === 0) SETUP =================================================================

# Install Packages (install once if needed)
# install.packages(c("readr","dplyr","ggplot2","ggtext","drc", "DescTools"))
library(readr);
library(dplyr);
library(ggplot2);


# === 1) LOAD DATA & QUICK SUMMARY =============================================

# Set working directory
#setwd("C:/OneDrive/KTH/Programs/Civil Engineer in Medical Technology 300 credits/
#      Courses/CM2018 Statistics for Medical Engineering 7.5 credits/HT25-1/
#      Seminars/Seminar 2/cm2018_ht25_seminar_2_b2/Data")

# # Load data (n=30)
data <- read_csv("Data_T3.csv", show_col_types = FALSE, n_max = 30)

# Rename columns
names(data) <- c("ID", "dose_mg", "GE_content")
dose_mg <- data$dose_mg
GE_content = data$GE_content

# Quick summary
summary(data)

# Group by dose_mg and calculate mean of GE_content
data_by_dose <- data %>% 
  group_by(dose_mg) %>%
  summarise(mean_GE_content = mean(GE_content))

# Show summary of data grouped by dose
summary(data_by_dose)


# === 2) PLOT ==================================================================

# Histogram of gastric emptying content
hist(GE_content,
     main = "Histogram of GE content",
     xlab = "GE content",
     col = "lightblue",
     border = "black")

# Boxplot of GE content by dose
boxplot(GE_content ~ dose_mg, data = data,
        main = "Boxplot of GE content by dose [mg]",
        xlab = "Dose [mg]",
        ylab = "GE content",
        col = "lightgreen")

# Scatterplot of dose [mg] vs GE content
plot(dose_mg, GE_content,
     main = "Scatterplot: Dose [mg] vs GE content",
     xlab = "dose_mg",
     ylab = "GE_content",
     pch = 19, col = "blue")

# Same scatterplot using ggplot (nicer style)
ggplot(data, aes(x = dose_mg, y = GE_content)) +
  geom_point(color = "blue") +
  labs(title = "ggplot: Dose [mg] vs GE content",
       x = "Dose [mg]", y = "GE content")

# Simple barplot of mean values
barplot(data_by_dose$mean_GE_content, names.arg = data_by_dose$dose_mg,
        main = "Mean GE content by dose",
        xlab = "Dose [mg]", ylab = "Mean GE content",
        col = "orange")



# === 3) LINEAR REGRESSION =====================================================

# Fitting Linear Models (linear regression)
#m <- lm(GE_content ~ dose_mg, data = data)

# Fit a simple linear model
m <- lm(GE_content ~ dose_mg, data = data)

# Generalized least squares for unequal variance
#m <- gls(GE_content ~ dose_mg, data = data)

summary(m)

# Solve analytically for the dose that gives 50% and 80% retention
b0 <- coef(m)[1]   # intercept
b1 <- coef(m)[2]   # slope

getDose <- function(target_GE) {
  (target_GE - b0) / b1
}

dose_50 <- getDose(50)
dose_80 <- getDose(80)

dose_50
dose_80

# Plot data + fitted line
plot(dose_mg, GE_content,
     xlab = "Dose (mg)",
     ylab = "Gastric retention (%)",
     main = "Dose–response (linear fit)")
abline(m, col = "red", lwd = 2)

# Mark the 50% and 80% target lines
abline(h = 50, col = "blue", lty = 2)
abline(h = 80, col = "green", lty = 2)

# ANOVA
#anova <- aov(GE_content ~ dose_mg, data = data)
#summary(anova)

#PostHocTest(anova, method ="lsd")
