library(ggplot2)
library(lmtest)

data <- read.csv("data/Data_T4.csv")

# ------------------- Visualise data -------------------------------------------
# Number of treated vs control
table(data$TRT)


# Histograms
non_treated = data$LeeIdx[data$TRT == 0]
plot(hist(non_treated), xlim = c(310,350), main = "Histogram of Control Group")

treated = data$LeeIdx[data$TRT == 1]
plot(hist(treated), xlim = c(310,350), main = "Histogram of DiGeMon-123 treated Group")


# Boxplot
boxplot(LeeIdx ~ TRT, data = data, names = c("Control", "Treatment"))


# Scatterplot + correlation for treated group
treated_all = data[data$TRT == 1, c("LeeIdx", "Time")]
plot(treated_all$Time, treated_all$LeeIdx, main = "Scatterplot of Treated Group", pch = 16, xlab = "Time (minutes)", ylab = "Lee Index")
cor(treated_all$Time, treated_all$LeeIdx) # Weak negative correlation -0.4129381

# Scatterplot + correlation for control group
non_treated_all = data[data$TRT == 0, c("LeeIdx", "Time")]
plot(non_treated_all$Time, non_treated_all$LeeIdx, main = "Scatterplot of Control Group", pch = 16, xlab = "Time (minutes)", ylab = "Volume of Gastric Content (%)")
cor(non_treated_all$Time, non_treated_all$LeeIdx)


# ggplot showing LeeIdx vs Time for both groups
# Treat data as categorical instead of numerical for ggplot
data$TRT_cat <- data$TRT
data$TRT_cat <- factor(data$TRT, levels = c(0,1), labels = c("Control", "Treatment"))

ggplot(data, aes(x = Time, y = LeeIdx, color = TRT_cat)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Lee index over time for control and treatment groups", 
       x = "Time (weeks)", y = "Lee index") + 
  theme_minimal()


# --------------------- Linear regression --------------------------------------
# Fit linear regression model with interaction between treatment and time
model <- lm(LeeIdx ~ TRT * Time, data = data)
summary(model)
plot(model)

# p-value for TRTTreatment:Time is 0.012 - Lee index decreases faster for treated rats than control

# Check for heteroscedasticity (also visualised through plot 1 and 3 from plot(model))
bptest(model) # p-value = 0.2644: homoscedasticity - good!

# Check if residuals are normally distributed (also visualised through plot 2 from plot(model))
shapiro.test(resid(model)) # p-value = 0.9098: residuals are normal - good!

# Check for independence of errors
dwtest(model) # p-value = 0.6167: residuals are not correlated - good!


# Present and discuss your findings. 
#- Significant TRT effect → DiGeMon-123 alters Lee index regardless of time.
#- Significant Time effect → Lee index changes over time.
#- Significant TRT:Time interaction → Treatment modifies the rate of change in Lee index over time.

#- Intercept (333.22): This is the estimated Lee index at Time = 0 for the control group (TRT = 0). It serves as the baseline.
#- TRT (-0.79, p = 0.782): The main effect of treatment is not statistically significant. At Time = 0, treated rats do not differ significantly in Lee index from controls.
#- Time (0.24, p = 0.433): Time alone does not significantly affect Lee index in the control group. There's no strong evidence that Lee index changes over time without treatment.
#- TRT:Time (-1.11, p = 0.012): This is the key finding. The interaction between treatment and time is statistically significant. It suggests that Lee index decreases over time in the treated group, at a rate of approximately 1.11 units per week more than the control group.
