library(ggplot2)
library(lmtest)

data <- read.csv("data/Data_T4.csv")

# Treat data as categorical instead of numerical
data$TRT <- factor(data$TRT, levels = c(0,1), labels = c("Control", "Treatment"))

# Visualise data
boxplot(LeeIdx ~ TRT, data = data, names = c("Control", "Treatment"))
table(data$TRT)

ggplot(data, aes(x = Time, y = LeeIdx, color = TRT)) + 
  geom_point() +
  geom_smooth(method = "lm") + # For checking linearity
  labs(title = "Lee index over time for control and treatment groups", x = "Time (weeks)", y = "Lee index") + 
  theme_minimal()

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



# Read the data
df = read.csv("data/Data_T4.csv", sep = ",")


df["LeeIdx"] = df[,2]
df["TRT"] = df[,3]
df["Time"] = df[,4]
df = df[,5:7]
df

# Investigate if the GLP-1 agonist has any effect on the Lee Index
# See histograms of the DiGeMon levels  
df
non_treated = df$LeeIdx[df$TRT == 0]

plot(hist(non_treated), xlim = c(310,350), main = "Histogram of Control Group")

treated = df$LeeIdx[df$TRT == 1]
plot(hist(treated_TRT), xlim = c(310,350), main = "Histogram of DiGeMon-123 treated Group")

# Treatment time vs Lee Index in the treated
# and non-treated groups

df
treated_all = df[df$TRT == 1, c("LeeIdx", "Time")]
plot(treated_all$Time, treated_all$LeeIdx, main = "Scatterplot of Treated Group", pch = 16, xlab = "Time (minutes)", ylab = "Lee Index")
cor(treated_all$Time, treated_all$LeeIdx)
# Weak negative correlation -0.4129381 could be relevant bc the y is

non_treated_all = df[df$TRT == 0, c("LeeIdx", "Time")]
plot(non_treated_all$Time, non_treated_all$LeeIdx, main = "Scatterplot of Control Group", pch = 16, xlab = "Time (minutes)", ylab = "Volume of Gastric Content (%)")
cor(non_treated_all$Time, non_treated_all$LeeIdx)

# Add categorical TRT data to the correlation 


# Find effect of time and TRT to the Lee Idx
model <- lm(LeeIdx ~ TRT, data = df)
summary(model)$coef
plot(model)

model <- lm(LeeIdx ~ Time, data = df)
summary(model)$coef
plot(model)

model <- lm(LeeIdx ~ TRT+Time, data = df) # Difference between these types of modelling
summary(model)$coef
plot(model)

model <- lm(LeeIdx ~ TRT*Time, data = df)
summary(model)$coef
plot(model)


# Present and discuss your findings. 
#- Significant TRT effect → DiGeMon-123 alters Lee index regardless of time.
#- Significant Time effect → Lee index changes over time.
#- Significant TRT:Time interaction → Treatment modifies the rate of change in Lee index over time.

#- Intercept (333.22): This is the estimated Lee index at Time = 0 for the control group (TRT = 0). It serves as the baseline.
#- TRT (-0.89, p = 0.755): The main effect of treatment is not statistically significant. At Time = 0, treated rats do not differ significantly in Lee index from controls.
#- Time (0.24, p = 0.432): Time alone does not significantly affect Lee index in the control group. There's no strong evidence that Lee index changes over time without treatment.
#- TRT:Time (-1.13, p = 0.011): This is the key finding. The interaction between treatment and time is statistically significant. It suggests that Lee index decreases over time in the treated group, at a rate of approximately 1.13 units per week more than the control group.

Estimate Std. Error     t value     Pr(>t) 
(Intercept) 333.2154618  1.9145668 174.0422197 7.802612e-92 TRT          -0.8945587  2.8565897  -0.3131562 7.551205e-01 Time          0.2428079  0.3071695   0.7904688 4.320022e-01 TRT:Time     -1.1274633  0.4305956  -2.6183807 1.088228e-02
