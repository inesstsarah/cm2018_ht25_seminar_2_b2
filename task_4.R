
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
plot(treated_all$Time, treated_all$LeeIdx, main = "Scatterplot of Treated Group", pch = 16, xlab = "Time (minutes)", ylab = "Volume of Gastric Content (%)")
cor(treated_all$Time, treated_all$LeeIdx)


non_treated_all = df[df$TRT == 0, c("LeeIdx", "Time")]
plot(non_treated_all$Time, non_treated_all$LeeIdx, main = "Scatterplot of Control Group", pch = 16, xlab = "Time (minutes)", ylab = "Volume of Gastric Content (%)")
cor(non_treated_all$Time, non_treated_all$LeeIdx)

# Find effect of time and TRT to the Lee Idx
model <- lm(LeeIdx ~ TRT, data = df)
summary(model)$coef
plot(model)

model <- lm(LeeIdx ~ Time, data = df)
summary(model)$coef
plot(model)

model <- lm(LeeIdx ~ TRT*Time, data = df)
summary(model)$coef
plot(model)

model

# Present and discuss your findings. 
#- Significant TRT effect → DiGeMon-123 alters Lee index regardless of time.
#- Significant Time effect → Lee index changes over time.
#- Significant TRT:Time interaction → Treatment modifies the rate of change in Lee index over time.

#- Intercept (333.22): This is the estimated Lee index at Time = 0 for the control group (TRT = 0). It serves as the baseline.
#- TRT (-0.89, p = 0.755): The main effect of treatment is not statistically significant. At Time = 0, treated rats do not differ significantly in Lee index from controls.
#- Time (0.24, p = 0.432): Time alone does not significantly affect Lee index in the control group. There's no strong evidence that Lee index changes over time without treatment.
#- TRT:Time (-1.13, p = 0.011): This is the key finding. The interaction between treatment and time is statistically significant. It suggests that Lee index decreases over time in the treated group, at a rate of approximately 1.13 units per week more than the control group.
