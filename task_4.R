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

