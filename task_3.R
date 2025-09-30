# TASK 3 Part2 ---> DRM version (Dose response model)

# The DRC (Dose response curve package) is a ready to use package that includes 
# special functions for dose-response data. Typically used in pharmacololgy!)

#_______________________________________________________________________________________
# Predict the dose that produces 80% gastric retention at 30 minutes.
#_______________________________________________________________________________________

# Install if needed:
#install.packages("drc")

library(drc)

# Load the data
data <- read.csv("Data_T3.csv")

# Showing the data 
head(data)

# FIT THE DATA TO THE DOSE-RESPONSE MODEL

# EXPLANATION OF FUNCTION:
# GE_content ~ Dose --> says that the GE content is dependent on dose
# data = data is just loading our data from Data_T3 to variable data 
# fct = LL.4 is the 4-parameters log-logistic function which is the           !!Im gonna see later on options for this function etc!!
# standard S-shape curve for dose response 

model <- drm(GE_content ~ Dose, data = data, fct = LL.4())

# FINDING THE DOSE FOR 80% RETENTION

# EXPLANATION OF FUNCTION
# ED()            -->  Effective dose function from drc package
# model           -->  We use the model we created 
# type = absolute -->  We want the absolute response level to be 80% not relatively

dose_80 <- ED(model, 80, type = "absolute", interval = "delta")
print(dose_80)

# Extracting data from function so that we can get the 95% CI as well

est <- dose_80[1, "Estimate"]
lwr <- dose_80[1, "Lower"]
upr <- dose_80[1, "Upper"]

# THE ANSWER:
print(paste("Dose for 80% gastric retention:", round(dose_80[1], 2), "mg"))
print(paste("95% CI: [", round(lwr, 2), "-", round(upr, 2), "] mg"))

# Dose for 80% gastric retention: 100.85 mg
# 95% CI: [ 13.83 - 187.86 ] mg              <--- OBS!! Very bad CI. High uncertainty

# PLOTTING THE DRC 
plot(model, type = "all", 
     xlab = "Dose (mg)", 
     ylab = "Gastric Retention (%)",
     main = "Dose-Response Curve")
abline(h = 80, col = "red", lty = 2)
abline(v = dose_80[1], col = "red", lty = 2)



