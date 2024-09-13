rm(list=ls())
load("~/longitudinalEHR/organized/step1_train_data.RData")
source("~/longitudinalEHR/organized/task1_fun.R")

### Take a look at the data
# The data is in long format
head(train_data)
## Missing data:
# Some patients have no measurements at all
tail(train_data)

### Variables in the models ###
outcome_var = "Glucose"
time_var = "time"
id_var = "ID"
VPM_variables = c("Age","Sex","SNP") # in the VP model
LM_fixedEffect_withoutTime_variables = c("Age","Sex","SNP") # in the LP model
LM_fixedEffect_withTime_variables = c(c("Age","Sex","SNP"),time_var) # in the LP model
LM_randomEffect_variables = "SNP" # in the LP model

# All methods: LME, VA-LME, JMVL-LY, JMVL-Liang, JMVL-G, IIRR, imputation_LME
long_data = train_data
colnames(long_data)[which(colnames(long_data)==id_var)] = "id"
colnames(long_data)[which(colnames(long_data)==outcome_var)] = "Y"
colnames(long_data)[which(colnames(long_data)==time_var)] = "time"

#########################
# run the mixed-effect model
print("start LME model")
time.start = Sys.time()
fit_standardLME = standard_LME_fun(long_data=long_data)
time.end = Sys.time()
time_standardLME = difftime(time.end,time.start,units = "secs")
fit_standardLME

print("start VA_LME model")
time.start = Sys.time()
fit_VA_LME = VA_LME_fun(long_data=long_data)
time.end = Sys.time()
time_VA_LME = difftime(time.end,time.start,units = "secs")
fit_VA_LME

print("start JMVL_LY model")
time.start = Sys.time()
fit_JMVL_LY = JMVL_LY_fun(long_data=long_data)
time.end = Sys.time()
time_JMVL_LY = difftime(time.end,time.start,units = "secs")
fit_JMVL_LY

print("start the IIRR-weighting model")
time.start = Sys.time()
fit_IIRR_weighting = IIRR_weighting_fun(long_data=long_data)
time.end = Sys.time()
time_IIRR_weighting = difftime(time.end,time.start,units = "secs")
fit_IIRR_weighting

print("start the JMVL_Liang model")
time.start = Sys.time()
fit_JMVL_Liang = JMVL_Liang_fun(long_data=long_data)
time.end = Sys.time()
time_JMVL_Liang = difftime(time.end,time.start,units = "secs")
fit_JMVL_Liang

print("start the JMVL_G model")
time.start = Sys.time()
fit_JMVL_G = JMVL_G_fun(long_data=long_data)
time.end = Sys.time()
time_JMVL_G = difftime(time.end,time.start,units = "secs")
fit_JMVL_G

print("start the impuation_LME model")
time.start = Sys.time()
fit_imputation_LME = imputation_LME_fun(long_data=long_data,imp_time_factor = 0.5)
fit_imputation_LME = imputation_LME_fun(long_data=long_data)
time.end = Sys.time()
time_imputation_LME = difftime(time.end,time.start,units = "secs")
fit_imputation_LME


