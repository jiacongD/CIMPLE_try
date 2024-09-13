rm(list=ls())
load("~/longitudinalEHR/organized_task2/train_data_task2.RData")
source("~/longitudinalEHR/organized_task2/task2_fun.R")
ls()

# define the variables
SM_variables = c("Age","Sex","SNP","Y")
SM_base_variables = c("Age","Sex","SNP")
LM_fixedEffect_withTime_variables = c("Age","Sex","SNP","time")
LM_fixedEffect_withoutTime_variables = c("Age","Sex","SNP")
LM_randomEffect_variables = c("SNP")

##### Run the analysis #####
# All methods: Cox_fun, Imputation_Cox, JMLD, VAImputation_Cox, VA_JMLD 
print("start the Cox model")
time_start = Sys.time()
alpha_hat_coxph = cox_fun(long_data, surv_data)
time_end = Sys.time()
time_coxph = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_coxph) = "coxph"
print(alpha_hat_coxph)

print("start JMLD")
time_start = Sys.time()
alpha_hat_JMLD = JMLD_fun(long_data, surv_data)
time_end = Sys.time()
time_JMLD = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_JMLD) = "JMLD"
print(alpha_hat_JMLD$surv_proc)

print("start Imputation_Cox")
time_start = Sys.time()
alpha_hat_imp_cox = Imputation_Cox_fun(long_data, surv_data,imp_time_factor = 0.5)
time_end = Sys.time()
time_imp_cox = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_imp_cox) = "Imputation_Cox"
print(alpha_hat_imp_cox)

print("start VAImputation_Cox")
time_start = Sys.time()
alpha_hat_VAimp_cox = VAImputation_Cox_fun(long_data, surv_data,imp_time_factor = 0.5)
time_end = Sys.time()
time_VAimp_cox = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_VAimp_cox) = "VAImputation_Cox"
print(alpha_hat_VAimp_cox)

print("start VA_JMLD")
time_start = Sys.time()
alpha_hat_VAJMLD = VA_JMLD_fun(long_data, surv_data)
time_end = Sys.time()
time_VAJMLD = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_VAJMLD) = "VA_JMLD"
print(alpha_hat_VAJMLD$surv_proc)
