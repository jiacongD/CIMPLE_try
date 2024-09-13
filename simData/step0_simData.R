
### Task 1: data simualation ###
library(lme4)
library(nleqslv)
library(mice)
library(MASS)
library(dplyr)

source("~/longitudinalEHR/codes2/association_AY/simData.R")
source("~/longitudinalEHR/codes2/association_AY/simData_fun.R")

seed = 2
set.seed(seed)
m=1000
time.start=1
time.end = 60

simData = simData_withIVPunmeasuredGamma_randomZ(m,time.start,time.end,seed)

train_data = simData$long_data
save(train_data,file="~/longitudinalEHR/organized/step1_sim_train_data.RData")


### Task 2: data simulation with unmeasured confounder ###
rm(list=ls())
source("~/longitudinalEHR/codes2/association_YD/simData.R")
source("~/longitudinalEHR/codes2/association_YD/simData_fun.R")

m=1000
time.start=1
time.end=60
seed = 2
set.seed(seed)

data_gene = simData_withIVPunmeasuredGamma_randomZ(m,time.start,time.end,seed)
long_data = data_gene$long_data
surv_data = data_gene$surv_data

dim(long_data)
dim(surv_data)
save(long_data,surv_data,file="~/longitudinalEHR/organized_task2/step1_sim_train_data.RData")

