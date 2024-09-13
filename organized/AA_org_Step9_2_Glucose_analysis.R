library(data.table)
library(dplyr)
library(lme4)
library(lmerTest)
# OS_VERSION <- readLines("/etc/debian_version")
# if (OS_VERSION == "bullseye/sid"){
#   library("lmerTest",lib.loc = "~/Rpackages/")
#   print(OS_VERSION)
# } else {
#   library("lmerTest",lib.loc = "~/Rpackages/4.15/")
#   print(OS_VERSION)
# }
library(nlme)
library(doParallel) 
library(foreach)
library(stringr)


source("~/longitudinalEHR/realData/codes/longehr_pr_direct_realData.R")
time.end = 5

# bootstrap seed
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
seed <- as.numeric(slurm_arrayid)

# variable = "HDL"
variable = "Glucose"

set.seed(seed)
out_file = paste0("~/longitudinalEHR/realData/results/Step9_2_",variable,"/result",seed,".csv")
if (file.exists(out_file)) {
  print("file exists")
  quit()
}

# read the data
load(paste0("~/longitudinalEHR/realData/data/Step9_1_",variable,"_data_combined.RData"))


#### create bootstrap sample ####
long_data_all = data_combined
colnames(long_data_all)[which(colnames(long_data_all)=="Value")] = "Y"
id_crosswalk = data.frame("Deid_ID"=unique(long_data_all$Deid_ID),"id"=1:length(unique(long_data_all$Deid_ID)))
long_data_all = merge(long_data_all,id_crosswalk,by="Deid_ID",all.x=TRUE)
long_data_all = long_data_all[,-c(1)]
## long_data0 should be the bootstrapped data
id_all = unique(long_data_all$id)
# set.seed(123)
# id_sub = sample(id_all,500,replace = FALSE)
id_sub = id_all
boot_length = length(id_sub)
set.seed(seed)
id_boot = sample(id_sub,boot_length,replace = TRUE)
id_boot = id_boot[order(id_boot)]
length(unique(id_boot))
length(id_boot)

# note that after bootstrap, we have overlapped ID, thus, we need to re-label id.
id_select = lapply(id_boot,function(x) which(long_data_all$id == x)) # need to make this faster
id_new = rep(1:boot_length,unlist(lapply(id_select,length)))
long_data0 = long_data_all[unlist(id_select),] 
long_data0$id = id_new
long_data0 = long_data0 %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(time,.by_group = TRUE)
length(unique(long_data0$id))

##### run the analysis #####
VPM_variables = c("Age","Sex","SNP")
LM_fixedEffect_withTime_variables = c("Age","Sex","time","SNP")
LM_fixedEffect_withoutTime_variables = c("Age","Sex","SNP")
LM_randomEffect_variables = "SNP"

#### functions ###### 
regression_fun = function(snp){

  long_data = long_data0[c("id","Y","Age","Sex","time",snp)]
  colnames(long_data)[which(colnames(long_data)==snp)] = "SNP"
  
  long_data$Y = as.numeric(long_data$Y)
  long_data$Age = as.numeric(long_data$Age)/10 # unit: 10 year
  long_data$time = as.numeric(long_data$time)/12 # unit: year
  long_data$Sex = as.numeric(long_data$Sex)-1

  ## round-up the time variable
  long_data$time = round(long_data$time,2)
  long_data$Age = round(long_data$Age,2)
  # remove duplicated records 
  sum(duplicated(long_data[,c("id","time")]))/nrow(long_data)
  long_data = long_data[!duplicated(long_data[,c("id","time")]),]
  
  lmm_long_data = long_data

  # run the mixed-effect model
  print("start LME model")
  lmm_long_data$Y = log(lmm_long_data$Y)
  time.start = Sys.time()
  lmm_standardLME = standard_LME_realData(lmm_long_data)
  time.end = Sys.time()
  lmm_time = difftime(time.end,time.start,units = "secs")
  lmm_standardLME$beta_hat
  lmm_standardLME$beta_sd
  lmm_result = t(data.frame("lmm" = c(lmm_standardLME$beta_hat[5],lmm_standardLME$beta_sd[5])))
  colnames(lmm_result)=c("Estimate","Std. Error")

  print("start conditional LME model")
  time.start = Sys.time()
  beta_hat_conditionalLME = conditional_LME_realData(lmm_long_data)
  time.end = Sys.time()
  conditionalLME_time = difftime(time.end,time.start,units = "secs")
  print(beta_hat_conditionalLME)
  conditionalLME_result = t(data.frame("conditionalLME" = c(beta_hat_conditionalLME[5],0)))
  colnames(conditionalLME_result)=c("Estimate","Std. Error")
  conditionalLME_result 

  print("start joint_LiangGamma model")
  time.start = Sys.time()
  beta_hat_jointLiangGamma_model = joint_LiangGamma_realData(lmm_long_data)
  time.end = Sys.time()
  jointLiangGamma_time = difftime(time.end,time.start,units = "secs")
  beta_hat_jointLiangGamma = beta_hat_jointLiangGamma_model$beta.hat
  gamma_hat_jointLiangGamma = beta_hat_jointLiangGamma_model$gamma.hat
  print(beta_hat_jointLiangGamma)
  jointLiangGamma_result = t(data.frame("jointLiangGamma" = c(beta_hat_jointLiangGamma[3],0)))
  colnames(jointLiangGamma_result)=c("Estimate","Std. Error")

  print("start joint_LY model")
  time.start = Sys.time()
  beta_hat_jointLY = joint_LY_realData(lmm_long_data)$beta_hat_LY
  # beta_hat_joint_LY = rep(0,length(beta_hat_conditionalLME))
  time.end = Sys.time()
  jointLY_time = difftime(time.end,time.start,units = "secs")
  print(beta_hat_jointLY)
  jointLY_result = t(data.frame("jointLY" = c(beta_hat_jointLY[3],0)))
  colnames(jointLY_result)=c("Estimate","Std. Error")

  print("start joint_weightedGEE model")
  time.start = Sys.time()
  beta_hat_weightedGEE = joint_weightedGEE_realData(lmm_long_data)$beta_hat_weightedGEE
  time.end = Sys.time()
  weightedGEE_time = difftime(time.end,time.start,units = "secs")
  print(beta_hat_weightedGEE)
  jointWeightedGEE_result = t(data.frame("jointWeightedGEE" = c(beta_hat_weightedGEE[5],0)))
  colnames(jointWeightedGEE_result)=c("Estimate","Std. Error")


  print("start imputation method")
  time.start = Sys.time()
  beta_hat_imp = imputation_realData_year(lmm_long_data)
  # beta_hat_imp = rep(0,length(beta_hat_conditionalLME))
  time.end = Sys.time()
  imp_time = difftime(time.end,time.start,units = "secs")
  imp_result = t(data.frame("imputation" = c(beta_hat_imp[5],0)))

  lanalysis_result = rbind(lmm_result,conditionalLME_result,jointLiangGamma_result,jointLY_result,jointWeightedGEE_result,imp_result)

  ##############################
  # collapes the data
  print("start the summary statistics approach") 
  collapsed_data0 = na.omit(long_data) 
  collapsed_data1 = collapsed_data0 %>%
                  group_by(id) %>%
                  mutate(min_Y = min(Y),mean_Y = mean(Y),median_Y = median(Y),max_Y = max(Y))
  # data transformation using log
  collapsed_data1$min_Y = log(collapsed_data1$min_Y)
  collapsed_data1$mean_Y = log(collapsed_data1$mean_Y)
  collapsed_data1$median_Y = log(collapsed_data1$median_Y)
  collapsed_data1$max_Y = log(collapsed_data1$max_Y)

  collapsed_data2 = collapsed_data1[,c("id","Age","Sex","SNP","min_Y","mean_Y","median_Y","max_Y")]

  time.start = Sys.time()
  lm_model_min = lm(as.formula(paste0("min_Y~Age+Sex+SNP")),data=collapsed_data2)
  time.end = Sys.time()
  lm_min_time = difftime(time.end,time.start,units = "secs")

  time.start = Sys.time()
  lm_model_mean = lm(as.formula(paste0("mean_Y~Age+Sex+SNP")),data=collapsed_data2)
  time.end = Sys.time()
  lm_mean_time = difftime(time.end,time.start,units = "secs")

  time.start = Sys.time()
  lm_model_median = lm(as.formula(paste0("median_Y~Age+Sex+SNP")),data=collapsed_data2)
  time.end = Sys.time()
  lm_median_time = difftime(time.end,time.start,units = "secs")

  time.start = Sys.time()
  lm_model_max = lm(as.formula(paste0("max_Y~Age+Sex+SNP")),data=collapsed_data2)
  time.end = Sys.time()
  lm_max_time = difftime(time.end,time.start,units = "secs")

  print("start summarizing the results")
  # save the point estimates and the standard errors
  summaryStat_result = t(data.frame("min" = summary(lm_model_min)$coefficients[4,1:2],
                                  "mean" = summary(lm_model_mean)$coefficients[4,1:2],
                                  "median" = summary(lm_model_median)$coefficients[4,1:2],
                                  "max" = summary(lm_model_max)$coefficients[4,1:2]))

  beta_hat_result0 = as.data.frame(rbind(lanalysis_result,summaryStat_result))
  beta_hat_result0$snp = paste0(beta_hat_result0$Estimate," (",beta_hat_result0$`Std. Error`,")")
  beta_hat_result = beta_hat_result0[,"snp",drop=FALSE]
  colnames(beta_hat_result) = snp

  # save the running time
  print("start summarizing the running time")
  (time_result = c(lmm_time,conditionalLME_time,jointLiangGamma_time,jointLY_time,weightedGEE_time,imp_time,
                  lm_min_time,lm_mean_time,lm_median_time,lm_max_time))
  
  results = list("beta_hat_result"=beta_hat_result,"time_result"=time_result)

  return(results)

}

gasparini_fun = function(snp){
    long_data = long_data0[c("id","Y","Age","Sex","time",snp)]
    colnames(long_data)[which(colnames(long_data)==snp)] = "SNP"

    long_data$Y = as.numeric(long_data$Y)
    long_data$Age = as.numeric(long_data$Age)/10 # unit: 10 year
    long_data$time = as.numeric(long_data$time)/12 # unit: year
    long_data$Sex = as.numeric(long_data$Sex)-1

  # round-up the time variable
    long_data$time = round(long_data$time,2)
    long_data$Age = round(long_data$Age,2)
  # remove duplicated records 
    sum(duplicated(long_data[,c("id","time")]))/nrow(long_data)
    long_data = long_data[!duplicated(long_data[,c("id","time")]),]
    
    lmm_long_data = long_data
    lmm_long_data$Y = log(lmm_long_data$Y)
    
    long_data = lmm_long_data
  
  # run the gasparini's method
    print("start gasparini's method")
    time.start = Sys.time()
    gasparini_fit = likelihood_longehr(long_data)
    time.end = Sys.time()
    gasparini_time = difftime(time.end,time.start,units = "secs")
    print(gasparini_fit)
    # gasparini_result = t(data.frame("JMVL-G" = c(gasparini_fit$beta_hat[5],gasparini_fit$beta_sd[5])))
    gasparini_result = paste0(gasparini_fit[5]," (",0,")")

    results = list("beta_hat_result"=gasparini_result,"time_result"=gasparini_time)
    return(results)
}

# # run the regression
# results = NULL
# for(i in 1:length(snps)){
#   (result_i = regression_fun(snp=snps[i]))
#   # extract the number before "(" in the string
#   (result_i_2= as.numeric(gsub("\\(.*","",result_i)))
#   print(i)
#   results = rbind(results,result_i_2)  
# }
# colnames(results) = colnames(result_i)
# rownames(results) = snps
# results

# #####  parallel running ##### 
# registerDoParallel(length(snps))  # Register the cluster
# results2 = foreach(snp=snps,.combine='c') %dopar% regression_fun(snp)
# list_names = c("beta_hat_result","time_result")
# beta_hat_results = do.call(cbind, results2[names(results2) == list_names[1]])
# time_results = do.call(cbind, results2[names(results2) == list_names[2]])

# beta_hat_results = t(apply(beta_hat_results,2,function(x){as.numeric(gsub("\\(.*","",x))}))
# rownames(beta_hat_results) = snps
# colnames(beta_hat_results) = rownames(results2[[1]])
# beta_hat_results

# time_results = t(time_results)
# rownames(time_results) = snps
# colnames(time_results) = rownames(results2[[1]])
# time_results

# outfolder = paste0("~/longitudinalEHR/realData/results/Step9_2_",variable)
# if (!file.exists(outfolder)) {dir.create(outfolder)}
# write.csv(beta_hat_results,paste0("~/longitudinalEHR/realData/results/Step9_2_",variable,"/result",seed,".csv"))
# write.csv(time_results,paste0("~/longitudinalEHR/realData/results/Step9_2_",variable,"/time",seed,".csv"))

##### run the gasparini's method #####

# registerDoParallel(length(snps))  # Register the cluster
# results2 = foreach(snp=snps,.combine='c') %dopar% gasparini_fun(snp)
# list_names = c("beta_hat_result","time_result")
# beta_hat_results = do.call(cbind, results2[names(results2) == list_names[1]])
# time_results = do.call(cbind, results2[names(results2) == list_names[2]])

# beta_hat_results = apply(beta_hat_results,2,function(x){as.numeric(gsub("\\(.*","",x))})
# names(beta_hat_results) = snps
# beta_hat_results

# time_results = as.data.frame(time_results)
# names(time_results) = snps
# time_results

# outfolder = paste0("~/longitudinalEHR/realData/results/Step9_2_",variable)
# if (!file.exists(outfolder)) {dir.create(outfolder)}

# write.csv(beta_hat_results,paste0("~/longitudinalEHR/realData/results/Step9_2_",variable,"/gasparini_result",seed,".csv"))
# write.csv(time_results,paste0("~/longitudinalEHR/realData/results/Step9_2_",variable,"/gasparini_time",seed,".csv"))

##### run the gasparini's method with v_i #####

registerDoParallel(length(snps))  # Register the cluster
# results2 = foreach(snp=snps,.combine='c') %dopar% gasparini_fun(snp)

results2 <- foreach(snp=snps, .combine='c') %dopar% {
  tryCatch({
    gasparini_fun(snp)
  }, error = function(e) {
    paste("Error with SNP:", snp, "-", e$message)
  })
}

error_messages = results2[sapply(results2, is.character)]  # Extract only the character elements
chr_strings = str_extract(error_messages, "chr[0-9XY]+\\.[0-9]+\\.[A-Z]\\.[A-Z]\\.[a-z0-9]+")
(snps_success = setdiff(snps,na.omit(chr_strings)))


list_names = c("beta_hat_result","time_result")
beta_hat_results = do.call(cbind, results2[names(results2) == list_names[1]])
time_results = do.call(cbind, results2[names(results2) == list_names[2]])

beta_hat_results = apply(beta_hat_results,2,function(x){as.numeric(gsub("\\(.*","",x))})
names(beta_hat_results) = snps_success
beta_hat_results

time_results = as.data.frame(time_results)
names(time_results) = snps_success
time_results

print("beta_hat_results")
print(beta_hat_results)

outfolder = paste0("~/longitudinalEHR/realData/results/Step9_2_",variable)
if (!file.exists(outfolder)) {dir.create(outfolder)}

write.csv(beta_hat_results,paste0("~/longitudinalEHR/realData/results/Step9_2_",variable,"/gasparini_withvi_result",seed,".csv"))
write.csv(time_results,paste0("~/longitudinalEHR/realData/results/Step9_2_",variable,"/gasparini_withvi_time",seed,".csv"))
