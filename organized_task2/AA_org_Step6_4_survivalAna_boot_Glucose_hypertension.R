library(data.table)
library(dplyr)
library(MASS)
library(dplyr)
library(lme4)
library(survival)
library("JM",lib.loc = "~/Rpackages/")
OS_VERSION <- readLines("/etc/debian_version")
if (OS_VERSION == "bullseye/sid"){
  library("JMbayes2",lib.loc = "~/Rpackages/")
  print(OS_VERSION)
} else {
  library("JMbayes2",lib.loc = "~/Rpackages/4.15/")
  print(OS_VERSION)
}
library(nleqslv)
library(emdbook)
library(mice)
source("/net/snowwhite/home/jiacong/longitudinalEHR/realData_Aim3/association_YD/codes/estimation_fun.R")



###
year_thres = 5
variable = "Glucose"

############## survival data ################
ExPRS_results2 = read.csv("/net/snowwhite/home/jiacong/longitudinalEHR/realData_Aim3/data/ExPRS_results2.csv")
phecodes = as.character(ExPRS_results2$phecode[which(ExPRS_results2$exposure==variable)])
phecodes
outcome_data_org = fread("/net/snowwhite/home/jiacong/longitudinalEHR/magic_data/Outcomes/HPI_4347_20200219_Phecodes.txt")

outcome = "401.1" # Essential hypertension
print(outcome)
surv_data = outcome_data_org[which(outcome_data_org$phecode==outcome),]
colnames(surv_data) = c("Deid_ID","phecode","DaysSinceBirth_phecode")
dim(surv_data)

(output_folder = paste0("/net/snowwhite/home/jiacong/longitudinalEHR/realData_Aim3/association_YD/results/Step6_4_",variable,"_",outcome))
if(!dir.exists(output_folder)){dir.create(output_folder)}

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
seed <- as.numeric(slurm_arrayid)
# if the following file exists, then end the program
if(file.exists(paste0(output_folder,"/alpha_hat_summaryStat",seed,".csv"))){
  print("file exists")
  quit()
}

############## longitudinal data ################
# get valid id
EHR_time_data = fread("/net/snowwhite/home/jiacong/longitudinalEHR/magic_data/HPI_4347_20200219_DaysSinceBirth_FirstLast_EHR_entries.txt", header = TRUE, sep = "\t")
sum(EHR_time_data$DaysInEHR<0) # 0 negative values, two entries
id_valid1 = EHR_time_data$Deid_ID[which(EHR_time_data$YearsInEHR>year_thres)]
EHR_time_data = EHR_time_data[which(EHR_time_data$YearsInEHR>year_thres),]
EHR_time_data = EHR_time_data[,c("Deid_ID","FirstDaySinceBirth")]

# Demographic data
file_loc = "/net/snowwhite/home/jiacong/longitudinalEHR/magic_data/"
Demo_data = fread(paste0(file_loc,"HPI_4347_20200219_Demographics.txt"), header = TRUE, sep = "\t")
Demo_data = Demo_data[which(Demo_data$Deid_ID %in% id_valid1),]
Demo_data1 = Demo_data[,c("Deid_ID","Age","Sex","RaceEthnicity4")]
Demo_data1$Sex = ifelse(Demo_data1$Sex=="M",1,ifelse(Demo_data1$Sex=="F",0,NA))
Demo_data1$Race = ifelse(Demo_data1$RaceEthnicity4=="Caucasian / Non-Hispanic","White",
                         ifelse(Demo_data1$RaceEthnicity4=="African American / Non-Hispanic","Black","Other"))
Demo_data2 = na.omit(Demo_data1[,c("Deid_ID","Sex","Race")])

# genetic data
geneticData_org = read.csv("~/longitudinalEHR/realData_Aim3/data/HPI-9573_Freeze6_LabWAS_SNP_genotypes_reformat.csv")
colnames(geneticData_org)[1] = "Deid_ID"
# requestedFile = read.csv("~/longitudinalEHR/realData_Aim3/data/LabWAS_findings2.csv")
# crosswalkFile = read.csv("~/longitudinalEHR/realData_Aim3/data/HPI-9573_requested_variants_hg19_to_hg38.csv")
# requestedFile$Freeze6_variantID=crosswalkFile$Freeze6_variantID[match(requestedFile$Variant,crosswalkFile$RequestedVariant)]
if(variable == "HDL"){
  geneticData = geneticData_org[,c(1,which(colnames(geneticData_org)=="chr9.104903458.G.A.rs2575876"))]
}
if(variable == "Glucose"){
  geneticData = geneticData_org[,c(1,which(colnames(geneticData_org)=="chr10.112998590.C.T.rs7903146"))]
}
colnames(geneticData)=c("Deid_ID","SNP")

# lab data
file_loc2 = "~/MGI_data/Exposure/"
if(variable=="Glucose"){lab_file_loc = paste0(file_loc2,"HPI_4347_20200219_GLUC_LOINC_2345-7.txt")}
if(variable=="HDL"){lab_file_loc = paste0(file_loc2,"HPI_4347_20200219_HDL_LOINC_2085-9.txt")}

# read data
lab_data = fread(lab_file_loc, header = TRUE, sep = "\t")    
if("VALUE" %in% colnames(lab_data)) colnames(lab_data)[colnames(lab_data)=="VALUE"] = "Value"
if("DeID_PatientID" %in% colnames(lab_data)) colnames(lab_data)[colnames(lab_data)=="DeID_PatientID"] = "Deid_ID"

lab_data = lab_data[which(lab_data$Deid_ID %in% id_valid1),] 
lab_data = lab_data[,c("Deid_ID","Value","DaysSinceBirth")]
lab_data$Value = as.numeric(lab_data$Value)
lab_data = na.omit(lab_data) %>%
           arrange(Deid_ID,DaysSinceBirth)
apply(is.na(lab_data),2,sum)
length(unique(lab_data$Deid_ID))

# keep the longitudinal measurements within 5 years after the first visit
# we filter out many patients, because they either do not have measurement, or have measurements after 5 years
data0 = merge(lab_data,EHR_time_data[,c("Deid_ID","FirstDaySinceBirth")],by="Deid_ID", all.x = TRUE)
data1 = data0 %>%
        mutate(CTime = FirstDaySinceBirth+365*5) %>%
        mutate(R = ifelse(is.na(Value),1,ifelse(DaysSinceBirth<=CTime,1,0))) %>%
        filter(R==1) %>%
        dplyr::select(-c("R"))
apply(is.na(data1),2,sum)
length(unique(data1$Deid_ID))
dim(data1)

# merge the survival information
data2 = merge(data1, surv_data, by="Deid_ID",all.x=TRUE)
length(unique(data2$Deid_ID))
dim(data2)

# disease onset time needs to be within 5 years after the first visit
data3 = data2 %>%
        mutate(DaysSinceBirth_phecode1=ifelse(is.na(DaysSinceBirth_phecode),NA,ifelse(DaysSinceBirth_phecode>CTime,NA,DaysSinceBirth_phecode))) %>% 
        mutate(S = ifelse(DaysSinceBirth<=DaysSinceBirth_phecode1|is.na(DaysSinceBirth_phecode1),1,0)) %>%
        filter(S==1) %>% # only leave records before the disease onset and those who didn't have the disease
        mutate(time = DaysSinceBirth-FirstDaySinceBirth) %>%
        filter(time>=0) %>%
        mutate(Age = FirstDaySinceBirth/365) %>%
        mutate(D=ifelse(is.na(DaysSinceBirth_phecode1), 365*5, DaysSinceBirth_phecode1-FirstDaySinceBirth)) %>%
        mutate(d=ifelse(is.na(DaysSinceBirth_phecode1), 0, 1)) %>%
        dplyr::select(Deid_ID,time,Age,Value,D,d)        
apply(is.na(data3),2,sum)
max(data3$D,na.rm=TRUE)/365

# merge the demographics
data4 = merge(data3,Demo_data2,by="Deid_ID",all.x=TRUE)
data4 = data4 %>%
        filter(Race=="White") %>%
        dplyr::select(Deid_ID,time,Age,Value,Sex,D,d)
apply(data4,2,function(x) sum(is.na(x)))
head(data4)

# merge the genetic data
data5 = merge(data4,geneticData,by="Deid_ID",all.x=TRUE)
# remove na's in genetic information
data5 = data5[complete.cases(data5$SNP), ]
apply(is.na(data5),2,mean)
dim(data5)
head(data5)

### Re-format the data and units ###
data5$time = as.numeric(data5$time/365) # unit: 1 years
data5$Age = as.numeric(data5$Age/10) # unit: 10 years
data5$Sex = as.numeric(data5$Sex)
data5$D = as.numeric(data5$D/365) # unit: 1 years
data5$d = as.numeric(data5$d)
data6 = data5[,c("Deid_ID","Age","Sex","SNP","time","Value","D","d")]
data6 = na.omit(data6)
head(data6)

# re-number ID
id_crosswalk = data.frame("Deid_ID"=unique(data6$Deid_ID),"id"=1:length(unique(data6$Deid_ID)))
data6$id = id_crosswalk$id[match(data6$Deid_ID,id_crosswalk$Deid_ID)]
data7 = data6[,c("id","Age","Sex","SNP","time","Value","D","d")]
colnames(data7)[which(colnames(data7)=="Value")] = "Y"
# data7$Y = qnorm((rank(data7$Y,na.last = "keep")-0.5)/sum(!is.na(data7$Y)))
data7$Y = log(data7$Y) # log-transform

############ prepare the dataset ##############
# data_full = data7 # on the full data

########## bootstrap ###########
id_all = unique(data7$id) # length=10364
set.seed(123)
# id_sub = sample(id_all,500,replace = FALSE)
id_sub = id_all
boot_length = length(id_sub)

# bootstrap
rep_boot = 300

set.seed(seed)

# bootstrap sample
id_boot = sample(id_sub,boot_length,replace = TRUE)
# note that after bootstrap, we have overlapped ID, thus, we need to re-label id.
id_select = lapply(id_boot,function(x) which(data7$id == x))
id_new = rep(1:boot_length,unlist(lapply(id_select,length)))

data_boot = data7[unlist(id_select),] 
data_boot$id = id_new
data_boot = data_boot %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(time,.by_group = TRUE)

data_full = data_boot
apply(data_full,2,class)

############################
## round-up the time variable
data_full$time = round(data_full$time,2)
data_full$Age = round(data_full$Age,2)
# remove duplicated records 
sum(duplicated(data_full[,c("id","time")]))/nrow(data_full)
data_full = data_full[!duplicated(data_full[,c("id","time")]),]

long_data = data_full[,c("id","Age","Sex","SNP","time","Y")]
head(long_data)

surv_data = data_full[,c("id","Age","Sex","SNP","D","d")]
surv_data = surv_data[!duplicated(surv_data),]
mean(surv_data$d,na.omit=TRUE)
head(surv_data)

#################### run the analysis #####################
SM_variables = c("Age","Sex","SNP","Y")
SM_base_variables = c("Age","Sex","SNP")
LM_fixedEffect_withTime_variables = c("Age","Sex","SNP","time")
LM_fixedEffect_withoutTime_variables = c("Age","Sex","SNP")
LM_randomEffect_variables = c("SNP")

# long_data_copy = long_data
# surv_data_copy = surv_data

###### run the analysis #####
print("start coxph")
time_start = Sys.time()
alpha_hat_coxph = coxph_fun_realData(long_data, surv_data)
time_end = Sys.time()
time_coxph = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_coxph) = "coxph"
print(alpha_hat_coxph)
write.csv(alpha_hat_coxph, paste0(output_folder,"/alpha_hat_coxph",seed,".csv"))

print("start JM")
time_start = Sys.time()
alpha_hat_joint = JM_fun_realData(long_data, surv_data)
time_end = Sys.time()
time_JM = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_JM) = "JM"
print(alpha_hat_joint)
write.csv(alpha_hat_joint$surv_proc, paste0(output_folder,"/alpha_hat_joint",seed,".csv"))
write.csv(alpha_hat_joint$long_proc, paste0(output_folder,"/beta_hat_joint",seed,".csv"))

print("start imputation")
time_start = Sys.time()
alpha_hat_imp = imputation_realData_year(long_data,surv_data)
time_end = Sys.time()
time_imp = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_imp) = "imputation"
print(alpha_hat_imp)
write.csv(alpha_hat_imp, paste0(output_folder,"/alpha_hat_imp",seed,".csv"))

print("start the summary statistics approach")
summaryStat_fit = summaryStat(long_data, surv_data)
print(summaryStat_fit$alpha_hat)
print(summaryStat_fit$time)
write.csv(summaryStat_fit$alpha_hat, paste0(output_folder,"/alpha_hat_summaryStat",seed,".csv"))

print("start the VA-imputation method")
time_start = Sys.time()
alpha_hat_VAimp = VAimputation_realData_year(long_data,surv_data)
time_end = Sys.time()
time_VAimp = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_VAimp) = "VAimp"
print(time_VAimp)
write.csv(alpha_hat_VAimp, paste0(output_folder,"/alpha_hat_VA",seed,".csv"))

print("start the VAJM method")
time_start = Sys.time()
alpha_hat_VAJM = VAJM_fun_realData(long_data, surv_data)
time_end = Sys.time()
time_VAJM = as.numeric(difftime(time_end, time_start, units = "secs")); names(time_VAJM) = "VAJM"
print(time_VAJM)
write.csv(alpha_hat_VAJM$surv_proc, paste0(output_folder,"/alpha_hat_VAJM",seed,".csv"))
write.csv(alpha_hat_VAJM$long_proc, paste0(output_folder,"/beta_hat_VAJM",seed,".csv"))

# combine all time
time_all = c(time_coxph,time_JM,time_imp,summaryStat_fit$time,time_VAimp,time_VAJM)
write.csv(time_all, paste0(output_folder,"/time_all",seed,".csv"))
