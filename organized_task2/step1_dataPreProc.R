library(dplyr)
library(data.table)


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

# take a subset of the data
id_all = unique(data7$id) # length=10364
set.seed(123)
id_sub = sample(id_all,2000,replace = FALSE)
data_sub = data7[which(data7$id %in% id_sub),]

id_crosswalk_sub = data.frame("id"=unique(data_sub$id),"id_sub"=1:length(unique(data_sub$id)))
data_sub$id = id_crosswalk_sub$id_sub[match(data_sub$id,id_crosswalk_sub$id)]
head(data_sub)

# remove observations with time==0
data_sub1 = data_sub[which(data_sub$time>0),]

## round-up the time variable
data_full = data_sub1
data_full$time = round(data_full$time,2)
data_full$Age = round(data_full$Age,2)
# remove duplicated records 
data_full = data_full[!duplicated(data_full[,c("id","time")]),]

long_data = data_full[,c("id","Age","Sex","SNP","time","Y")]
head(long_data)

surv_data = data_full[,c("id","Age","Sex","SNP","D","d")]
surv_data = surv_data[!duplicated(surv_data),]
mean(surv_data$d,na.omit=TRUE)
head(surv_data)

length(unique(long_data$id))
length(unique(surv_data$id))

# save data
save(long_data,surv_data,file="~/longitudinalEHR/organized_task2/train_data_task2.RData")

