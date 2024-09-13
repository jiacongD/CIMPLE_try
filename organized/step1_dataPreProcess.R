rm(list=ls())
set.seed(2024)
library(dplyr)
library(tidyr)
load("~/longitudinalEHR/organized/AA_org_Step9_1_Glucose_data_combined.RData")
train_data = data_combined[,c(1:5,7)]
colnames(train_data) = c("ID","Glucose","time","Age","Sex","SNP")
apply(train_data, 2, function(x) mean(is.na(x)))

# Take a subset
id_all = unique(train_data$ID)
id_subset = sample(id_all, 2000, replace = FALSE)
id_crosswalk = data.frame(ID = id_subset, id_new = 1:length(unique(id_subset)))
train_data = merge(train_data, id_crosswalk, by = "ID")
head(train_data)

train_data_subset = train_data[train_data$ID %in% id_subset,c("id_new","Glucose","time","Age","Sex","SNP")]
colnames(train_data_subset)[1] = "ID"
train_data = train_data_subset
train_data = train_data %>% arrange(ID, time)

train_data_org = train_data
train_data = as.data.frame(sapply(train_data_org,as.numeric))
train_data$Glucose = log(train_data$Glucose)
train_data$Age = train_data$Age/10
train_data$time = round(train_data$time/12,2)
train_data$Sex = train_data$Sex-1
# remove duplicated rows
head(train_data)
max(na.omit(train_data$time))

train_data = train_data[!duplicated(train_data[,c("ID","time")]),]

# remove if the time is greater than 60 and time is 0, keep NA
train_data = train_data %>% 
                filter((time<=60)%>% replace_na(TRUE)) %>% 
                filter((time!=0)%>% replace_na(TRUE))

# save the data
save(train_data, file = "~/longitudinalEHR/organized/step1_train_data.RData")
