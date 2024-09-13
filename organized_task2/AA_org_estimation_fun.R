
coxph_fun_realData = function(long_data, surv_data){
  
  long_data2 = long_data %>%
    # mutate(time=time*10) %>%
    group_by(id) %>%
    mutate(time0 = ifelse(is.na(lag(time)),0,lag(time) )) %>%
    filter(time!=0) %>%
    filter(time0<=time)
  long_data2$d0 = surv_data$d[match(long_data2$id,surv_data$id)]
  long_data2 = na.omit(long_data2) %>%
               group_by(id) %>%
               mutate(d=ifelse(time<max(time,na.rm = TRUE),0,d0))
  model_formula = as.formula(paste("Surv(time0, time, d) ~ ",paste(SM_variables,collapse="+")))
  model = coxph(model_formula, data = long_data2) 
  (alpha.hat = summary(model)$coef[,1])
  
  return(alpha.hat)
}

JM_fun_realData = function(long_data, surv_data){

  # remove patients with no Y measurement
  zeroRecords_id = long_data[is.na(long_data$Y),"id"]
  long_data = long_data[!long_data$id %in% zeroRecords_id,]
  surv_data = surv_data[!surv_data$id %in% zeroRecords_id,]  

  # longitudinal submodel, try nlminb optimizer first, if there is an error, then use optim optimizer
  lmeFit_fixedformula = as.formula(paste("Y ~ ",paste(LM_fixedEffect_withTime_variables,collapse="+")))
  lmeFit_randomformula = as.formula(paste("~",paste(LM_randomEffect_variables,collapse="+"),"|id"))
  control_optim = lmeControl(msMaxIter = 1000, msMaxEval = 1000, opt='optim')
  control_nlminb = lmeControl(msMaxIter = 1000, msMaxEval = 1000, opt='nlminb')
  lmeFit = tryCatch({
    lmeFit = lme(lmeFit_fixedformula,random = lmeFit_randomformula, data = long_data, control=control_optim)
  }, error = function(e) {
    print(paste0("Error with optim: ",e))
    lmeFit = lme(lmeFit_fixedformula,random = lmeFit_randomformula, data = long_data, control=control_nlminb)
  })
  print(lmeFit)

  # survival submodel
  coxFit_formula = as.formula(paste("Surv(D, d) ~ ",paste(SM_base_variables,collapse="+")))
  coxFit = coxph(coxFit_formula, data = surv_data, x = TRUE)
  summary(coxFit)

  # jointFit = JMbayes2::jm(coxFit,lmeFit,time_var="time", n_chains=1L,n_iter=11000L,n_burnin=1000L) 
  jointFit = JMbayes2::jm(coxFit,lmeFit,time_var="time", n_chains=1L) 

  print(summary(jointFit))
  
  surv_proc = unlist(coef(jointFit))
  long_proc = unlist(fixef(jointFit))

  results = list("long_proc" = long_proc,
                 "surv_proc" = surv_proc)

  return(results)
}

imputation_realData_year = function(long_data,surv_data){
    
  time_unit = year_unit=1
  month_unit=12/6
  data = long_data
  # monthly-grouped data
  data = data %>%
      mutate(time_6month = time*time_unit*2,time_roundup_6month = ceiling(time_6month))%>%
      group_by(id,time_roundup_6month) %>%
      mutate(Y_mean = mean(Y,na.rm = TRUE),Age_mean = mean(Age,na.rm=TRUE)) %>%
      dplyr::select(id,Y_mean,Age_mean,Sex,SNP,time_roundup_6month)
  data1 = data[!duplicated(data),]
  colnames(data1) = c("id","Y","Age","Sex","SNP","time")

  # create a new dataset with the same columns in data1. 
  # For a given id, the time is from 1 to max_time. 
  # If there is no value of Y at a given time, then Y is NA.
  max_time <- max(data1$time, na.rm = TRUE)
  id_all = rep(unique(data1$id),each = max_time)
  time_all = rep(as.numeric(1:max_time),length(unique(data1$id)))
  Age_all = rep(unique(data1[,c("id","Age")])$Age,each=max_time)
  Sex_all = rep(unique(data1[,c("id","Sex")])$Sex,each=max_time)
  SNP_all = rep(unique(data1[,c("id","SNP")])$SNP,each=max_time)
  data2 = data.frame(id = id_all,time = time_all,Sex = Sex_all, Age = Age_all,SNP = SNP_all)
  # delete time after the disease time
  data2$D = surv_data$D[match(data2$id,surv_data$id)]
  head(data2)
  
  sum(data2$time<=ceiling(data2$D*time_unit*month_unit))
  ######
  data2 = data2[data2$time<=ceiling(data2$D*time_unit*month_unit),c("id","time","Age","Sex","SNP")]
  
  # join data1 to data2 by id and time
  data3 = left_join(data2,data1,by = c("id","time","Age","Sex","SNP"))
  data3$time = data3$time/2/time_unit # convert time back to the original scale

  df_full2 = data3
  dim(df_full2)

  # empty imputation
  imp0 = mice(as.matrix(df_full2),maxit=0) 
  predM = imp0$predictorMatrix 
  impM = imp0$method

  #specify predictor matrix and method
  predM1 = predM 
  predM1["Y","id"] = -2 
  predM1["Y",LM_fixedEffect_withTime_variables] = 1 #fixed x effects imputation 
  impM1 = impM 
  impM1["Y"]="2l.lmer"

  # multilevel imputation 
  imp1=mice(as.matrix(df_full2), m=5, 
          predictorMatrix=predM1, method=impM1, maxit=5 )

  # # cox-ph model
  # long_data_imp = complete(imp1)
  # long_data_imp2 = long_data_imp %>%
  #   group_by(id) %>%
  #   mutate(time0 = ifelse(is.na(lag(time)),0,lag(time) ))
  # long_data_imp2$d0 = surv_data$d[match(long_data_imp2$id,surv_data$id)]
  # long_data_imp2 = na.omit(long_data_imp2) %>%
  #              group_by(id) %>%
  #              mutate(d=ifelse(time<max(time,na.rm = TRUE),0,d0))
               
  # model_formula = as.formula(paste("Surv(time0, time, d) ~ ",paste(SM_variables,collapse="+")))
  # model = coxph(model_formula, data = long_data_imp2) 
  # (alpha.hat = summary(model)$coef[,1])

  # fit the cox-ph model
  coxph_imp = function(data_imp){
    # cox-ph model
    long_data_imp2 = data_imp %>%
    group_by(id) %>%
      mutate(time0 = ifelse(is.na(lag(time)),0,lag(time) ))
    long_data_imp2$d0 = surv_data$d[match(long_data_imp2$id,surv_data$id)]
    long_data_imp2 = na.omit(long_data_imp2) %>%
                group_by(id) %>%
                mutate(d=ifelse(time<max(time,na.rm = TRUE),0,d0))
    model_formula = as.formula(paste("Surv(time0, time, d) ~ ",paste(SM_variables,collapse="+")))
    model = coxph(model_formula, data = long_data_imp2) 
    (alpha.hat = summary(model)$coef[,1])
  }
    fit = lapply(1:5, function(i) coxph_imp(complete(imp1, action=i)))
    alpha.hat = sapply(seq_along(fit[[1]]), function(i) mean(sapply(fit, `[`, i)))
    names(alpha.hat) = names(fit[[1]])

  return(alpha.hat)
}

YDimputation_fun = function(long_data, surv_data){
   # round up the time 
  long_data$time = ceiling(long_data$time)
  long_data = long_data[!duplicated(long_data[,c("id","time")]),]
  # create the full data
  id_all = unique(long_data$id)
  time.min = min(long_data$time,na.rm=TRUE)
  time.max = max(long_data$time,na.rm=TRUE)
  ni = time.max-time.min+1

  # create the full data
  base_data = long_data[,c("Z","X")]
  base_data = base_data[!duplicated(base_data),]
  base_Z = base_data[,1]
  base_X = base_data[,2]
  df_full = data.frame(id = rep(id_all,each=ni),
                       time = rep(1:ni,length(id_all)),
                       Z = rep(base_Z, each = ni),
                       X = rep(base_X,each = ni))
  # df_full2 = merge(df_full,long_data, by = c("id","time","Z","X"),all=TRUE)
  long_data_observed = na.omit(long_data)
  df_full2 = merge(df_full,long_data_observed, by = c("id","time","Z","X"),all=TRUE)

  # empty imputation
  imp0 = mice(as.matrix(df_full2),maxit=0) 
  predM = imp0$predictorMatrix 
  impM = imp0$method
  
  #specify predictor matrix and method
  predM1 = predM 
  predM1["Y","id"] = -2 
  predM1["Y",c("time","Z","X")] = 1 #fixed x effects imputation 
  impM1 = impM 
  impM1["Y"]="2l.lmer"
  
  # multilevel imputation 
  imp1=mice(as.matrix(df_full2), m=5, 
            predictorMatrix=predM1, method=impM1, maxit=5 ) 
  long_data_imp = complete(imp1) # extract the complete data

  # cox-ph model
  long_data_imp2 = long_data_imp %>%
    group_by(id) %>%
    mutate(time0 = ifelse(is.na(lag(time)),0,lag(time) ))
  long_data_imp2$d = surv_data$d[match(long_data_imp2$id,surv_data$id)]
  model = coxph(Surv(time0, time, d) ~ Z+X+Y, data = long_data_imp2) 
  (alpha.hat = summary(model)$coef[,1])

  #   # fit the cox-ph model
  # coxph_imp = function(data_imp){
  #   # cox-ph model
  #   long_data_imp2 = data_imp %>%
  #     group_by(id) %>%
  #     mutate(time0 = ifelse(is.na(lag(time)),0,lag(time) ))
  #   long_data_imp2$d0 = surv_data$d[match(long_data_imp2$id,surv_data$id)]
  #   long_data_imp2 = na.omit(long_data_imp2) %>%
  #               group_by(id) %>%
  #               mutate(d=ifelse(time<max(time,na.rm = TRUE),0,d0))
  #   model = coxph(Surv(time0, time, d) ~ Z+X+Y, data = long_data_imp2) 
  #   (alpha.hat = summary(model)$coef[,1])
  # }
  # fit = lapply(1:5, function(i) coxph_imp(complete(imp1, action=i)))
  # alpha.hat = sapply(seq_along(fit[[1]]), function(i) mean(sapply(fit, `[`, i)))
  # names(alpha.hat) = names(fit[[1]])
  
  return(alpha.hat)
}

summaryStat = function(long_data,surv_data){
  long_data1 = na.omit(long_data) %>%
              group_by(id) %>%
              summarise(Y_min=log(min(exp(Y))), Y_max=log(max(exp(Y))),Y_mean=log(mean(exp(Y))),Y_median=log(median(exp(Y))))
  surv_data1 = merge(surv_data,long_data1,by="id")

  model_formula = as.formula(paste("Surv(time=D, event=d) ~ ",paste(SM_variables[1:3],collapse="+"),"+Y_min"))
  time_start = Sys.time()
  model = coxph(model_formula,data=surv_data1)
  time_end = Sys.time()
  time_diff_min = as.numeric(difftime(time_end, time_start, units = "secs"))
  alpha.hat.min = summary(model)$coef[,1]

  model_formula = as.formula(paste("Surv(time=D, event=d) ~ ",paste(SM_variables[1:3],collapse="+"),"+Y_mean"))
  time.start = Sys.time()
  model = coxph(model_formula,data=surv_data1)
  time.end = Sys.time()
  time_diff_mean = as.numeric(difftime(time_end, time_start, units = "secs"))
  alpha.hat.mean = summary(model)$coef[,1]

  model_formula = as.formula(paste("Surv(time=D, event=d) ~ ",paste(SM_variables[1:3],collapse="+"),"+Y_median"))
  time.start = Sys.time()
  model = coxph(model_formula,data=surv_data1)
  time.end = Sys.time()
  time_diff_median = as.numeric(difftime(time_end, time_start, units = "secs"))
  alpha.hat.median = summary(model)$coef[,1]

  model_formula = as.formula(paste("Surv(time=D, event=d) ~ ",paste(SM_variables[1:3],collapse="+"),"+Y_max"))
  time.start = Sys.time()
  model = coxph(model_formula,data=surv_data1)
  time.end = Sys.time()
  time_diff_max = as.numeric(difftime(time_end, time_start, units = "secs"))
  alpha.hat.max = summary(model)$coef[,1]

  alpha_hat = rbind(alpha.hat.min,alpha.hat.mean,alpha.hat.median,alpha.hat.max)
  time_diff = c(time_diff_min,time_diff_mean,time_diff_median,time_diff_max)
  names(time_diff) = c("min","mean","median","max")

  results = list("alpha_hat" = alpha_hat,
                 "time_diff" = time_diff)
  
  return(results)
  
}

VAimputation_realData_year = function(long_data,surv_data){
    
  time_unit = year_unit=1
  month_unit=12/6
  data = long_data
  # monthly-grouped data
  data = data %>%
      mutate(time_6month = time*time_unit*2,time_roundup_6month = ceiling(time_6month))%>%
      group_by(id,time_roundup_6month) %>%
      mutate(Y_mean = mean(Y,na.rm = TRUE),Age_mean = mean(Age,na.rm=TRUE)) %>%
      dplyr::select(id,Y_mean,Age_mean,Sex,SNP,time_roundup_6month)
  data1 = data[!duplicated(data),]
  colnames(data1) = c("id","Y","Age","Sex","SNP","time")

  # create a new dataset with the same columns in data1. 
  # For a given id, the time is from 1 to max_time. 
  # If there is no value of Y at a given time, then Y is NA.
  max_time <- max(data1$time, na.rm = TRUE)
  id_all = rep(unique(data1$id),each = max_time)
  time_all = rep(as.numeric(1:max_time),length(unique(data1$id)))
  Age_all = rep(unique(data1[,c("id","Age")])$Age,each=max_time)
  Sex_all = rep(unique(data1[,c("id","Sex")])$Sex,each=max_time)
  SNP_all = rep(unique(data1[,c("id","SNP")])$SNP,each=max_time)
  data2 = data.frame(id = id_all,time = time_all,Sex = Sex_all, Age = Age_all,SNP = SNP_all)
  # delete time after the disease time
  data2$D = surv_data$D[match(data2$id,surv_data$id)]
  sum(data2$time<=ceiling(data2$D*time_unit*month_unit))

  ######
  data2 = data2[data2$time<=ceiling(data2$D*time_unit*month_unit),c("id","time","Age","Sex","SNP")]
  head(data2)

  # join data1 to data2 by id and time
  data3 = left_join(data2,data1,by = c("id","time","Age","Sex","SNP"))
  data3$time = data3$time/2/time_unit # convert time back to the original scale

  df_full2 = data3
  dim(df_full2)

  # insert the predictor Ni(t)
  df_full3 = df_full2 
  df_full3$Ni = NA
  for (i in 1:nrow(df_full3)){
    id = df_full3$id[i]
    time = df_full3$time[i]
    df_full3$Ni[i] = sum(!is.na(df_full3$Y[df_full3$id==id & df_full3$time<=time]))
  }
  df_full3 = df_full3[,c("id","time","Age","Sex","SNP","Ni","Y")]
  head(df_full3)

  # empty imputation
  imp0 = mice(as.matrix(df_full3),maxit=0) 
  predM = imp0$predictorMatrix 
  impM = imp0$method

  #specify predictor matrix and method
  predM1 = predM 
  predM1["Y","id"] = -2 
  predM1["Y",c(LM_fixedEffect_withTime_variables,"Ni")] = 1 #fixed x effects imputation 
  impM1 = impM 
  impM1["Y"]="2l.lmer"

  # multilevel imputation 
  imp1=mice(as.matrix(df_full3), m=5, 
          predictorMatrix=predM1, method=impM1, maxit=5 )

  # cox-ph model
  long_data_imp = complete(imp1)
  long_data_imp2 = long_data_imp %>%
    group_by(id) %>%
    mutate(time0 = ifelse(is.na(lag(time)),0,lag(time) ))
  long_data_imp2$d0 = surv_data$d[match(long_data_imp2$id,surv_data$id)]
  long_data_imp2 = na.omit(long_data_imp2) %>%
               group_by(id) %>%
               mutate(d=ifelse(time<max(time,na.rm = TRUE),0,d0))
               
  model_formula = as.formula(paste("Surv(time0, time, d) ~ ",paste(SM_variables,collapse="+")))
  model = coxph(model_formula, data = long_data_imp2) 
  (alpha.hat = summary(model)$coef[,1])

  # fit the cox-ph model
  coxph_imp = function(data_imp){
    # cox-ph model
    long_data_imp2 = data_imp %>%
    group_by(id) %>%
      mutate(time0 = ifelse(is.na(lag(time)),0,lag(time) ))
    long_data_imp2$d0 = surv_data$d[match(long_data_imp2$id,surv_data$id)]
    long_data_imp2 = na.omit(long_data_imp2) %>%
                group_by(id) %>%
                mutate(d=ifelse(time<max(time,na.rm = TRUE),0,d0))
    model_formula = as.formula(paste("Surv(time0, time, d) ~ ",paste(SM_variables,collapse="+")))
    model = coxph(model_formula, data = long_data_imp2) 
    (alpha.hat = summary(model)$coef[,1])
  }
    fit = lapply(1:5, function(i) coxph_imp(complete(imp1, action=i)))
    alpha.hat = sapply(seq_along(fit[[1]]), function(i) mean(sapply(fit, `[`, i)))
    names(alpha.hat) = names(fit[[1]])

  return(alpha.hat)
}

VAJM_fun_realData = function(long_data, surv_data){

  # remove patients with no Y measurement
  zeroRecords_id = long_data[is.na(long_data$Y),"id"]
  long_data = long_data[!long_data$id %in% zeroRecords_id,]
  surv_data = surv_data[!surv_data$id %in% zeroRecords_id,]  

  # add Ni(t) as a predictor
  long_data$Ni = NA
  for (i in 1:nrow(long_data)){
    id = long_data$id[i]
    time = long_data$time[i]
    long_data$Ni[i] = sum(!is.na(long_data$Y[long_data$id==id & long_data$time<=time]))
  }

  # longitudinal submodel, try nlminb optimizer first, if there is an error, then use optim optimizer
  lmeFit_fixedformula = as.formula(paste("Y ~ ",paste(LM_fixedEffect_withTime_variables,collapse="+"),"+Ni"))
  lmeFit_randomformula = as.formula(paste("~",paste(LM_randomEffect_variables,collapse="+"),"|id"))
  control_optim = lmeControl(msMaxIter = 1000, msMaxEval = 1000, opt='optim')
  control_nlminb = lmeControl(msMaxIter = 1000, msMaxEval = 1000, opt='nlminb')
  lmeFit = tryCatch({
    lmeFit = lme(lmeFit_fixedformula,random = lmeFit_randomformula, data = long_data, control=control_optim)
  }, error = function(e) {
    print(paste0("Error with optim: ",e))
    lmeFit = lme(lmeFit_fixedformula,random = lmeFit_randomformula, data = long_data, control=control_nlminb)
  })
  print(lmeFit)

  # survival submodel
  coxFit_formula = as.formula(paste("Surv(D, d) ~ ",paste(SM_base_variables,collapse="+")))
  coxFit = coxph(coxFit_formula, data = surv_data, x = TRUE)
  summary(coxFit)

  # jointFit = JMbayes2::jm(coxFit,lmeFit,time_var="time", n_chains=1L,n_iter=11000L,n_burnin=1000L) 
  jointFit = JMbayes2::jm(coxFit,lmeFit,time_var="time", n_chains=1L) 

  print(summary(jointFit))
  
  surv_proc = unlist(coef(jointFit))
  long_proc = unlist(fixef(jointFit))

  results = list("long_proc" = long_proc,
                 "surv_proc" = surv_proc)

  return(results)
}

