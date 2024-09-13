library(MASS)
library(dplyr)

simData_noIVP_randomZ = function(m,time.start,time.end,seed){
  
  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b = c(0,0)
  std_b = matrix(0,2,2)
  diag(std_b) = c(1,4) # control the variance of b
  b = as.matrix(mvrnorm(n=m,mu=mean_b,Sigma=std_b))
  beta = matrix(c("intercept"=-2, "Z"=-0.5, "X"=0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma = matrix(c("intercept"=-5, "Z"=1, "X"=0.5, "b"=0.5),ncol=1)
  
  ## survival process
  alpha = matrix(c("intercept" = -5, "Z"=-1, "X"=0.5, "Y"=1),ncol=1)
  
  set.seed(seed)
  
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit_data = cbind(data_base_covariates,b[,2])
  rates = exp(cbind(1,visit_data) %*% gamma)
#   visits_org = lapply(rates, visitingPattern,maxTime=time.end)
  visits_org = lapply(rates, function(x){seq(time.start,time.end,by=6)})
  
  visits_time = visits_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {
      Y_i = visits_time_i=NA
    } else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b[i,,drop=FALSE]) # Time has random effect
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y_i = eta_i+error
    }
    id = rep(i,length(Y_i))
    base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b[i,]),each=length(Y_i)),nrow=length(Y_i))
    df_i = cbind(id, base_i, visits_time_i, Y_i)
    df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b0","b1","time","Y")
  (df_full = as.data.frame(df_full))
  
  #### survival process ####
  generate_D = function(covars,J){
    meanY_fun = function(t){
      meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*covars$Z
      # meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*t
    }
    
    hazard_fun = function(t){
      meanY = meanY_fun(t)
      exp(alpha[1,]+alpha[2,]*covars$Z+alpha[3,]*covars$X+alpha[4,]*meanY)
    }
    
    surv_fun <- function(t){
      Lambda = integrate(f = hazard_fun, lower = 0, upper = t)$value
      exp(-Lambda)
    }
    
    u=runif(1,min=0,max=1)
    times = seq(0,J+2, by=0.1)
    times[1] = 1e-3
    surv_prob = unlist(lapply(times, surv_fun))
    D = times[nearest(x=surv_prob,u=u)]
    d = ifelse(D>J,0,1)
    D = min(D,(J))
    
    return(list("D"=D,"d"=d))
  }
  
  ###### create the disease outcome ######
  data.id = df_full[,c("id","Z","X","b0","b1")]
  data.id = data.id[!duplicated(data.id),]
  
  D = d = NULL
  for(i in 1:nrow(data.id)){
    covars = data.id[i,]
    temp = generate_D(covars,J)
    D = c(D,temp$D)
    d = c(d,temp$d)
  }
  D
  mean(d)
  mean(D)

  surv_data = data.id[,c("id","Z","X"),drop=FALSE]
  surv_data$D = D
  surv_data$d = d

   #### longitudinal data ###
  head(df_full)
  apply(is.na(df_full),2,mean)
  df_full = df_full[,c("id","Z","X","time","Y")]
  rep_times = visits_length
  rep_times[which(rep_times==0)] = 1
  df_full$maxTime = rep(D,rep_times)

  long_data_observed = na.omit(df_full)
  long_data_observed1 = long_data_observed[which(long_data_observed$time<=long_data_observed$maxTime),c("id","Z","X","time" ,"Y")]

  id_zeroObs = setdiff(1:m,unique(long_data_observed1$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = long_data_observed1
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(long_data_observed1,df_zeroObs)
  }
  long_data = df_full[,c("id","Z","X","time","Y")]  
  
  #### output ####
  zeroRecords_prop = mean(is.na(long_data[!duplicated(long_data$id),"Y"])) # proportion of patients with no Y measurement
  print(paste0("The proportion of patients with no longitudinal measurements is ", zeroRecords_prop))

  data = list(long_data = long_data,
              surv_data = surv_data)
  
  return(data)
  
}

simData_withIVPmeasured_randomZ = function(m,time.start,time.end,seed){
  
  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b = c(0,0)
  std_b = matrix(0,2,2)
  diag(std_b) = c(1,4) # control the variance of b
  b = as.matrix(mvrnorm(n=m,mu=mean_b,Sigma=std_b))
  beta = matrix(c("intercept"=-2, "Z"=-0.5, "X"=0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma = matrix(c("intercept"=-2.2, "Z"=0.5, "X"=0.5,  "b"=0),ncol=1)
  
  ## survival process
  alpha = matrix(c(c("intercept" = -5, "Z"=-1, "X"=0.5, "Y"=1)),ncol=1)
  
  set.seed(seed)
  
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit_data = cbind(data_base_covariates,b[,2])
  rates = exp(cbind(1,visit_data) %*% gamma)
  visits_org = lapply(rates, visitingPattern,maxTime=time.end)
  # visits_org = lapply(rates, function(x){time.start:time.end})
  
  visits_time = visits_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {
      Y_i = visits_time_i=NA
    } else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b[i,,drop=FALSE]) # Time has random effect
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y_i = eta_i+error
    }
    id = rep(i,length(Y_i))
    base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b[i,]),each=length(Y_i)),nrow=length(Y_i))
    df_i = cbind(id, base_i, visits_time_i, Y_i)
    df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b0","b1","time","Y")
  (df_full = as.data.frame(df_full))
  
  #### survival process ####
  generate_D = function(covars,J){
    meanY_fun = function(t){
      meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*covars$Z
      # meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*t
    }
    
    hazard_fun = function(t){
      meanY = meanY_fun(t)
      exp(alpha[1,]+alpha[2,]*covars$Z+alpha[3,]*covars$X+alpha[4,]*meanY)
    }
    
    surv_fun <- function(t){
      Lambda = integrate(f = hazard_fun, lower = 0, upper = t)$value
      exp(-Lambda)
    }
    
    u=runif(1,min=0,max=1)
    times = seq(0,J+2, by=0.1)
    times[1] = 1e-3
    surv_prob = unlist(lapply(times, surv_fun))
    D = times[nearest(x=surv_prob,u=u)]
    d = ifelse(D>J,0,1)
    D = min(D,(J))
    
    return(list("D"=D,"d"=d))
  }
  
  ###### create the disease outcome ######
  data.id = df_full[,c("id","Z","X","b0","b1")]
  data.id = data.id[!duplicated(data.id),]

  D = d = NULL
  for(i in 1:nrow(data.id)){
    covars = data.id[i,]
    temp = generate_D(covars,J)
    D = c(D,temp$D)
    d = c(d,temp$d)
  }
  D
  mean(d)
  mean(D)
  
  surv_data = data.id[,c("id","Z","X"),drop=FALSE]
  surv_data$D = D
  surv_data$d = d
  
  #### longitudinal data ###
  head(df_full)
  apply(is.na(df_full),2,mean)
  df_full = df_full[,c("id","Z","X","time","Y")]
  rep_times = visits_length
  rep_times[which(rep_times==0)] = 1
  df_full$maxTime = rep(D,rep_times)

  long_data_observed = na.omit(df_full)
  long_data_observed1 = long_data_observed[which(long_data_observed$time<=long_data_observed$maxTime),c("id","Z","X","time" ,"Y")]

  id_zeroObs = setdiff(1:m,unique(long_data_observed1$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = long_data_observed1
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(long_data_observed1,df_zeroObs)
  }
  long_data = df_full[,c("id","Z","X","time","Y")]  

#   temp = long_data %>%
#     group_by(id) %>%
#     summarise(n=n())
#   temp$n
    
  #### output ####
  zeroRecords_prop = mean(is.na(long_data[!duplicated(long_data$id),"Y"])) # proportion of patients with no Y measurement
  print(paste0("The proportion of patients with no longitudinal measurements is ", zeroRecords_prop))

  data = list(long_data = long_data,
              surv_data = surv_data)
  
  return(data)
  
}

simData_withIVPthres_randomZ = function(m,time.start,time.end,seed){
  
  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b = c(0,0)
  std_b = matrix(0,2,2)
  diag(std_b) = c(1,4) # control the variance of b
  b = as.matrix(mvrnorm(n=m,mu=mean_b,Sigma=std_b))
  beta = matrix(c("intercept"=-2, "Z"=-0.5, "X"=0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma = matrix(c("intercept"=-1, "Z"=1, "X"=0.5,  "b"=0.5),ncol=1)
  
  ## survival process
  alpha = matrix(c(c("intercept" = -5, "Z"=-1, "X"=0.5, "Y"=1)),ncol=1)
  
  set.seed(seed)
  
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit_data = cbind(data_base_covariates,b[,2])
  rates = exp(cbind(1,visit_data) %*% gamma)
  # visits_org = lapply(rates, visitingPattern,maxTime=time.end)
  visits_org = lapply(rates, function(x){time.start:time.end})
  
  visits_time = visits_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {
      Y_i = visits_time_i=NA
    } else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b[i,,drop=FALSE]) # Time has random effect
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y_i = eta_i+error
    }
    id = rep(i,length(Y_i))
    base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b[i,]),each=length(Y_i)),nrow=length(Y_i))
    df_i = cbind(id, base_i, visits_time_i, Y_i)
    df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b0","b1","time","Y")
  (df_full = as.data.frame(df_full))
  
  ### thresholding the data ###
  thres = quantile(df_full$Y,0.8)
  df_full$Y = ifelse(df_full$Y<thres,NA,df_full$Y)
  
  #### survival process ####
  generate_D = function(covars,J){
    meanY_fun = function(t){
      meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*covars$Z
      # meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*t
    }
    
    hazard_fun = function(t){
      meanY = meanY_fun(t)
      exp(alpha[1,]+alpha[2,]*covars$Z+alpha[3,]*covars$X+alpha[4,]*meanY)
    }
    
    surv_fun <- function(t){
      Lambda = integrate(f = hazard_fun, lower = 0, upper = t)$value
      exp(-Lambda)
    }
    
    u=runif(1,min=0,max=1)
    times = seq(0,J+2, by=0.1)
    times[1] = 1e-3
    surv_prob = unlist(lapply(times, surv_fun))
    D = times[nearest(x=surv_prob,u=u)]
    d = ifelse(D>J,0,1)
    D = min(D,(J))
    
    return(list("D"=D,"d"=d))
  }
  
  ###### create the disease outcome ######
  data.id = df_full[,c("id","Z","X","b0","b1")]
  data.id = data.id[!duplicated(data.id),]

  D = d = NULL
  for(i in 1:nrow(data.id)){
    covars = data.id[i,]
    temp = generate_D(covars,J)
    D = c(D,temp$D)
    d = c(d,temp$d)
  }
  D
  mean(d)
  mean(D)
  
  surv_data = data.id[,c("id","Z","X"),drop=FALSE]
  surv_data$D = D
  surv_data$d = d
  
 #### longitudinal data ###
  head(df_full)
  apply(is.na(df_full),2,mean)
  df_full = df_full[,c("id","Z","X","time","Y")]
  rep_times = visits_length
  rep_times[which(rep_times==0)] = 1
  df_full$maxTime = rep(D,rep_times)

  long_data_observed = na.omit(df_full)
  long_data_observed1 = long_data_observed[which(long_data_observed$time<=long_data_observed$maxTime),c("id","Z","X","time" ,"Y")]

  id_zeroObs = setdiff(1:m,unique(long_data_observed1$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = long_data_observed1
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(long_data_observed1,df_zeroObs)
  }
  long_data = df_full[,c("id","Z","X","time","Y")]  

#   temp = long_data %>%
#     group_by(id) %>%
#     summarise(n=n())
#   temp$n
  
  #### output ####
  zeroRecords_prop = mean(is.na(long_data[!duplicated(long_data$id),"Y"])) # proportion of patients with no Y measurement
  print(paste0("The proportion of patients with no longitudinal measurements is ", zeroRecords_prop))

  data = list(long_data = long_data,
              surv_data = surv_data)
  
  return(data)
  
}

simData_withIVPprevious_randomZ = function(m,time.start,time.end,seed){
  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b = c(0,0)
  std_b = matrix(0,2,2)
  diag(std_b) = c(1,4) # control the variance of b
  b = as.matrix(mvrnorm(n=m,mu=mean_b,Sigma=std_b))
  beta = matrix(c("intercept"=-2, "Z"=-0.5, "X"=0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma = matrix(c("intercept"=-2.2, "Z"=0, "X"=0,  "Y_prev"=1),ncol=1)
  
  ## survival process
  alpha = matrix(c("intercept" = -5, "Z"=-1, "X"=0.5, "Y"=1),ncol=1)
  
  set.seed(seed)
  
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit_data = cbind(data_base_covariates,b[,2])
  rates = exp(cbind(1,visit_data) %*% gamma)
  # visits_org = lapply(rates, visitingPattern,maxTime=time.end)
  visits_org = lapply(rates, function(x){0:time.end})
  
  visits_time = visits_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {
      Y_i = visits_time_i=NA
    } else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b[i,,drop=FALSE]) # Time has random effect
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y_i = eta_i+error
    }
    id = rep(i,length(Y_i))
    base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b[i,]),each=length(Y_i)),nrow=length(Y_i))
    df_i = cbind(id, base_i, visits_time_i, Y_i)
    df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b0","b1","time","Y")
  (df_full = as.data.frame(df_full))
  
  #### survival process ####
  generate_D = function(covars,J){
    meanY_fun = function(t){
      meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*covars$Z
      # meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*t
    }
    
    hazard_fun = function(t){
      meanY = meanY_fun(t)
      exp(alpha[1,]+alpha[2,]*covars$Z+alpha[3,]*covars$X+alpha[4,]*meanY)
    }
    
    surv_fun <- function(t){
      Lambda = integrate(f = hazard_fun, lower = 0, upper = t)$value
      exp(-Lambda)
    }
    
    u=runif(1,min=0,max=1)
    times = seq(0,J+2, by=0.1)
    times[1] = 1e-3
    surv_prob = unlist(lapply(times, surv_fun))
    D = times[nearest(x=surv_prob,u=u)]
    d = ifelse(D>J,0,1)
    D = min(D,(J))
    
    return(list("D"=D,"d"=d))
  }
  
  ###### create the disease outcome ######
  data.id = df_full[,c("id","Z","X","b0","b1")]
  data.id = data.id[!duplicated(data.id),]
  
  D = d = NULL
  for(i in 1:nrow(data.id)){
    covars = data.id[i,]
    temp = generate_D(covars,J)
    D = c(D,temp$D)
    d = c(d,temp$d)
  }
  D
  mean(d)
  mean(D)
  
  surv_data = data.id[,c("id","Z","X"),drop=FALSE]
  surv_data$D = D
  surv_data$d = d
  
  #### visiting pattern ####
  visit_data = as.matrix(df_full[,c("Z","X","Y")])
  rates = exp(cbind(1,visit_data) %*% gamma)
  
  visits_union = 0:time.end
  visits =  NULL
  for(i in 1:m){
    rates_i = rates[which(df_full$id==i)]
    visits_temp = visitingPattern3(rate=rates_i,maxTime = time.end)
    visits_temp2 = ifelse(visits_union %in% visits_temp, 1, NA)
    visits = c(visits, visits_temp2)
  }
  df_full$R = visits
  df_full$Y_obs = df_full$Y*df_full$R
  # remove time=0
  df_full = df_full[which(df_full$time!=0),]
  
 #### longitudinal data ###
  head(df_full)
  apply(is.na(df_full),2,mean)
  df_full = df_full[,c("id","Z","X","time","Y_obs")]
  colnames(df_full) = c("id","Z","X","time","Y")
  rep_times = visits_length-1
  rep_times[which(rep_times==0)] = 1
  df_full$maxTime = rep(D,rep_times)

  long_data_observed = na.omit(df_full)
  long_data_observed1 = long_data_observed[which(long_data_observed$time<=long_data_observed$maxTime),c("id","Z","X","time" ,"Y")]

  id_zeroObs = setdiff(1:m,unique(long_data_observed1$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = long_data_observed1
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(long_data_observed1,df_zeroObs)
  }
  long_data = df_full[,c("id","Z","X","time","Y")]  

#   temp = long_data %>%
#     group_by(id) %>%
#     summarise(n=n())
#   temp$n
  
  #### output ####
  zeroRecords_prop = mean(is.na(long_data[!duplicated(long_data$id),"Y"])) # proportion of patients with no Y measurement
  print(paste0("The proportion of patients with no longitudinal measurements is ", zeroRecords_prop))

  data = list(long_data = long_data,
              surv_data = surv_data)
  
  return(data)
  
}

simData_withIVPunmeasuredGamma_randomZ = function(m,time.start,time.end,seed){
  
  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b = c(0,0)
  std_b = matrix(0,2,2)
  diag(std_b) = c(1,4) # control the variance of b
  b = as.matrix(mvrnorm(n=m,mu=mean_b,Sigma=std_b))
  beta = matrix(c("intercept"=-2, "Z"=-0.5, "X"=0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma = matrix(c("intercept"=-3.5, "Z"=1, "X"=1),ncol=1)
  
  ## survival process
  alpha = matrix(c(c("intercept" = -5, "Z"=-1, "X"=0.5, "Y"=1)),ncol=1)
  
  set.seed(seed)
  
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit_data = cbind(Z,X)
  rates_base = exp(cbind(1,visit_data) %*% gamma)
  theta0 = matrix(c("b0"=0,"b1"=0.2),ncol=1)
  # eta = rgamma(m,shape=1/exp(b1%*%theta0),rate = 1/exp(b1%*%theta0))
  eta = rgamma(m,shape=1/exp(b%*%theta0),rate = 1)
  mean(eta);var(eta)
  rates= eta*rates_base

  visits_org = lapply(rates, visitingPattern,maxTime=time.end)
  # visits_org = lapply(rates, function(x){time.start:time.end})
  
  visits_time = visits_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {
      Y_i = visits_time_i=NA
    } else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b[i,,drop=FALSE]) # Time has random effect
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y_i = eta_i+error
    }
    id = rep(i,length(Y_i))
    base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b[i,]),each=length(Y_i)),nrow=length(Y_i))
    df_i = cbind(id, base_i, visits_time_i, Y_i)
    df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b0","b1","time","Y")
  (df_full = as.data.frame(df_full))
  
  #### survival process ####
  generate_D = function(covars,J){
    meanY_fun = function(t){
      meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*covars$Z
      # meanY = beta[1,]+covars$Z*beta[2,]+covars$X*beta[3,]+t*beta_t+covars$b0+covars$b1*t
    }
    
    hazard_fun = function(t){
      meanY = meanY_fun(t)
      exp(alpha[1,]+alpha[2,]*covars$Z+alpha[3,]*covars$X+alpha[4,]*meanY)
    }
    
    surv_fun <- function(t){
      Lambda = integrate(f = hazard_fun, lower = 0, upper = t)$value
      exp(-Lambda)
    }
    
    u=runif(1,min=0,max=1)
    times = seq(0,J+2, by=0.1)
    times[1] = 1e-3
    surv_prob = unlist(lapply(times, surv_fun))
    D = times[nearest(x=surv_prob,u=u)]
    d = ifelse(D>J,0,1)
    D = min(D,(J))
    
    return(list("D"=D,"d"=d))
  }
  
  ###### create the disease outcome ######
  data.id = df_full[,c("id","Z","X","b0","b1")]
  data.id = data.id[!duplicated(data.id),]

  D = d = NULL
  for(i in 1:nrow(data.id)){
    covars = data.id[i,]
    temp = generate_D(covars,J)
    D = c(D,temp$D)
    d = c(d,temp$d)
  }
  D
  mean(d)
  mean(D)
  
  surv_data = data.id[,c("id","Z","X"),drop=FALSE]
  surv_data$D = D
  surv_data$d = d
  
  
 #### longitudinal data ###
  head(df_full)
  apply(is.na(df_full),2,mean)
  df_full = df_full[,c("id","Z","X","time","Y")]
  rep_times = visits_length
  rep_times[which(rep_times==0)] = 1
  df_full$maxTime = rep(D,rep_times)

  long_data_observed = na.omit(df_full)
  long_data_observed1 = long_data_observed[which(long_data_observed$time<=long_data_observed$maxTime),c("id","Z","X","time" ,"Y")]

  id_zeroObs = setdiff(1:m,unique(long_data_observed1$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = long_data_observed1
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(long_data_observed1,df_zeroObs)
  }
  long_data = df_full[,c("id","Z","X","time","Y")]  

#   temp = long_data %>%
#     group_by(id) %>%
#     summarise(n=n())
#   temp$n
  
  #### output ####
  zeroRecords_prop = mean(is.na(long_data[!duplicated(long_data$id),"Y"])) # proportion of patients with no Y measurement
  print(paste0("The proportion of patients with no longitudinal measurements is ", zeroRecords_prop))

  data = list(long_data = long_data,
              surv_data = surv_data)
  
  return(data)
  
}

simData_mixture_randomZ = function(m,time.start,time.end,seed){
    set.seed(seed)

    id_class = apply(rmultinom(m, size=1, prob=rep(1/5,5)),2,function(x){which(x==1)}) 
    n1 = sum(id_class==1)
    n2 = sum(id_class==2)
    n3 = sum(id_class==3)
    n4 = sum(id_class==4)
    n5 = sum(id_class==5)

    model1 = simData_noIVP_randomZ(n1,time.start,time.end,seed)
    long_data1 = model1$long_data
    surv_data1 = model1$surv_data

    model2 = simData_withIVPmeasured_randomZ(n2,time.start,time.end,seed)
    long_data2 = model2$long_data
    surv_data2 = model2$surv_data
    long_data2$id = long_data2$id + n1
    surv_data2$id = surv_data2$id + n1

    model3 = simData_withIVPunmeasuredGamma_randomZ(n3,time.start,time.end,seed)
    long_data3 = model3$long_data
    surv_data3 = model3$surv_data
    long_data3$id = long_data3$id + n1 + n2
    surv_data3$id = surv_data3$id + n1 + n2

    model4 = simData_withIVPprevious_randomZ(n4,time.start,time.end,seed)
    long_data4 = model4$long_data
    surv_data4 = model4$surv_data
    long_data4$id = long_data4$id + n1 + n2 + n3
    surv_data4$id = surv_data4$id + n1 + n2 + n3

    model5 = simData_withIVPthres_randomZ(n5,time.start,time.end,seed)
    long_data5 = model5$long_data
    surv_data5 = model5$surv_data
    long_data5$id = long_data5$id + n1 + n2 + n3 + n4
    surv_data5$id = surv_data5$id + n1 + n2 + n3 + n4

    long_data = rbind(long_data1,long_data2,long_data3,long_data4,long_data5)
    surv_data = rbind(surv_data1,surv_data2,surv_data3,surv_data4,surv_data5)
    data = list(long_data = long_data,
                surv_data = surv_data)
    
    return(data)

}
