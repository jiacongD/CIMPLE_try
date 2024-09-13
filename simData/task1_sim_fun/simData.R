library(MASS)
library(dplyr)

# m=100
# time.start=1
# time.end=60
# seed=1

simData_noIVP_randomZ = function(m,time.start,time.end,seed){
  
  set.seed(seed)

  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b1 = c(0,0)
  std_b1 = matrix(0,2,2)
  diag(std_b1) = c(1,4) # control the variance of b
  b1 = as.matrix(mvrnorm(n=m,mu=mean_b1,Sigma=std_b1))
  beta = matrix(c("intercept"=-2, "Z"=-0.5, "X"=0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma1 = matrix(c("intercept"=-1, "Z"=1, "X"=0.5, "b"=0.5),ncol=1)
    
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit1_data = cbind(data_base_covariates,b1[,2])
  rates1 = exp(cbind(1,visit1_data) %*% gamma1)
  visits1_org = lapply(rates1, function(x){seq(time.start,time.end,by=6)})
  
  visits_time = visits1_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {
      next 
    } else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b1[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b1[i,,drop=FALSE]) # Time has random effect
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y1_i = eta_i+error
    }
      id = rep(i,length(Y1_i))
      base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b1[i,]),each=length(Y1_i)),nrow=length(Y1_i))
      df_i = cbind(id, base_i, visits_time_i, Y1_i)
      df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b_0","b_1","time","Y")
  (df_full = as.data.frame(df_full))
  
  ###
  df_observed = df_full[,c("id","Z","X","time","Y")]
  id_zeroObs = setdiff(1:m,unique(df_observed$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = df_observed
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(df_observed,df_zeroObs)
  }

  #### longitudinal data ###
  long_data = df_full[,c("id","Z","X","time","Y")]
  
  id_crosswalk = data.frame("id" = unique(long_data$id),"newid"= 1:length(unique(long_data$id)))
  long_data$id = id_crosswalk$newid[match(long_data$id,id_crosswalk$id)]
  
  data = list(long_data = long_data)
  
  return(data)
  
}

simData_withIVPmeasured_randomZ = function(m,time.start,time.end,seed){
  
  set.seed(seed)

  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b1 = c(0,0)
  std_b1 = matrix(0,2,2)
  diag(std_b1) = c(1,4) # control the variance of b
  b1 = as.matrix(mvrnorm(n=m,mu=mean_b1,Sigma=std_b1))
  beta = matrix(c(-2,-0.5,0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma1 = matrix(c(-2.2,  0.5, 0.5,   0),ncol=1)
    
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit1_data = cbind(data_base_covariates,b1[,2])
  rates1 = exp(cbind(1,visit1_data) %*% gamma1)
  visits1_org = lapply(rates1, visitingPattern,maxTime=time.end)
  # visits1_org = lapply(rates1, function(x){time.start:time.end})
  
  visits_time = visits1_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {next} else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b1[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b1[i,,drop=FALSE]) # Time has random effect      
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y1_i = eta_i+error
    }
      id = rep(i,length(Y1_i))
      base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b1[i,]),each=length(Y1_i)),nrow=length(Y1_i))
      df_i = cbind(id, base_i, visits_time_i, Y1_i)
      df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b_0","b_1","time","Y")
  (df_full = as.data.frame(df_full))

  #### 
  df_observed = df_full[,c("id","Z","X","time","Y")]
  id_zeroObs = setdiff(1:m,unique(df_observed$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = df_observed
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(df_observed,df_zeroObs)
  }
  
  #### longitudinal data ###
  long_data = df_full[,c("id","Z","X","time","Y")]
  
  id_crosswalk = data.frame("id" = unique(long_data$id),"newid"= 1:length(unique(long_data$id)))
  long_data$id = id_crosswalk$newid[match(long_data$id,id_crosswalk$id)]
  
  data = list(long_data = long_data)
  
  return(data)
  
}

simData_withIVPprevious_randomZ = function(m,time.start,time.end,seed){

  set.seed(seed)

  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b1 = c(0,0)
  std_b1 = matrix(0,2,2)
  diag(std_b1) = c(1,4) # control the variance of b
  b1 = as.matrix(mvrnorm(n=m,mu=mean_b1,Sigma=std_b1))
  beta = matrix(c("intercept"=-2, "Z"=-0.5, "X"=0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma1 = matrix(c("intercept"=-2.2, "Z"=0, "X"=0, "Y_prev"=1),ncol=1)
    
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit1_data = cbind(data_base_covariates,b1[,2])
  rates1 = exp(cbind(1,visit1_data) %*% gamma1)
  # visits1_org = lapply(rates1, visitingPattern,maxTime=time.end)
  visits1_org = lapply(rates1, function(x){0:time.end})
  
  visits_time = visits1_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {next} else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b1[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b1[i,,drop=FALSE]) # Time has random effect
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y1_i = eta_i+error
    }
      id = rep(i,length(Y1_i))
      base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b1[i,]),each=length(Y1_i)),nrow=length(Y1_i))
      df_i = cbind(id, base_i, visits_time_i, Y1_i)
      df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b_0","b_1","time","Y")
  (df_full = as.data.frame(df_full))
  
  #### visiting pattern ####
  visit_data = as.matrix(df_full[,c("Z","X","Y")])
  rates = exp(cbind(1,visit_data) %*% gamma1)
  
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
  
  ###
  df_observed = df_full[,c("id","Z","X","time","Y_obs")]
  colnames(df_observed) = c("id","Z","X","time","Y")
  id_zeroObs = setdiff(1:m,unique(df_observed$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = df_observed
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(df_observed,df_zeroObs)
  }
  
  #### longitudinal data ###
  long_data = na.omit(df_full[,c("id","Z","X","time","Y")])
  colnames(long_data) = c("id","Z","X","time","Y")
  
  id_crosswalk = data.frame("id" = unique(long_data$id),"newid"= 1:length(unique(long_data$id)))
  long_data$id = id_crosswalk$newid[match(long_data$id,id_crosswalk$id)]
  
  data = list(long_data = long_data)
  
  return(data)
  
}

simData_withIVPthres_randomZ = function(m,time.start,time.end,seed){
  
  set.seed(seed)

  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b1 = c(0,0)
  std_b1 = matrix(0,2,2)
  diag(std_b1) = c(1,4) # control the variance of b
  b1 = as.matrix(mvrnorm(n=m,mu=mean_b1,Sigma=std_b1))
  beta = matrix(c(-2,-0.5,0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma1 = matrix(c(-1,  1, 0.5,   2),ncol=1)
    
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit1_data = cbind(data_base_covariates,b1[,2])
  rates1 = exp(cbind(1,visit1_data) %*% gamma1)
  # visits1_org = lapply(rates1, visitingPattern,maxTime=time.end)
  visits1_org = lapply(rates1, function(x){time.start:time.end})
  
  visits_time = visits1_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {next} else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b1[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b1[i,,drop=FALSE]) # Time has random effect
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y1_i = eta_i+error
    }
      id = rep(i,length(Y1_i))
      base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b1[i,]),each=length(Y1_i)),nrow=length(Y1_i))
      df_i = cbind(id, base_i, visits_time_i, Y1_i)
      df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b_0","b_1","time","Y")
  (df_full = as.data.frame(df_full))
  
  ### thresholding the data ###
  thres = quantile(df_full$Y,0.8)
  df_full = df_full[which(df_full$Y>=thres),]
  
  ####
  df_observed = df_full[,c("id","Z","X","time","Y")]
  id_zeroObs = setdiff(1:m,unique(df_observed$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = df_observed
  } else{
    print("have zero-observation individuals")
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(df_observed,df_zeroObs)
  }
  dim(df_full)

  #### longitudinal data ###
  long_data = df_full[,c("id","Z","X","time","Y")]
  
  id_crosswalk = data.frame("id" = unique(long_data$id),"newid"= 1:length(unique(long_data$id)))
  long_data$id = id_crosswalk$newid[match(long_data$id,id_crosswalk$id)]
  
  temp = long_data %>%
    group_by(id) %>%
    summarise(n=n())
  temp$n
  
  data = list(long_data = long_data)
  
  return(data)
  
}

# eta=exp(b_i*gamma)
simData_withIVPunmeasuredGamma_randomZ = function(m,time.start,time.end,seed){
  
  set.seed(seed)

  #### set parameters
  # longitudinal model
  ## marker 1
  mean_b1 = c(0,0)
  std_b1 = matrix(0,2,2)
  diag(std_b1) = c(1,4) # control the variance of b
  b1 = as.matrix(mvrnorm(n=m,mu=mean_b1,Sigma=std_b1))
  beta = matrix(c(-2,-0.5,0.5),ncol=1)
  beta_t = 0.1
  
  ## visiting process
  gamma1 = matrix(c("intercept" = -3.5,  "Z"=1, "X"=1),ncol=1)
    
  id=1:m
  J = time.end
  Z=rbinom(m,1,prob=0.5)
  # Z=rnorm(m,1,1)
  X=rnorm(m,0,1)
  data_base_covariates = cbind(Z,X)
  
  # visiting process
  visit1_data = cbind(data_base_covariates)
  rates1_base = exp(cbind(1,visit1_data) %*% gamma1)
  theta0 = matrix(c("b0"=0,"b1"=0.2),ncol=1)
  # eta = rgamma(m,shape=1/exp(b1%*%theta0),rate = 1/exp(b1%*%theta0))
  eta = rgamma(m,shape=1/exp(b1%*%theta0),rate = 1)
  mean(eta);var(eta)
  rates1= eta*rates1_base
  visits1_org = lapply(rates1, visitingPattern,maxTime=time.end)
  # visits1_org = lapply(rates1, function(x){time.start:time.end})
  
  visits_time = visits1_org
  visits_length = unlist(lapply(visits_time,function(x){length(x)}))
  
  df = NULL
  for(i in 1:m){
    if(is.null(visits_time[[i]])) {next} else{
      visits_time_i = matrix(visits_time[[i]],ncol=1)
      time_inv_effect = beta[1,]+data_base_covariates[i,c("Z","X")] %*% beta[-1,]
      time_effect = beta_t*ft(visits_time_i)
      random_effect = cbind(1,data_base_covariates[i,1]) %*% t(b1[i,,drop=FALSE]) # Z has random effect
      # random_effect = cbind(1,visits_time_i) %*% t(b1[i,,drop=FALSE]) # Time has random effect      
      error = matrix(rnorm(n=length(visits_time_i),0,1),ncol=1)
      eta_i = as.numeric(time_inv_effect)+time_effect+as.numeric(random_effect)
      Y1_i = eta_i+error
    }
      id = rep(i,length(Y1_i))
      base_i = matrix(rep(c(as.numeric(data_base_covariates[i,]),b1[i,]),each=length(Y1_i)),nrow=length(Y1_i))
      df_i = cbind(id, base_i, visits_time_i, Y1_i)
      df = rbind(df,df_i)
  }
  df_full = df
  colnames(df_full) = c("id","Z","X","b_0","b_1","time","Y")
  (df_full = as.data.frame(df_full))

  ####
  df_observed = df_full[,c("id","Z","X","time","Y")]
  id_zeroObs = setdiff(1:m,unique(df_observed$id))
  if(length(id_zeroObs)==0) {
    print("no zero-observation individuals")
    df_full = df_observed
  } else{
    print(paste0("the proportion of zero-observation individuals is ", length(id_zeroObs)/m))
    df_zeroObs = data.frame("id"=id_zeroObs,
                          "Z" = Z[match(id_zeroObs,1:m)],
                          "X" = X[match(id_zeroObs,1:m)],
                          "time"=NA,
                          "Y"=NA)
    df_full = rbind(df_observed,df_zeroObs)
  }

  dim(df_full)
  
  #### longitudinal data ###
  long_data = df_full[,c("id","Z","X","time","Y")]
  
  id_crosswalk = data.frame("id" = unique(long_data$id),"newid"= 1:length(unique(long_data$id)))
  long_data$id = id_crosswalk$newid[match(long_data$id,id_crosswalk$id)]
  
  data = list(long_data = long_data)
  
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

    model1 = simData_noIVP_randomZ(n1,time.start,time.end,seed)$long_data

    model2 = simData_withIVPmeasured_randomZ(n2,time.start,time.end,seed)$long_data
    # re-number id in model2, add n1 to all id
    model2$id = model2$id + n1

    model3 = simData_withIVPunmeasuredGamma_randomZ(n3,time.start,time.end,seed)$long_data
    # re-number id in model3, add n1+n2 to all id
    model3$id = model3$id + n1 + n2

    model4 = simData_withIVPprevious_randomZ(n4,time.start,time.end,seed)$long_data
    # re-number id in model4, add n1+n2+n3 to all id
    model4$id = model4$id + n1 + n2 + n3

    model5 = simData_withIVPthres_randomZ(n5,time.start,time.end,seed)$long_data
    # re-number id in model5, add n1+n2+n3+n4 to all id
    model5$id = model5$id + n1 + n2 + n3 + n4

    long_data = rbind(model1,model2,model3,model4,model5)
    data = list(long_data = long_data)
    
    return(data)

}
