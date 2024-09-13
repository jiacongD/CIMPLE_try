library(lme4)
library(nleqslv)
library(mice)
library(pracma)
library(merlin)

# One biomarker

standard_LME_realData = function(long_data){
  # long_data$Sex = as.factor(long_data$Sex)
  lme_model_formula = paste("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+"),"+(1+",paste(LM_randomEffect_variables,collapse = "+"),"|id)")
  control_params <- lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb',kkt=FALSE,tol=0.2,maxit=20000))
  lme_model = lmer(lme_model_formula,data=long_data,REML = FALSE,
                   control=lmerControl(optCtrl=control_params))
  summary(lme_model)
  (beta_hat = summary(lme_model)$coef[,1])
  # lme_model_fixed_formula = as.formula(paste("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+")))
  # lme_model_random_formula = as.formula(paste("~1+",paste(LM_randomEffect_variables,collapse = "+"),"|id"))
  # lme_model = lme(lme_model_fixed_formula,data=na.omit(long_data),random=lme_model_random_formula, method="REML",control = lmeControl(opt = "optim",msMaxIter = 5000, msMaxEval = 5000))
  # (beta_hat = coef(summary(lme_model))[,1])

  result = list(beta_hat = beta_hat, 
                beta_sd = summary(lme_model)$coef[,2])
  
  return(result)
}

conditional_LME_realData = function(long_data){
  long_data = long_data %>%
    group_by(id) %>%
    mutate(n_visits = 1:length(id))

  long_data$Sex = as.factor(long_data$Sex)
  lme_model_formula = paste("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+"),"+n_visits","+(1+",paste(LM_randomEffect_variables,collapse = "+"),"|id)")
  lme_model = lmer(lme_model_formula,data=long_data,REML = FALSE,
                   control=lmerControl(optCtrl=list(optimizer ='optimx', optCtrl=list(method='L-BFGS-B'),maxfun=50000)))
  (beta_hat = summary(lme_model)$coef[,1])

  # lme_model_fixed_formula = as.formula(paste("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+"),"+n_visits"))
  # lme_model_random_formula = as.formula(paste("~1+",paste(LM_randomEffect_variables,collapse = "+"),"|id"))
  # lme_model = lme(lme_model_fixed_formula,data=na.omit(long_data),random=lme_model_random_formula, method="REML",control = lmeControl(opt = "optim",msMaxIter = 1000, msMaxEval = 1000))
  # beta_hat = coef(summary(lme_model))[,1]

  return(beta_hat)
}

joint_Liang_realData = function(long_data){
  # log-normal model
  
  data = long_data
  # data$time = round(data$time,1)
  
  # variables in the visiting process
  V= data[,c("id",VPM_variables)] # row = n
  VV = as.matrix(V[!duplicated(V),-1]) # row=m
  V = as.matrix(V[,VPM_variables]) # variables affecting the visiting process
  # variables in the longitudinal process
  X = data[,c("id",LM_fixedEffect_withoutTime_variables)]
  XX = as.matrix(X[!duplicated(X),-1])
  X = as.matrix(X[,LM_fixedEffect_withoutTime_variables])
  # variables with random effect in the longitudinal process
  Q=data[,LM_randomEffect_variables,drop=FALSE] # variables with random effect

  temp = data %>%
    group_by(id) %>%
    summarize(C = max(time),ni=n())
  C = rep(max(na.omit(data$time)),nrow(VV)) # censoring time
  ni=temp$ni # number of observations per subject
  id.unique = unique(data$id)
  id = data$id
  Time = data$time
  Time.unique = sort(unique(Time),decreasing = FALSE)
  y = data$Y
  # J = time.end
  
  gamma.est.fun = function(V,VV){
    # data,data_base,vp_pred
    m = nrow(VV)
    VV = as.matrix(VV)
    
    gamma_function = function(gamma){
      gamma = matrix(gamma,ncol=1)
      # compute the average value
      exp.gamma = exp(as.matrix(VV) %*% gamma)
      
      # V.denom = sapply(Time, function(t){sum((t<=C)*exp.gamma)})
      # V.numer = apply(VV,2,function(x){
      #   sapply(Time,function(t){sum( (t<=C)*exp.gamma*x )})
      # })
      V.denom = rep(sum(na.omit(exp.gamma)),length(Time))
      V.numer = matrix(rep(apply(na.omit(as.vector(exp.gamma)) * VV,2,sum),length(Time)),nrow=length(Time),byrow = TRUE)

      V.weightedAve = V.numer/V.denom
      (obj = apply(V-V.weightedAve,2,mean))
    }
    
    gamma = rep(0.1,ncol(V))
    gamma.solve = nleqslv(x=gamma,fn=gamma_function,control=list(trace=1))
    (gamma.hat = gamma.solve$x)
    
    return(gamma.hat)
  }
  
  (gamma.hat = gamma.est.fun(V,VV))
  names(gamma.hat) = colnames(V)
  print(gamma.hat)
  
  # Lambda_estimation
  exp_gamma = as.vector(exp(VV %*%gamma.hat))
  denom_Lambda = rep(sum(exp_gamma),length(Time.unique))
  Time_freqtable = table(Time)
  Time_freqtable1 = Time_freqtable[match(as.character(Time.unique),names(Time_freqtable))]
  Lambda_t_est = as.numeric(Time_freqtable1/denom_Lambda)
  Lambda_C_est= rep(sum(Lambda_t_est[Time.unique<=C[1]]),length(C))
  
  # sigma^2 estimation (log-normal)
  exp_gamma_Lambda = exp_gamma*Lambda_C_est
  sigma_hat_sq = max( log(sum(ni^2-ni)/sum(exp_gamma_Lambda^2)),0)
  # sigma_hat_sq = 1
  print(sigma_hat_sq)

  # estimate the conditional expectation
  tmp = log(Lambda_C_est)+VV%*%gamma.hat
  E_eta = rep(NA,nrow(VV))
  # Numerical integration approach 1
  for(i in 1:nrow(VV)){
    # (ni_i = ni[i])
    (ni_i = ni[i])
    (tmp_i = tmp[i])
    E_eta_numer_fun = function(y){
       (log(y)-tmp_i)* exp(-1/(2*sigma_hat_sq)*(log(y)-tmp_i)^2 + (ni_i-2+1/2)*log(y)-y)
    }
    E_eta_denom_fun = function(x){
        exp(-1/(2*sigma_hat_sq)*(log(x)-tmp_i)^2 + (ni_i-2+1/2)*log(x) - x)
    }

    # try integrate function E_eta_denom_fun using tryCatch function, if the output has error message, then set output to be NA 
    tryCatch({
      E_eta_denom_i = integrate(E_eta_denom_fun, lower=0, upper=200)$value
      E_eta_numer_i = integrate(E_eta_numer_fun, lower=0, upper=200)$value
      E_eta[i] = E_eta_numer_i / E_eta_denom_i
    }, error = function(e) {
      print(i)
    })
  }
  #  E_eta[which(is.na(E_eta))] = 1 # if the integral is infinite, set the expectation to be 1.
  
  # # Numerical integration approach 2
  # cc = gaussLaguerre(15)
  # for(i in 1:nrow(VV)){
  #   (ni_i = ni[i])
  #   (tmp_i = tmp[i])
  #   E_eta_numer_ij = cc$w*(cc$x^(ni_i-2+1/2)*(log(cc$x)-tmp_i)*exp(-1/(2*sigma_hat_sq)*(log(cc$x)-tmp_i)^2))
  #   E_eta_denom_ij = cc$w*(cc$x^(ni_i-2+1/2)*exp(-1/(2*sigma_hat_sq)*(log(cc$x)-tmp_i)^2))
  #   # sum E_eta_numer_ij if it is finite
  #   E_eta_numer_i = sum(E_eta_numer_ij[is.finite(E_eta_numer_ij)]) 
  #   E_eta_denom_i = sum(E_eta_denom_ij[is.finite(E_eta_denom_ij)]) 
  #   E_eta[i] = E_eta_numer_i / E_eta_denom_i
  # }

  # B estimation
  Q0 = cbind(1,Q)
  # Q0 = as.matrix(Q,ncol=1)
  Bhat_i = as.vector(E_eta+sigma_hat_sq/2)
  # Bhat_long = sapply(id,function(i){Bhat_i[id.unique==i]})
  id_freqtable = table(id)
  id_freqtable1 = id_freqtable[match(as.character(id.unique),names(id_freqtable))]
  Bhat_long = rep(Bhat_i, id_freqtable1)
  Bhat = Bhat_long*Q0
  Bhat_base0 = cbind(data$id,Bhat)
  Bhat_base = Bhat_base0[!duplicated(Bhat_base0),-1,drop=FALSE]

  # estimating longitudinal parameters
  denom = rep(sum(ni/Lambda_C_est),length(Time))
  # X_bar, XX has only time-invariant variables
  numer_X = apply(XX,2,function(x){
    rep(sum(ni*x/Lambda_C_est),length(Time))
  })
  numer_B = apply(Bhat_base,2,function(x){
    rep(sum(na.omit(ni*x/Lambda_C_est)),length(Time))
  })
  
  
  if (sigma_hat_sq != 0 ){
    X_bar = numer_X/denom
    B_bar = numer_B/denom
    design_M = as.matrix(cbind(V,Bhat))
    design_Mbar = as.matrix(cbind(X_bar, B_bar))
    Liang.function = function(beta){
      tmp1 = (design_M - design_Mbar)*as.vector(data$Y-design_M%*%beta)
      value = apply(na.omit(tmp1),2,sum)
      return(value)
    }
    beta = as.matrix(rep(0.1,ncol(design_M)),ncol=1)
    beta.solve = nleqslv(x=beta,fn=Liang.function,control=list(trace=1))
    (beta.hat = beta.solve$x[1:ncol(V)])
  }
  if (sigma_hat_sq == 0 ){
    X_bar = numer_X/denom
    B_bar = numer_B/denom
    design_M = as.matrix(cbind(V))
    design_Mbar = as.matrix(cbind(numer_X/denom))
    Liang.function = function(beta){
      tmp1 = (design_M - design_Mbar)*as.vector(data$Y-design_M%*%beta)
      value = apply(na.omit(tmp1),2,sum)
      return(value)
    }
    beta = as.matrix(rep(0,ncol(design_M)),ncol=1)
    beta.solve = nleqslv(x=beta,fn=Liang.function,control=list(trace=1))
    (beta.hat = beta.solve$x[1:ncol(V)])
  }
  names(beta.hat) = colnames(X)
  (beta.hat)

  results = list(gamma.hat = gamma.hat,
                 beta.hat = beta.hat)
  
  return(results)
}

joint_LiangGamma_realData = function(long_data){
  # log-normal model
  
  data = long_data
  # data$time = round(data$time,1)
  
  # variables in the visiting process
  V= data[,c("id",VPM_variables)] # row = n
  VV = as.matrix(V[!duplicated(V),-1]) # row=m
  V = as.matrix(V[,VPM_variables]) # variables affecting the visiting process
  # variables in the longitudinal process
  X = data[,c("id",LM_fixedEffect_withoutTime_variables)]
  XX = as.matrix(X[!duplicated(X),-1])
  X = as.matrix(X[,LM_fixedEffect_withoutTime_variables])
  # variables with random effect in the longitudinal process
  Q=data[,LM_randomEffect_variables,drop=FALSE] # variables with random effect
  dNit = as.numeric(!is.na(data$Y))

  temp = as.data.frame(data) %>%
    group_by(id) %>%
    summarize(C = max(time),ni=sum(!is.na(Y)))
  C = rep(max(na.omit(data$time)),nrow(VV)) # censoring time
  ni=temp$ni # number of observations per subject
  id.unique = unique(data$id)
  id = data$id
  Time = data$time
  Time.unique = sort(unique(Time),decreasing = FALSE)
  y = data$Y
  # J = time.end
  
  gamma.est.fun = function(V,VV){
    # data,data_base,vp_pred
    m = nrow(VV)
    VV = as.matrix(VV)
    
    gamma_function = function(gamma){
      gamma = matrix(gamma,ncol=1)
      # compute the average value
      exp.gamma = exp(as.matrix(VV) %*% gamma)
      
      V.denom = rep(sum(na.omit(exp.gamma)),length(Time))
      V.numer = matrix(rep(apply(na.omit(as.vector(exp.gamma)) * VV,2,sum),length(Time)),nrow=length(Time),byrow = TRUE)

      V.weightedAve = V.numer/V.denom
      (obj = apply(V-V.weightedAve,2,mean))
    }
    
    gamma = rep(0.1,ncol(V))
    gamma.solve = nleqslv(x=gamma,fn=gamma_function,control=list(trace=1))
    (gamma.hat = gamma.solve$x)
    
    return(gamma.hat)
  }
  
  (gamma.hat = gamma.est.fun(V,VV))
  names(gamma.hat) = colnames(V)
  print(gamma.hat)
  
  # Lambda_estimation
  exp_gamma = as.vector(exp(VV %*%gamma.hat))
  denom_Lambda = rep(sum(exp_gamma),length(Time.unique))
  Time_freqtable = table(Time)
  Time_freqtable1 = Time_freqtable[match(as.character(Time.unique),names(Time_freqtable))]
  Lambda_t_est = as.numeric(Time_freqtable1/denom_Lambda)
  Lambda_C_est= rep(sum(Lambda_t_est[Time.unique<=C[1]]),length(C))
  
   # sigma^2 estimation
  exp_gamma_Lambda = exp_gamma*Lambda_C_est
  sigma_hat_sq = max( (sum(ni^2-ni-exp_gamma_Lambda^2)/sum(exp_gamma_Lambda^2)),0)
  print(sigma_hat_sq)

  # B estimation
  denom = rep(sum(ni/Lambda_C_est),length(Time))
  
  Exp_eta = ((1+ni*sigma_hat_sq)/(1+exp_gamma_Lambda*sigma_hat_sq)-1)
  id_freqtable = table(id)
  id_freqtable1 = id_freqtable[match(as.character(id.unique),names(id_freqtable))]
  Exp_eta_long = rep(Exp_eta, id_freqtable1)
  # if time is not included in Q
  if ("time"%in%colnames(Q)){
    # remove time column in Q
    Q_noTime = Q[,-which(colnames(Q)=="time"),drop=FALSE]
    Q0_noTime = cbind(1,Q_noTime)
    Bhat = Exp_eta_long*Q0_noTime
    Bhat_base0 = cbind(data$id,Bhat)
    Bhat_base = Bhat_base0[!duplicated(Bhat_base0),-1]
    numer_B1 = apply(Bhat_base,2,function(x){
      rep(sum(na.omit(ni*x/Lambda_C_est)),length(Time))
    })
    numer_BTime = unlist(lapply(Time,function(x){sum(na.omit(ni*x*Exp_eta/Lambda_C_est))} ))
    B_bar = cbind(numer_B1/denom, numer_BTime/denom)
    apply(na.omit(B_bar),2,sd)
  } else{
    Q0 = cbind(1,Q)
    Bhat = Exp_eta_long*Q0
    Bhat_base0 = cbind(data$id,Bhat)
    Bhat_base = Bhat_base0[!duplicated(Bhat_base0),-1]
    numer_B = apply(Bhat_base,2,function(x){
      rep(sum(na.omit(ni*x/Lambda_C_est)),length(Time))
    })
    B_bar = numer_B/denom
  }
  Bhat = Exp_eta_long*cbind(1,Q)

  # Q0 = cbind(1,Q)
  # Bhat_i = ((1+ni*sigma_hat_sq)/(1+exp_gamma_Lambda*sigma_hat_sq)-1)
  # id_freqtable = table(id)
  # id_freqtable1 = id_freqtable[match(as.character(id.unique),names(id_freqtable))]
  # Exp_eta_long = rep(Bhat_i, id_freqtable1)
  # Bhat = Exp_eta_long*Q0
  # Bhat_base0 = cbind(data$id,Bhat)
  # Bhat_base = Bhat_base0[!duplicated(Bhat_base0),-1]
  # # estimating longitudinal parameters
  # numer_B = apply(Bhat_base,2,function(x){
  #   rep(sum(na.omit(ni*x/Lambda_C_est)),length(Time))
  # })

  # X estimation
  # X_bar, XX has only time-invariant variables
  numer_X = apply(XX,2,function(x){
    rep(sum(ni*x/Lambda_C_est),length(Time))
  })  
  X_bar = numer_X/denom

  if (sigma_hat_sq != 0 ){
    design_M = as.matrix(cbind(X,Bhat))
    design_Mbar = as.matrix(cbind(X_bar, B_bar))
    Liang.function = function(beta){
      tmp1 = (design_M - design_Mbar)*as.vector(data$Y-design_M%*%beta)
      value = apply(na.omit(tmp1),2,sum)
      return(value)
    }
    beta = as.matrix(rep(0.1,ncol(design_M)),ncol=1)
    beta.solve = nleqslv(x=beta,fn=Liang.function,control=list(trace=1))
    (beta.hat = beta.solve$x[1:ncol(X)])
  }

  apply(na.omit(design_M),2,sd)

  if (sigma_hat_sq == 0 ){
    design_M = as.matrix(cbind(X))
    design_Mbar = as.matrix(cbind(X_bar))
    Liang.function = function(beta){
      tmp1 = (design_M - design_Mbar)*as.vector(data$Y-design_M%*%beta)
      value = apply(na.omit(tmp1),2,sum)
      return(value)
    }
    beta = as.matrix(rep(0,ncol(design_M)),ncol=1)
    beta.solve = nleqslv(x=beta,fn=Liang.function,control=list(trace=1))
    (beta.hat = beta.solve$x[1:ncol(X)])
  }
  names(beta.hat) = colnames(X)
  (beta.hat)

  results = list(gamma.hat = gamma.hat,
                 beta.hat = beta.hat)
  
  return(results)
}

joint_LY_realData = function(long_data){
  
  data = long_data
  # data$time = round(data$time,1)
  
  # variables in the visiting process
  V= data[,c("id",VPM_variables)] # row = n
  VV = as.matrix(V[!duplicated(V),-1]) # row=m
  V = as.matrix(V[,VPM_variables]) # variables affecting the visiting process
  # variables in the longitudinal process
  X = data[,c("id",LM_fixedEffect_withoutTime_variables)]
  XX = as.matrix(X[!duplicated(X),-1])
  X = as.matrix(X[,LM_fixedEffect_withoutTime_variables])

  temp = data %>%
    group_by(id) %>%
    summarize(C = max(time),ni=n())
  C = rep(max(na.omit(data$time)),nrow(VV)) # censoring time
  ni=temp$ni # number of observations per subject
  id.unique = unique(data$id)
  id = data$id
  Time = data$time
  Time.unique = sort(unique(Time),decreasing = FALSE)
  y = data$Y
  # J = time.end
  dNit = as.numeric(!is.na(data$Y))
  
 gamma.est.fun = function(V,VV){
    # data,data_base,vp_pred
    m = nrow(VV)
    VV = as.matrix(VV)
    
    gamma_function = function(gamma){
      gamma = matrix(gamma,ncol=1)
      # compute the average value
      exp.gamma = exp(as.matrix(VV) %*% gamma)
      
      # V.denom = sapply(Time, function(t){sum((t<=C)*exp.gamma)})
      # V.numer = apply(VV,2,function(x){
      #   sapply(Time,function(t){sum( (t<=C)*exp.gamma*x )})
      # })
      V.denom = rep(sum(na.omit(exp.gamma)),length(Time))
      V.numer = matrix(rep(apply(na.omit(as.vector(exp.gamma)) * VV,2,sum),length(Time)),nrow=length(Time),byrow = TRUE)

      V.weightedAve = V.numer/V.denom
      (obj = apply(V-V.weightedAve,2,mean))
    }
    
    gamma = rep(0,ncol(V))
    gamma.solve = nleqslv(x=gamma,fn=gamma_function,control=list(trace=1))
    (gamma.hat = gamma.solve$x)
    
    return(gamma.hat)
  }
  
  (gamma.hat = gamma.est.fun(V,VV))
  names(gamma.hat) = colnames(V)
  print(gamma.hat)
  
  # gamma.hat = c(0.5,0.5)
  exp.gamma = as.vector(exp(VV %*%gamma.hat))
  y.data = data.frame(id,Time,y)
  y.star.fun = function(t){
      y.data1 = y.data %>%
                mutate(time_diff = abs(Time-t)) %>%
                group_by(id)%>%
                summarize(y_star = y[which.min(time_diff)][1])
      # y.data1 = y.data %>%
      #           mutate(time_diff = abs(Time-t)) %>%
      #           group_by(id)%>%
      #           mutate(y_star = y[which.min(time_diff)][1])
      return(y.data1$y_star)
  } 

  length(Time.unique)
  y.bar.numer = sapply(Time.unique,function(u){ sum(na.omit(y.star.fun(u)*exp.gamma)) })
  y.bar.denom = rep(sum(exp.gamma),length(Time.unique))
  y.bar.ordered = y.bar.numer/y.bar.denom
  y.bar = y.bar.ordered[match(Time,Time.unique)]
  
  # I should not include time in the formula cuz time is handled in the baseline.
  X.bar.ordered = t(sapply(Time.unique,function(t){
    X.bar.numer = c(colSums(XX*exp.gamma*(t<=C)))
    X.bar.denom = sum(exp.gamma*(t<=C))
    X.bar.numer/X.bar.denom
  }))

  X.bar.denom = rep(sum(exp.gamma),length(Time.unique))
  # replicate X.numer.ave to match the length of Time.unique
  X.bar.numer = matrix(rep(apply(XX*exp.gamma,2,sum),length(Time.unique)),nrow=length(Time.unique),byrow = TRUE)
  X.bar.ordered = X.bar.numer/X.bar.denom
  X.bar = X.bar.ordered[match(Time,Time.unique),]

  ###### LY method ######
  LY.function = function(beta){
    res = as.vector((y-y.bar)-(X-X.bar) %*% beta)
    (obj = apply(X-X.bar,2,function(x){sum(na.omit(x*res))}))
    return(obj)
  }
  
  beta = matrix(rep(0, ncol(X)),ncol=1)
  beta_hat_LY.solve = nleqslv(x=beta,fn=LY.function,control=list(trace=1))
  (beta_hat_LY = beta_hat_LY.solve$x)
  names(beta_hat_LY) = colnames(X)

  results = list(gamma.hat = gamma.hat,
                 beta_hat_LY = beta_hat_LY)
  
  return(results)
}

joint_weightedGEE_realData = function(long_data){
  
  data = long_data
  # data$time = round(data$time,1)
  
  # variables in the visiting process
  V= data[,c("id",VPM_variables)] # row = n
  VV = as.matrix(V[!duplicated(V),-1]) # row=m
  V = as.matrix(V[,VPM_variables]) # variables affecting the visiting process
  # variables in the longitudinal process
  X = data[,c("id",LM_fixedEffect_withoutTime_variables)]
  XX = as.matrix(X[!duplicated(X),-1])
  X = as.matrix(X[,LM_fixedEffect_withoutTime_variables])

  temp = data %>%
    group_by(id) %>%
    summarize(C = max(time),ni=n())
  C = rep(max(na.omit(data$time)),nrow(VV)) # censoring time
  ni=temp$ni # number of observations per subject
  id.unique = unique(data$id)
  id = data$id
  Time = data$time
  Time.unique = sort(unique(Time),decreasing = FALSE)
  y = data$Y
  # J = time.end
  dNit = as.numeric(!is.na(data$Y))
  
 gamma.est.fun = function(V,VV){
    # data,data_base,vp_pred
    m = nrow(VV)
    VV = as.matrix(VV)
    
    gamma_function = function(gamma){
      gamma = matrix(gamma,ncol=1)
      # compute the average value
      exp.gamma = exp(as.matrix(VV) %*% gamma)
      
      # V.denom = sapply(Time, function(t){sum((t<=C)*exp.gamma)})
      # V.numer = apply(VV,2,function(x){
      #   sapply(Time,function(t){sum( (t<=C)*exp.gamma*x )})
      # })
      V.denom = rep(sum(na.omit(exp.gamma)),length(Time))
      V.numer = matrix(rep(apply(na.omit(as.vector(exp.gamma)) * VV,2,sum),length(Time)),nrow=length(Time),byrow = TRUE)

      V.weightedAve = V.numer/V.denom
      (obj = apply(V-V.weightedAve,2,mean))
    }
    
    gamma = rep(0,ncol(V))
    gamma.solve = nleqslv(x=gamma,fn=gamma_function,control=list(trace=1))
    (gamma.hat = gamma.solve$x)
    
    return(gamma.hat)
  }
  
  (gamma.hat = gamma.est.fun(V,VV))
  names(gamma.hat) = colnames(V)
  print(gamma.hat)
    
  ###### weightedGEE with the intercept and time ######
  weights = as.vector(exp(V %*%gamma.hat))
  X.GEE = as.matrix(data[,c(LM_fixedEffect_withTime_variables)])
  X.GEE = cbind(1,X.GEE)
  Y.function = function(beta){
    res = as.vector(1/weights*(y-(X.GEE %*% beta)))
    (obj = apply(X.GEE,2,function(x){sum(na.omit(x*res))}))
    return(obj)
  }
  beta = matrix(rep(0.1, ncol(X.GEE)),ncol=1)
  beta_hat_weightedGEE.solve = nleqslv(x=beta,fn=Y.function,control=list(trace=0))
  (beta_hat_weightedGEE = beta_hat_weightedGEE.solve$x)
  names(beta_hat_weightedGEE) = colnames(X.GEE)

  results = list(gamma.hat = gamma.hat,
                 beta_hat_weightedGEE = beta_hat_weightedGEE)
  
  return(results)
}

joint_combined_realData = function(long_data){
  
  data = long_data
  # data$time = round(data$time,1)
  
  # variables in the visiting process
  V= data[,c("id",VPM_variables)] # row = n
  VV = as.matrix(V[!duplicated(V),-1]) # row=m
  V = as.matrix(V[,VPM_variables]) # variables affecting the visiting process
  # variables in the longitudinal process
  X = data[,c("id",LM_fixedEffect_withoutTime_variables)]
  XX = as.matrix(X[!duplicated(X),-1])
  X = as.matrix(X[,LM_fixedEffect_withoutTime_variables])

  temp = data %>%
    group_by(id) %>%
    summarize(C = max(time),ni=n())
  C = rep(max(na.omit(data$time)),nrow(VV)) # censoring time
  ni=temp$ni # number of observations per subject
  id.unique = unique(data$id)
  id = data$id
  Time = data$time
  Time.unique = sort(unique(Time),decreasing = FALSE)
  y = data$Y
  # J = time.end
  dNit = as.numeric(!is.na(data$Y))
  
 gamma.est.fun = function(V,VV){
    # data,data_base,vp_pred
    m = nrow(VV)
    VV = as.matrix(VV)
    
    gamma_function = function(gamma){
      gamma = matrix(gamma,ncol=1)
      # compute the average value
      exp.gamma = exp(as.matrix(VV) %*% gamma)
      
      # V.denom = sapply(Time, function(t){sum((t<=C)*exp.gamma)})
      # V.numer = apply(VV,2,function(x){
      #   sapply(Time,function(t){sum( (t<=C)*exp.gamma*x )})
      # })
      V.denom = rep(sum(na.omit(exp.gamma)),length(Time))
      V.numer = matrix(rep(apply(na.omit(as.vector(exp.gamma)) * VV,2,sum),length(Time)),nrow=length(Time),byrow = TRUE)

      V.weightedAve = V.numer/V.denom
      (obj = apply(V-V.weightedAve,2,mean))
    }
    
    gamma = rep(0,ncol(V))
    gamma.solve = nleqslv(x=gamma,fn=gamma_function,control=list(trace=1))
    (gamma.hat = gamma.solve$x)
    
    return(gamma.hat)
  }
  
  (gamma.hat = gamma.est.fun(V,VV))
  names(gamma.hat) = colnames(V)
  print(gamma.hat)
  
  # gamma.hat = c(0.5,0.5)
  exp.gamma = as.vector(exp(VV %*%gamma.hat))
  y.data = data.frame(id,Time,y)
  y.star.fun = function(t){
      y.data1 = y.data %>%
                mutate(time_diff = abs(Time-t)) %>%
                group_by(id)%>%
                summarize(y_star = y[which.min(time_diff)][1])
      # y.data1 = y.data %>%
      #           mutate(time_diff = abs(Time-t)) %>%
      #           group_by(id)%>%
      #           mutate(y_star = y[which.min(time_diff)][1])
      return(y.data1$y_star)
  } 

  length(Time.unique)
  y.bar.numer = sapply(Time.unique,function(u){ sum(na.omit(y.star.fun(u)*exp.gamma)) })
  y.bar.denom = rep(sum(exp.gamma),length(Time.unique))
  y.bar.ordered = y.bar.numer/y.bar.denom
  y.bar = y.bar.ordered[match(Time,Time.unique)]
  
  
  # I should not include time in the formula cuz time is handled in the baseline.
  X.bar.ordered = t(sapply(Time.unique,function(t){
    X.bar.numer = c(colSums(XX*exp.gamma*(t<=C)))
    X.bar.denom = sum(exp.gamma*(t<=C))
    X.bar.numer/X.bar.denom
  }))

  X.bar.denom = rep(sum(exp.gamma),length(Time.unique))
  # replicate X.numer.ave to match the length of Time.unique
  X.bar.numer = matrix(rep(apply(XX*exp.gamma,2,sum),length(Time.unique)),nrow=length(Time.unique),byrow = TRUE)
  X.bar.ordered = X.bar.numer/X.bar.denom
  X.bar = X.bar.ordered[match(Time,Time.unique),]

  ###### LY method ######
  LY.function = function(beta){
    res = as.vector((y-y.bar)-(X-X.bar) %*% beta)
    (obj = apply(X-X.bar,2,function(x){sum(na.omit(x*res))}))
    return(obj)
  }
  
  beta = matrix(rep(0, ncol(X)),ncol=1)
  beta_hat_LY.solve = nleqslv(x=beta,fn=LY.function,control=list(trace=1))
  (beta_hat_LY = beta_hat_LY.solve$x)
  names(beta_hat_LY) = colnames(X)
  
  ###### weightedGEE with the intercept and time ######
  weights = as.vector(exp(V %*%gamma.hat))
  X.GEE = as.matrix(data[,c(LM_fixedEffect_withTime_variables)])
  X.GEE = cbind(1,X.GEE)
  Y.function = function(beta){
    res = as.vector(1/weights*(y-(X.GEE %*% beta)))
    (obj = apply(X.GEE,2,function(x){sum(na.omit(x*res))}))
    return(obj)
  }
  beta = matrix(rep(0.1, ncol(X.GEE)),ncol=1)
  beta_hat_weightedGEE.solve = nleqslv(x=beta,fn=Y.function,control=list(trace=0))
  (beta_hat_weightedGEE = beta_hat_weightedGEE.solve$x)
  names(beta_hat_weightedGEE) = colnames(X.GEE)

  results = list(gamma.hat = gamma.hat,
                 beta_hat_LY = beta_hat_LY,
                 beta_hat_weightedGEE = beta_hat_weightedGEE)
  
  return(results)
}

merlin_fun = function(long_data){

  # calculate the inter-visit time
  long_data1 = long_data %>%
      group_by(id) %>%
      mutate(time0 = lag(time))%>%
      mutate(time0 = ifelse(is.na(time0),0,time0))%>%
      mutate(s = time - time0) %>%
      # create a column "delta" where the last observation for each id is 0 and others are 1
      mutate(delta = ifelse(row_number() == n(), 0, 1)) |> as.data.frame()

  control_list = list(intmethod=c("ghermite"),maxit=1000)
  # control_list = list(intmethod=c("sobol"),maxit=1000) # this optimization doesn't work
  # control_list = list(intmethod=c("halton"),maxit=1000) # this optimization doesn't work

  # Use Gaspirini's approach, use the "ghermite" (default) optimization method
  merlin_LongModel = as.formula(paste0("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+"),"+M1[id]"))
  merlin_VisitModel = as.formula(paste0("Surv(s, delta)~",paste(VPM_variables,collapse = "+"),"+M1[id]*1"))
  merlin_model = merlin(model = list(merlin_LongModel,merlin_VisitModel), 
        family = c("gaussian","exponential"), 
        levels = "id",
        timevar = c("time","s"),
        data = long_data1,
        control=control_list)
  print(summary(merlin_model))
  merlin_summary = summary(merlin_model)
  para_est_LM = merlin_summary$coefftable[1:(length(LM_fixedEffect_withTime_variables)+3),]
  para_est_VPM = merlin_summary$coefftable[-(1:(length(LM_fixedEffect_withTime_variables)+3)),]
  beta_hat = para_est_LM[c("_cons",LM_fixedEffect_withTime_variables),"Estimate"]
  names(beta_hat)[1] = "(Intercept)"
  beta_sd = para_est_LM[c("_cons",LM_fixedEffect_withTime_variables),"Std. Error"]
  names(beta_sd)[1] = "(Intercept)"

  return(list(beta_hat = beta_hat,
              beta_sd = beta_sd))
}

merlin_withvi_fun = function(long_data){

  # calculate the inter-visit time
  long_data1 = long_data %>%
      group_by(id) %>%
      mutate(time0 = lag(time))%>%
      mutate(time0 = ifelse(is.na(time0),0,time0))%>%
      mutate(s = time - time0) %>%
      # create a column "delta" where the last observation for each id is 0 and others are 1
      mutate(delta = ifelse(row_number() == n(), 0, 1)) |> as.data.frame()

  control_list = list(intmethod=c("ghermite"),maxit=1000)
  # control_list = list(intmethod=c("sobol"),maxit=1000) # this optimization doesn't work
  # control_list = list(intmethod=c("halton"),maxit=1000) # this optimization doesn't work

  # Use Gaspirini's approach, use the "ghermite" (default) optimization method
  merlin_LongModel = as.formula(paste0("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+"),"+M1[id]+M2[id]*1"))
  merlin_VisitModel = as.formula(paste0("Surv(s, delta)~",paste(VPM_variables,collapse = "+"),"+M1[id]*1"))
  merlin_model = merlin(model = list(merlin_LongModel,merlin_VisitModel), 
        family = c("gaussian","exponential"), 
        levels = "id",
        timevar = c("time","s"),
        data = long_data1,
        control=control_list)
  print(summary(merlin_model))
  merlin_summary = summary(merlin_model)
  para_est_LM = merlin_summary$coefftable[1:(length(LM_fixedEffect_withTime_variables)+3),]
  para_est_VPM = merlin_summary$coefftable[-(1:(length(LM_fixedEffect_withTime_variables)+3)),]
  beta_hat = para_est_LM[c("_cons",LM_fixedEffect_withTime_variables),"Estimate"]
  names(beta_hat)[1] = "(Intercept)"
  beta_sd = para_est_LM[c("_cons",LM_fixedEffect_withTime_variables),"Std. Error"]
  names(beta_sd)[1] = "(Intercept)"

  return(list(beta_hat = beta_hat,
              beta_sd = beta_sd))
}


imputation_realData_year = function(long_data){

    data = long_data
    # monthly-grouped data
    data = data %>%
        mutate(time_6month = time*2,time_roundup_6month = ceiling(time_6month))%>%
        group_by(id,time_roundup_6month) %>%
        mutate(Y_mean = mean(Y,na.rm = TRUE),Age_mean = mean(Age,na.rm=TRUE)) %>%
        select(id,Y_mean,Age_mean,Sex,time_roundup_6month,SNP)
    data1 = data[!duplicated(data),]
    colnames(data1) = c("id","Y","Age","Sex","time","SNP")

    # create a new dataset with the same columns in data1. For a given id, the time is from 1 to max_time. If there is no value of Y at a given time, then Y is NA.
    max_time <- max(data1$time, na.rm = TRUE)
    id_all = rep(unique(data1$id),each = max_time)
    time_all = rep(as.numeric(1:max_time),length(unique(data1$id)))
    Age_all = rep(unique(data1[,c("id","Age")])$Age,each=max_time)
    Sex_all = rep(unique(data1[,c("id","Sex")])$Sex,each=max_time)
    SNP_all = rep(unique(data1[,c("id","SNP")])$SNP,each=max_time)
    data2 = data.frame(id = id_all,time = time_all,Sex = Sex_all, Age = Age_all, SNP = SNP_all)
    
    # join data1 to data2 by id and time
    data3 = left_join(data2,data1,by = c("id","time","Age","Sex","SNP"))
    data3$time = data3$time/2

    df_full2 = data3
    dim(df_full2)
    head(df_full2)
    cor(na.omit(df_full2[,-1]))

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

    # # multilevel analysis 
    # data_imp = complete(imp1)
    # data_imp$Sex = as.factor(data_imp$Sex)
    # lme_model_formula = paste("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+"),"+(1+",paste(LM_randomEffect_variables,collapse = "+"),"|id)")
    # lme_model = lmer(lme_model_formula,data=data_imp,REML = TRUE,
    #                 control=lmerControl(optCtrl=list(optimizer ='optimx', optCtrl=list(method='L-BFGS-B'),maxfun=50000)))
    # (beta_hat = summary(lme_model)$coef[,1])
    lme_imp = function(data_imp){
      data_imp$Sex = as.factor(data_imp$Sex)
      lme_model_formula = paste("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+"),"+(1+",paste(LM_randomEffect_variables,collapse = "+"),"|id)")
      lme_model = lmer(lme_model_formula,data=data_imp,REML = TRUE,
                      control=lmerControl(optCtrl=list(optimizer ='optimx', optCtrl=list(method='L-BFGS-B'),maxfun=50000)))
      (beta_hat = summary(lme_model)$coef[,1])
    }

    fit = lapply(1:5, function(i) lme_imp(complete(imp1, action=i)))
    beta_hat = sapply(seq_along(fit[[1]]), function(i) mean(sapply(fit, `[`, i)))
    names(beta_hat) = names(fit[[1]])


    # lme_model_fixed_formula = as.formula(paste("Y ~",paste(LM_fixedEffect_withTime_variables,collapse = "+")))
    # lme_model_random_formula = as.formula(paste("~1+",paste(LM_randomEffect_variables,collapse = "+"),"|id"))
    # lme_model = lme(lme_model_fixed_formula,data=data_imp,random=lme_model_random_formula, method="REML",control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))

    # (beta_hat = coef(summary(lme_model))[,1])

  return(beta_hat)
}






