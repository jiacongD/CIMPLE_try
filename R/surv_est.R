cox_fun <- function(long_data, surv_data) {
  long_data2 <- long_data %>%
    group_by(id) %>%
    mutate(time0 = ifelse(is.na(lag(time)), 0, lag(time))) %>%
    filter(time != 0) %>%
    filter(time0 < time)
  long_data2$d0 <- surv_data$d[match(long_data2$id, surv_data$id)]
  long_data2 <- na.omit(long_data2) %>%
    group_by(id) %>%
    mutate(d = ifelse(time < max(time, na.rm = TRUE), 0, d0))
  model_formula <- as.formula(paste("Surv(time0, time, d) ~ ", paste(SM_variables, collapse = "+")))
  model <- coxph(model_formula, data = long_data2)
  alpha.hat <- summary(model)$coef[, 1]

  return(alpha.hat)
}


surv_est_1 <- function(long_data, surv_data, model, LM_fixedEffect_withTime_variables, LM_randomEffect_variables, SM_base_variables) {
  if (model == "JMLD") {
    # remove patients with no Y measurement
    zeroRecords_id <- long_data[is.na(long_data$Y), "id"]
    long_data <- long_data[!long_data$id %in% zeroRecords_id, ]
    surv_data <- surv_data[!surv_data$id %in% zeroRecords_id, ]

    # longitudinal submodel, try nlminb optimizer first, if there is an error, then use optim optimizer
    lmeFit_fixedformula <- as.formula(paste("Y ~ ", paste(LM_fixedEffect_withTime_variables, collapse = "+")))
    lmeFit_randomformula <- as.formula(paste("~", paste(LM_randomEffect_variables, collapse = "+"), "|id"))
    control_optim <- lmeControl(msMaxIter = 1000, msMaxEval = 1000, opt = "optim")
    control_nlminb <- lmeControl(msMaxIter = 1000, msMaxEval = 1000, opt = "nlminb")
    lmeFit <- tryCatch(
      {
        lmeFit <- lme(lmeFit_fixedformula, random = lmeFit_randomformula, data = long_data, control = control_optim)
      },
      error = function(e) {
        print(paste0("Error with optim: ", e))
        lmeFit <- lme(lmeFit_fixedformula, random = lmeFit_randomformula, data = long_data, control = control_nlminb)
      }
    )
    # print(lmeFit)

    # survival submodel
    coxFit_formula <- as.formula(paste("Surv(D, d) ~ ", paste(SM_base_variables, collapse = "+")))
    coxFit <- coxph(coxFit_formula, data = surv_data, x = TRUE)
    # summary(coxFit)

    # jointFit = JMbayes2::jm(coxFit,lmeFit,time_var="time", n_chains=1L,n_iter=11000L,n_burnin=1000L)
    jointFit <- JMbayes2::jm(coxFit, lmeFit, time_var = "time", n_chains = 1L)

    # print(summary(jointFit))

    surv_proc <- unlist(coef(jointFit)) # change the name of the output to alpha
    long_proc <- unlist(fixef(jointFit))

    results <- list(
      "long_proc" = long_proc,
      "surv_proc" = surv_proc
    )

    return(results)
  } else if (model == "VA_JMLD") {
    # remove patients with no Y measurement
    zeroRecords_id <- long_data[is.na(long_data$Y), "id"]
    long_data <- long_data[!long_data$id %in% zeroRecords_id, ]
    surv_data <- surv_data[!surv_data$id %in% zeroRecords_id, ]

    # add Ni(t) as a predictor
    long_data$Ni <- NA
    for (i in 1:nrow(long_data)) {
      id <- long_data$id[i]
      time <- long_data$time[i]
      long_data$Ni[i] <- sum(!is.na(long_data$Y[long_data$id == id & long_data$time <= time]))
    }

    # longitudinal submodel, try nlminb optimizer first, if there is an error, then use optim optimizer
    lmeFit_fixedformula <- as.formula(paste("Y ~ ", paste(LM_fixedEffect_withTime_variables, collapse = "+"), "+Ni"))
    lmeFit_randomformula <- as.formula(paste("~", paste(LM_randomEffect_variables, collapse = "+"), "|id"))
    control_optim <- lmeControl(msMaxIter = 1000, msMaxEval = 1000, opt = "optim")
    control_nlminb <- lmeControl(msMaxIter = 1000, msMaxEval = 1000, opt = "nlminb")
    lmeFit <- tryCatch(
      {
        lmeFit <- lme(lmeFit_fixedformula, random = lmeFit_randomformula, data = long_data, control = control_optim)
      },
      error = function(e) {
        print(paste0("Error with optim: ", e))
        lmeFit <- lme(lmeFit_fixedformula, random = lmeFit_randomformula, data = long_data, control = control_nlminb)
      }
    )
    print(lmeFit)

    # survival submodel
    coxFit_formula <- as.formula(paste("Surv(D, d) ~ ", paste(SM_base_variables, collapse = "+")))
    coxFit <- coxph(coxFit_formula, data = surv_data, x = TRUE)
    summary(coxFit)

    # jointFit = JMbayes2::jm(coxFit,lmeFit,time_var="time", n_chains=1L,n_iter=11000L,n_burnin=1000L)
    jointFit <- JMbayes2::jm(coxFit, lmeFit, time_var = "time", n_chains = 1L)

    print(summary(jointFit))

    surv_proc <- unlist(coef(jointFit))
    long_proc <- unlist(fixef(jointFit))

    results <- list(
      "long_proc" = long_proc,
      "surv_proc" = surv_proc
    )

    return(results)
  }
}


surv_est_2 <- function(long_data, surv_data, model, imp_time_factor = NULL) {
  if (model == "Imputation_Cox") {
    data <- long_data

    if (is.null(imp_time_factor)) {
      imp_time_factor <- 1
      data3 <- data
    } else {
      # If the imp_time_factor is specified and is not 1
      # imp_time_factor = 0.5
      # group the data based on the time factor
      data <- data %>%
        mutate(time_new = ceiling(time / imp_time_factor)) %>%
        group_by(id, time_new) %>%
        mutate(Y_mean = mean(Y, na.rm = TRUE)) %>%
        # mutate(Y_mean = ifelse(is.numeric(Y_mean),Y_mean,NA)) %>%
        dplyr::select(id, Y_mean, time_new, setdiff(names(data), c("id", "Y", "time"))) %>%
        ungroup()
      colnames(data)[2:3] <- c("Y", "time")

      # create a new dataset with the same columns in data1. For a given id, the time is from 1 to max_time. If there is no value of Y at a given time, then Y is NA.
      max_time <- max(data$time, na.rm = TRUE)
      id_all <- rep(unique(data$id), each = max_time)
      time_all <- rep(as.numeric(1:max_time), length(unique(data$id)))
      data_base <- data %>%
        dplyr::select(-Y, -time) %>%
        unique()
      data2 <- data_base %>%
        slice(rep(1:n(), each = max_time)) %>%
        group_by(id) %>%
        mutate(time = rep(1:max_time))
      # delete time after the disease time
      data2$D <- surv_data$D[match(data2$id, surv_data$id)]
      data2 <- data2[data2$time <= ceiling(data2$D / imp_time_factor), setdiff(names(data2), c("Y"))]

      # join the two datasets
      data3 <- left_join(data2, data, by = setdiff(names(data), c("Y")))
      data3$time <- data3$time * imp_time_factor # convert time back to the original scale
    }

    print(paste0("imp factor: ", imp_time_factor))

    df_full <- data3
    dim(df_full)

    #### Start imputation ####
    # empty imputation
    imp0 <- mice(as.matrix(df_full), maxit = 0)
    predM <- imp0$predictorMatrix
    impM <- imp0$method

    # specify predictor matrix and method
    predM1 <- predM
    predM1["Y", "id"] <- -2
    predM1["Y", LM_fixedEffect_withTime_variables] <- 1 # fixed x effects imputation
    impM1 <- impM
    impM1["Y"] <- "2l.lmer"

    # multilevel imputation
    imp1 <- mice(as.matrix(df_full),
      m = 5,
      predictorMatrix = predM1, method = impM1, maxit = 5
    )

    # fit the cox-ph model
    coxph_imp <- function(data_imp) {
      # cox-ph model
      long_data_imp2 <- data_imp %>%
        group_by(id) %>%
        mutate(time0 = ifelse(is.na(lag(time)), 0, lag(time)))
      long_data_imp2$d0 <- surv_data$d[match(long_data_imp2$id, surv_data$id)]
      long_data_imp2 <- na.omit(long_data_imp2) %>%
        group_by(id) %>%
        mutate(d = ifelse(time < max(time, na.rm = TRUE), 0, d0)) %>%
        filter(time0 < time)
      model_formula <- as.formula(paste("Surv(time0, time, d) ~ ", paste(SM_variables, collapse = "+")))
      model <- coxph(model_formula, data = long_data_imp2)
      (alpha.hat <- summary(model)$coef[, 1])
    }
    fit <- lapply(1:5, function(i) coxph_imp(complete(imp1, action = i)))
    alpha.hat <- sapply(seq_along(fit[[1]]), function(i) mean(sapply(fit, `[`, i)))
    names(alpha.hat) <- names(fit[[1]])

    return(alpha.hat)
  } else if (model == "VAImputation_Cox") {
    data <- long_data

    if (is.null(imp_time_factor)) {
      imp_time_factor <- 1
      data3 <- data
    } else {
      # If the imp_time_factor is specified and is not 1
      # group the data based on the time factor
      data <- data %>%
        mutate(time_new = ceiling(time / imp_time_factor)) %>%
        group_by(id, time_new) %>%
        mutate(Y_mean = mean(Y, na.rm = TRUE)) %>%
        # mutate(Y_mean = ifelse(is.numeric(Y_mean),Y_mean,NA)) %>%
        dplyr::select(id, Y_mean, time_new, setdiff(names(data), c("id", "Y", "time"))) %>%
        ungroup()
      colnames(data)[2:3] <- c("Y", "time")

      # create a new dataset with the same columns in data1. For a given id, the time is from 1 to max_time. If there is no value of Y at a given time, then Y is NA.
      max_time <- max(data$time, na.rm = TRUE)
      id_all <- rep(unique(data$id), each = max_time)
      time_all <- rep(as.numeric(1:max_time), length(unique(data$id)))
      data_base <- data %>%
        dplyr::select(-Y, -time) %>%
        unique()
      data2 <- data_base %>%
        slice(rep(1:n(), each = max_time)) %>%
        group_by(id) %>%
        mutate(time = rep(1:max_time))
      # delete time after the disease time
      data2$D <- surv_data$D[match(data2$id, surv_data$id)]
      data2 <- data2[data2$time <= ceiling(data2$D / imp_time_factor), setdiff(names(data2), c("Y"))]

      # join the two datasets
      data3 <- left_join(data2, data, by = setdiff(names(data), c("Y")))
      data3$time <- data3$time * imp_time_factor # convert time back to the original scale
    }

    # insert the predictor Ni(t)
    df_full <- data3
    df_full$Ni <- NA
    for (i in 1:nrow(df_full)) {
      id <- df_full$id[i]
      time <- df_full$time[i]
      df_full$Ni[i] <- sum(!is.na(df_full$Y[df_full$id == id & df_full$time <= time]))
    }
    df_full <- df_full %>%
      as.data.frame() %>%
      dplyr::select(all_of(c(colnames(data), "Ni")))
    head(df_full)

    #### Start imputation ####
    # empty imputation
    imp0 <- mice(as.matrix(df_full), maxit = 0)
    predM <- imp0$predictorMatrix
    impM <- imp0$method

    # specify predictor matrix and method
    predM1 <- predM
    predM1["Y", "id"] <- -2
    predM1["Y", c(LM_fixedEffect_withTime_variables, "Ni")] <- 1 # fixed x effects imputation
    impM1 <- impM
    impM1["Y"] <- "2l.lmer"

    # multilevel imputation
    imp1 <- mice(as.matrix(df_full),
      m = 5,
      predictorMatrix = predM1, method = impM1, maxit = 5
    )

    # # cox-ph model
    # long_data_imp = complete(imp1)
    # long_data_imp2 = long_data_imp %>%
    #   group_by(id) %>%
    #   mutate(time0 = ifelse(is.na(lag(time)),0,lag(time) ))
    # long_data_imp2$d0 = surv_data$d[match(long_data_imp2$id,surv_data$id)]
    # long_data_imp2 = na.omit(long_data_imp2) %>%
    #             group_by(id) %>%
    #             mutate(d=ifelse(time<max(time,na.rm = TRUE),0,d0))

    # model_formula = as.formula(paste("Surv(time0, time, d) ~ ",paste(SM_variables,collapse="+")))
    # model = coxph(model_formula, data = long_data_imp2)
    # (alpha.hat = summary(model)$coef[,1])

    # fit the cox-ph model
    coxph_imp <- function(data_imp) {
      # cox-ph model
      long_data_imp2 <- data_imp %>%
        group_by(id) %>%
        mutate(time0 = ifelse(is.na(lag(time)), 0, lag(time)))
      long_data_imp2$d0 <- surv_data$d[match(long_data_imp2$id, surv_data$id)]
      long_data_imp2 <- na.omit(long_data_imp2) %>%
        group_by(id) %>%
        filter(time0 < time) %>%
        mutate(d = ifelse(time < max(time, na.rm = TRUE), 0, d0))
      model_formula <- as.formula(paste("Surv(time0, time, d) ~ ", paste(SM_variables, collapse = "+")))
      model <- coxph(model_formula, data = long_data_imp2)
      (alpha.hat <- summary(model)$coef[, 1])
    }
    fit <- lapply(1:5, function(i) coxph_imp(complete(imp1, action = i)))
    alpha.hat <- sapply(seq_along(fit[[1]]), function(i) mean(sapply(fit, `[`, i)))
    names(alpha.hat) <- names(fit[[1]])

    return(alpha.hat)
  }
}
