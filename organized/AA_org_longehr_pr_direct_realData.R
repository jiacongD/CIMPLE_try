library(survival)
library(nlme)
library(dplyr)
library(Matrix)
library(merlin)

# source("~/longitudinalEHR/merlin/selfcodes/gauher.R")

sum_rmna = function(x){
    sum(x[!is.na(x)])
}

likelihood_longehr = function(long_data){
    long_data$time = ifelse(is.na(long_data$time),5,long_data$time)
    long_data = long_data %>% filter(time!=0)
    
    # pre-process the data
    vp_data_censor = long_data %>% 
        group_by(id) %>% 
        # summarise(Z=first(Z),X=first(X),time=5)
        summarise(Age = first(Age),Sex=first(Sex),SNP = first(SNP),time=5)
    vp_data1 = long_data[,c("id",VPM_variables,"time")]
    vp_data = merge(vp_data1,vp_data_censor,by=c("id",VPM_variables,"time"),all=TRUE)
    vp_data = vp_data[order(vp_data$id,vp_data$time),] %>%
            group_by(id) %>%
            mutate(d = ifelse(time==time.end,0,1))%>%
            mutate(time0 = lag(time))%>%
            mutate(time0 = ifelse(is.na(time0),0,time0),
                    s = time-time0) %>%
                    dplyr::select(id,all_of(VPM_variables),s,d) %>%
            filter(s!=0) |>na.omit()
    long_data_censor = long_data %>% 
        group_by(id) %>% 
        # summarise(Z=first(Z),X=first(X),time=60,Y=NA)
        summarise(Age = first(Age),Sex=first(Sex),SNP = first(SNP),time=5,Y=NA)
    long_data = merge(long_data,long_data_censor,by=c("id",LM_fixedEffect_withTime_variables),all=TRUE)
    long_data = long_data[order(long_data$id,long_data$time),]
    long_data = long_data[which(!is.na(long_data$time)),] %>%
            select(id,all_of(LM_fixedEffect_withTime_variables),Y.x) %>%
            rename(Y=Y.x)
    
    # reformat the longitudinal data
    dim(long_data)
    dim(vp_data)

    # fit an lme object
    # lme_model_fixed_formula = as.formula(paste("Y~Z+X+time"))
    # lme_model_random_formula = as.formula(paste("~1|id"))
    lme_model_fixed_formula = as.formula(paste0("Y~",paste0(LM_fixedEffect_withTime_variables,collapse="+")))
    lme_model_random_formula = as.formula(paste("~1|id"))
    lme_model = lme(lme_model_fixed_formula,data=na.omit(long_data),random=lme_model_random_formula, method="ML",control = lmeControl(opt='optim',msMaxIter = 5000, msMaxEval = 5000))
    lmeObject = lme_model

    control=NULL
    control$verbose = TRUE
    control$tol1 = 1e-3 # control zero division
    control$tol2 = 1e-3 # tolerance rate for relative theta change
    control$tol3 = 1e-3 # tolerance rate for log-likelihood change

    initial_values = NULL
    initial_values$betas = coef(summary(lme_model))[,1]
    initial_values$sigma2_e = (lme_model$sigma)^2
    initial_values$sigma2_b = 0.5^2
    initial_values$rho = 0

    # vp_model = survreg(Surv(s, d) ~ Z + X, data=vp_data, dist="exponential") 
    vp_model_formula = as.formula(paste0("Surv(s, d)~",paste(VPM_variables,collapse="+")))
    vp_model = survreg(vp_model_formula, data=vp_data, dist="exponential") 
    # Note that when setting dist="exponential", the link function is "log" of the expected of the survival time, not the rate.
    # therefore, I need to take the negative of the coefficients
    summary(vp_model)
    initial_values$gammas = -coef(summary(vp_model))

    # input
    long_data = long_data
    vp_data = vp_data
    long_variable = LM_fixedEffect_withoutTime_variables
    vp_variable = VPM_variables   
    time_variable = "time"
    # control 
    control$GHk=10
    control$typeGH = "simple"
    control$maxiter = 150

    y=long_data$Y
    X = cbind(1,long_data[,c(LM_fixedEffect_withTime_variables)]) |>as.matrix()
    s = vp_data$s
    d = vp_data$d
    W = cbind(1,vp_data[,c(VPM_variables)]) |>as.matrix()
    id = long_data$id

    # sample size settings
    ncx = ncol(X)
    ncw = if (is.null(W)) 0 else ncol(W)
    N = length(y)
    n = length(unique(long_data$id))
    ni = as.vector(tapply(id, id, length))

    # crossproducts and others
    ncz = 1

    # Gauss-Hermite quadrature rule components
    gh_rule = statmod::gauss.quad(control$GHk)
    b = gh_rule$nodes
    wGH = gh_rule$weights
    b = as.matrix(expand.grid(rep(list(b), ncz)))
    wGH = as.matrix(expand.grid(rep(list(wGH), ncz)))  
    wGH = 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
    b = sqrt(2) * b

    # initial values
    betas = as.vector(initial_values$betas)
    sigma2_e = initial_values$sigma2_e
    gammas = as.vector(initial_values$gammas)
    sigma2_b = initial_values$sigma2_b
    rho = initial_values$rho

    ############# direct solve the score function ##########
    iter = control$maxiter
    Y_mat = matrix(0, iter + 1, ncx + 1)
    S_mat <- matrix(0, iter + 1, ncw )
    B_mat = matrix(0, iter + 1, ncz)
    Rho_mat = matrix(0, iter + 1, 1)
    lgLik = numeric(iter + 1)
    conv <- TRUE

    # Start the iteration
    it=1
    print(initial_values)

    for(it in 1:iter){
        # save parameter values in matrix
        Y_mat[it, ] <- c(betas, sigma2_e)
        S_mat[it, ] <- c(gammas)
        B_mat[it,] <- sigma2_b
        Rho_mat[it] <- rho

        bN = rep(1,N) %*% t(b)
        bn = rep(1,n) %*% t(b)

        # linear predictors
        gamma_eval = function(gammas){
            rho_b = rep(rho,N) %*% t(b)
            # longitudinal process
            eta_yx = as.vector(X %*% betas)
            mu_y = eta_yx+rho_b
            log_p_yb_N = dnorm(y, mu_y, sqrt(sigma2_e), log =TRUE)
            log_p_yb = rowsum(log_p_yb_N, id,na.rm=TRUE) 

            # visiting process
            eta_sw = as.vector(W %*% gammas)
            mu_s = eta_sw+bN
            log_p_sb_N = mu_s*d-exp(mu_s)*s
            log_p_sb = rowsum(log_p_sb_N, id)

            # random effect
            log_p_b_i = matrix(dnorm(b, rep(0, ncz), sd=sqrt(sigma2_b), log=TRUE),nrow=1)
            log_p_b = apply(log_p_b_i,2,function(x){rep(x,n)})
            
            # joint probability
            p_ysb = exp(log_p_yb+log_p_sb+log_p_b)
            p_ys = c(p_ysb %*% wGH)
            p_ys[which(p_ys==0)] = 1e-10
            p_bys = p_ysb/p_ys

            # partial derivative
            partial_gamma = d-exp(as.vector(W%*%gammas)+bN)*s
            partial_gamma = rowsum(partial_gamma, id)

            # take the product of each columns of p_bs and partial_gamma
            W_unique = unique(cbind(id,W))[,-1]
            (f = crossprod(W_unique,(partial_gamma*p_bys) %*% wGH))
            return(f)
        }
        gammas_solve = nleqslv::nleqslv(x=initial_values$gammas,gamma_eval,control = list(trace=0))
        (gammas_new = gammas_solve$x)

        # longitudinal model
        beta_eval = function(betasrho){
            betas = betasrho[1:ncx]
            rho = betasrho[ncx+1]
            rho_b = rep(rho,N) %*% t(b)

            # longitudinal process
            eta_yx = as.vector(X %*% betas)
            mu_y = eta_yx+rho_b
            log_p_yb_N = dnorm(y, mu_y, sqrt(sigma2_e), log =TRUE)
            log_p_yb = rowsum(log_p_yb_N, id,na.rm=TRUE) 

            # visiting process
            eta_sw = as.vector(W %*% gammas_new)
            mu_s = eta_sw+bN
            log_p_sb_N = mu_s*d-exp(mu_s)*s
            log_p_sb = rowsum(log_p_sb_N, id)

            # random effect
            log_p_b_i = matrix(dnorm(b, rep(0, ncz), sd=sqrt(sigma2_b), log=TRUE),nrow=1)
            log_p_b = apply(log_p_b_i,2,function(x){rep(x,n)})
            
            # joint probability
            p_ysb = exp(log_p_yb+log_p_sb+log_p_b)
            p_ys = c(p_ysb %*% wGH)
            p_ys[which(p_ys==0)] = 1e-10
            p_bys = p_ysb/p_ys

            p_bys_N = do.call(rbind, lapply(1:nrow(p_bys), 
                    function(i) matrix(rep(p_bys[i, ], each = ni[i]), nrow = ni[i], byrow = FALSE)))  

            # partial derivatives
            res = y-mu_y
            temp = as.vector((-res*p_bys_N/sigma2_e)%*%wGH)
            score_betas = colSums(X*temp,na.rm=TRUE)

            temp1 = as.vector((-res*bN*p_bys_N/sigma2_e)%*%wGH)
            score_rho = sum(temp1,na.rm=TRUE)

            (f = c(score_betas,score_rho))
            return(f)
        }

        betasrho_solve_fun = function(betas,rho){
            initial_guesses <- list(c(betas,rho),c(betas,0),c(betas,0.05),c(betas,-0.05),c(betas,1.5))

            results <- lapply(initial_guesses, function(guess) {
                tryCatch({
                    nleqslv::nleqslv(guess, beta_eval, method = "Newton")
                }, error = function(e) {
                    NULL
                })
            })

            # Filter successful results
            successful_results <- Filter(function(res) !is.null(res) && res$termcd == 1, results)

            # Choose the best result (e.g., with the lowest residual)
            if (length(successful_results) > 0) {
                best_result <- successful_results[[which.min(sapply(successful_results, function(res) sum(res$fvec^2)))]]
            } else {
                best_result <- NULL
                print("No successful solution found.")
            }

            return(best_result)
        }

        betasrho = c(betas,rho)
        # betasrho_solve = nleqslv::nleqslv(x=betasrho,beta_eval,control = list(trace=1))
        # betasrho_solve
        betasrho_solve = betasrho_solve_fun(betas,rho)
        if(is.null(betasrho_solve)){
            betas = NA
            stop("\n\nnot converged!\n")
        } else {
            betasrho_new = betasrho_solve$x
        }
        # (betasrho_new = betasrho_solve$x)
        betas_new = betasrho_new[1:ncx]
        rho_new = betasrho_new[ncx+1]

        # update the variance components
        # longitudinal process
        eta_yx = as.vector(X %*% betas_new)
        mu_y = eta_yx+rep(rho_new,N) %*% t(b)
        log_p_yb_N = dnorm(y, mu_y, sqrt(sigma2_e), log =TRUE)
        log_p_yb = rowsum(log_p_yb_N, id,na.rm=TRUE) 

        # visiting process
        eta_sw = as.vector(W %*% gammas_new)
        mu_s = eta_sw+bN
        log_p_sb_N = mu_s*d-exp(mu_s)*s
        log_p_sb = rowsum(log_p_sb_N, id)

        # random effect
        log_p_b_i = matrix(dnorm(b, rep(0, ncz), sd=sqrt(sigma2_b), log=TRUE),nrow=1)
        log_p_b = apply(log_p_b_i,2,function(x){rep(x,n)})
        
        # joint probability
        p_ysb = exp(log_p_yb+log_p_sb+log_p_b)
        p_ys = c(p_ysb %*% wGH)
        # which(p_ys==0)
        p_ys[which(p_ys==0)] = 1e-10
        p_bys = p_ysb/p_ys

        post_b = p_bys %*% (b * wGH) 
        post_b2 = p_bys %*% (b^2 * wGH) 
        post_vb = post_b2 - post_b^2
        post_b_N = rep(post_b,ni)
        post_b2_N = rep(post_b2,ni)

        res1 = y-eta_yx
        (sigma2_e_new = sum_rmna(sum_rmna(res1^2)-2*rho_new*sum_rmna(res1*post_b_N)+rho^2*sum_rmna(post_b2_N))/N)
        if(sigma2_e_new<=0){
            # annother approach: directly integrate b
            res1 = y-eta_yx
            nna_index = which(!is.na(res1))
            res_sq = ((res1-bN)*(res1-bN))[nna_index,]
            id_nna = id[nna_index]
            id_select = which(rownames(p_ysb) %in% unique(id_nna))
            p_ysb_nna = p_ysb[rownames(p_ysb) %in% unique(id_nna),]
            p_ys_nna = p_ys[id_select]
            (sigma2_e_new = sum((rowsum(res_sq,id_nna)*p_ysb_nna/p_ys_nna)%*%wGH)/dim(res_sq)[1])
        }

        (sigma2_b_new = sum_rmna(post_b^2+post_vb)/n)

        # update parameter values
        betas <- betasrho_new[1:ncx]
        rho <- betasrho_new[ncx+1]
        gammas <- gammas_new
        sigma2_e <- sigma2_e_new
        sigma2_b <- sigma2_b_new

        # compute log-likelihood
        log_p_ys <- log(p_ys)
        lgLik[it] <- sum(log_p_ys[is.finite(log_p_ys)], na.rm = TRUE)

        # print results if verbose
        if (control$verbose) {
            cat("\n\niter:", it, "\n")
            cat("log-likelihood:", lgLik[it], "\n")
            cat("betas:", round(betas, 4), "\n")
            cat("sigma2_e:", round(sigma2_e, 4), "\n")
            cat("gammas:", round(gammas, 4), "\n")
            cat("sigma2_b:", round(sigma2_b, 4), "\n")
            cat("rho:", round(rho, 4), "\n")
        }

        # check convergence
        if(it>1){
            if (lgLik[it] > lgLik[it - 1]) {
                thetas1 = c(Y_mat[it - 1, ], S_mat[it - 1, ], B_mat[it - 1, ], Rho_mat[it - 1])
                thetas2 = c(Y_mat[it, ], S_mat[it, ], B_mat[it, ], Rho_mat[it])
                check1 <- max(abs(thetas2 - thetas1) / (abs(thetas1) + control$tol1)) < control$tol2
                check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * (abs(lgLik[it - 1]) + control$tol3)
                if (check1 || check2) {
                # if (check1) {
                    conv <- FALSE # why set this to FALSE?
                    if (control$verbose)
                        conv = FALSE
                        cat("\n\nconverged!\n")
                    break
                }
            }
        }

    }

    # if there's no converge message, print "do not converge" and exit the function
    if (it == control$maxiter) {
        betas = NA
        stop("\n\nnot converged! Increase maxiter. \n")
    }

    betas
    return(betas)
}





