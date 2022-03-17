# Background functions

gaussian_unknown_mean_var_BFF_sequential <- function(data_new,N_prev,D_prev,M_prev,mu_0,sigma2_0,alpha_0,beta_0,alpha=39,beta=1.8,t=0){
  lambda_estimate <- optimise(lambda_posterior_unk_mu_sigma2_sequential_2_beta_prior_log,c(0,1),data_new=data_new,N_prev=N_prev,D_prev=D_prev,M_prev=M_prev,mu_0=mu_0,sigma2_0=sigma2_0,alpha_0=alpha_0,beta_0=beta_0,alpha=alpha,beta=beta,maximum=TRUE)
  N_new <- (lambda_estimate$maximum)*N_prev + (data_new)
  D_new <- (lambda_estimate$maximum)*D_prev + 1
  M_new <- (lambda_estimate$maximum)*M_prev + (data_new^2)

  mu_estimate <- (N_new+(mu_0/sigma2_0))/(D_new+(1/sigma2_0))
  sigma2_estimate <- (beta_0+0.5*((mu_0^2/sigma2_0)+M_new-((N_new+(mu_0/sigma2_0))^2/(D_new+(1/sigma2_0)))))/((0.5*D_new)+alpha_0+1)
  return(list(lambda_estimate=lambda_estimate$maximum,mu_estimate=mu_estimate,sigma2_estimate=sigma2_estimate,N_new=N_new,D_new=D_new,M_new=M_new))
}


lambda_posterior_unk_mu_sigma2_sequential_2_beta_prior_log <- function(lambda,data_new,N_prev,D_prev,M_prev,mu_0,sigma2_0,alpha_0,beta_0,alpha=39,beta=1.8){


  output_value <- (log(lambda)*{alpha-1})+(log(1-lambda)*{beta-1})+
    -0.5*log((lambda*D_prev)+1+(1/sigma2_0))+
    lgamma((0.5*lambda*D_prev)+0.5+alpha_0)+
    (log(beta_0+0.5*((mu_0^2/sigma2_0)+(lambda*M_prev)+(data_new^2)-((((lambda*N_prev)+data_new+(mu_0/sigma2_0))^2)/((lambda*D_prev)+1+(1/sigma2_0)))))*(-((0.5*lambda*D_prev)+0.5+(alpha_0))))+
    (log(1/(2*pi))*(0.5*lambda*D_prev+0.5))
  return(output_value)
}


lambda_posterior_unk_mu_sigma2_sequential_2_beta_prior_unnorm <- function(lambda,data_new,N_prev,D_prev,M_prev,mu_0,sigma2_0,alpha_0,beta_0,alpha=39,beta=1.8){


  if(length(lambda)>1){
    test_values <- seq(0.7,0.95,length.out = 20)
    rescaler <- max(sapply(test_values,lambda_posterior_unk_mu_sigma2_sequential_2_beta_prior_log,data_new=data_new,N_prev=N_prev,D_prev=D_prev,M_prev=M_prev,mu_0=mu_0,sigma2_0=sigma2_0,alpha_0=alpha_0,beta_0=beta_0,alpha=alpha,beta=beta))
    results <- sapply(lambda,lambda_posterior_unk_mu_sigma2_sequential_2_beta_prior_log,data_new=data_new,N_prev=N_prev,D_prev=D_prev,M_prev=M_prev,mu_0=mu_0,sigma2_0=sigma2_0,alpha_0=alpha_0,beta_0=beta_0,alpha=alpha,beta=beta)
    exp_results <- exp(results-rescaler)
    exp_results[is.infinite(exp_results)] <- 1e300
  }
  else{
    test_values <- seq(0.7,0.95,length.out = 20)
    rescaler <- max(sapply(test_values,lambda_posterior_unk_mu_sigma2_sequential_2_beta_prior_log,data_new=data_new,N_prev=N_prev,D_prev=D_prev,M_prev=M_prev,mu_0=mu_0,sigma2_0=sigma2_0,alpha_0=alpha_0,beta_0=beta_0,alpha=alpha,beta=beta))
    results <- lambda_posterior_unk_mu_sigma2_sequential_2_beta_prior_log(lambda=lambda,data_new=data_new,N_prev=N_prev,D_prev=D_prev,M_prev=M_prev,mu_0=mu_0,sigma2_0=sigma2_0,alpha_0=alpha_0,beta_0=beta_0,alpha=alpha,beta=beta)
    exp_results <- exp(results-rescaler)
    exp_results[is.infinite(exp_results)] <- 1e300
  }

  return(exp_results)
}


p_value_calculator_lambda <- function(FUNC,v,...){
  # Find normalising constant:
  normalising_constant <- integrate(FUNC,0,1,...)$value
  p_val<-integrate(FUNC, lower = 0, upper = v,...)$value/normalising_constant

  return(p_val)

}


BFF_predictive_p_value_calculator_approx <- function(x_new,muhat,N_value,D_value,M_value,mu_0_prior_value,sigma2_0_prior_value,sigma2_alpha_0_prior_value,sigma2_beta_0_prior_value){
  khat <- (1/sigma2_0_prior_value)+D_value
  alphahat <- (D_value/2) + sigma2_alpha_0_prior_value
  betahat <- sigma2_beta_0_prior_value + (0.5*(((mu_0_prior_value^2)/sigma2_0_prior_value)+M_value-((N_value+(mu_0_prior_value/sigma2_0_prior_value))^2/(D_value+(1/sigma2_0_prior_value)))))
  precision_val <- (alphahat*khat)/((1+khat)*betahat)
  df_val <- 2*alphahat

  if(precision_val<0){
    precision_val<-1e-9
  }

  find_maximum_post_pred <- muhat
  if(x_new<find_maximum_post_pred){
    lower_p_value <- pstp(x_new,mu=muhat,tau=precision_val,nu=df_val)

    calculated_p_value <- 2*lower_p_value
  }
  else{
    upper_p_value <- pstp(x_new,mu=muhat,tau=precision_val,nu=df_val,lower.tail = F)

    calculated_p_value <- 2*upper_p_value
  }
  return(calculated_p_value)
}


threshold_p_value_detector <- function(p_values,threshold_val,burn_out){
  anom_p_vals <- which(p_values<threshold_val)
  difference_vals <- diff(anom_p_vals)
  not_anom<-which(difference_vals>burn_out)
  true_anom_p_vals <- anom_p_vals[c(1,not_anom+1)]
  return(true_anom_p_vals)
}


threshold_p_value_detector_seq <- function(i,p_value,current_anom,threshold_val,burn_out){
  is_anom <- p_value<threshold_val
  if(length(current_anom)==0){
    out_grace <- T
  }else{
    out_grace <- (i-max(current_anom))>burn_out
  }
  if((is_anom) & out_grace){
    return(c(current_anom,i))
  }else{
    return(current_anom)
  }
}


poisson_unknown_rate_BFF_sequential <- function(data_new,N_prev,D_prev,F_prev,alpha_0,beta_0,alpha,beta){
  lambda_estimate <- optimise(lambda_posterior_poisson_unk_rate_sequential_beta_prior_log,c(0,1),data_new=data_new,N_prev=N_prev,D_prev=D_prev,F_prev=F_prev,alpha_0=alpha_0,beta_0=beta_0,alpha=alpha,beta=beta,maximum=TRUE)
  N_new <- (lambda_estimate$maximum)*N_prev + (data_new)
  D_new <- (lambda_estimate$maximum)*D_prev + 1
  F_new <- (lambda_estimate$maximum)*F_prev + (log(factorial(data_new)))
  gamma_estimate <- (alpha_0+N_new-1)/(D_new+beta_0)
  return(list(lambda_estimate=lambda_estimate$maximum,gamma_estimate=gamma_estimate,N_new=N_new,D_new=D_new,F_new=F_new))
}


lambda_posterior_poisson_unk_rate_sequential_beta_prior_log <- function(lambda,data_new,N_prev,D_prev,F_prev,alpha_0,beta_0,alpha,beta){
  output_value <- (log(lambda)*{alpha-1})+(log(1-lambda)*{beta-1})+
    lgamma(alpha_0+(lambda*N_prev)+data_new)-
    ((alpha_0+(lambda*N_prev)+data_new)*log(beta_0+(lambda*D_prev)+1))-
    (lambda*F_prev+lfactorial(data_new))
  return(output_value)
}



performance_large_changepoint <- function(estimated_cp,true_cp,datasize,D_period){
  false_positives <- setdiff(estimated_cp,true_cp+rep(c(0:(D_period)),each=length(true_cp)))
  true_positives <- setdiff(estimated_cp,false_positives)
  true_positives_refined <- c()
  AR1 <- c()
  for(i in c(1:length(true_cp))){
    possible_cp <- true_positives[(true_positives>=true_cp[i]) & (true_positives<=(true_cp[i]+D_period))]
    if(!is.null(possible_cp)){
      AR1 <- c(AR1,(possible_cp[which.min(possible_cp-true_cp[i])]-true_cp[i]))
      true_positives_refined <- c(true_positives_refined,possible_cp[which.min(possible_cp-true_cp[i])])
    }
  }
  AR1 <- mean(AR1,na.rm=T)
  FP <- length(false_positives)
  TotalPos <- (length(true_positives_refined)+length(false_positives))
  AR0 <- mean(diff(sort(false_positives)),na.rm=T)
  CCD <- length(true_positives_refined)/length(true_cp)
  DNF <- length(true_positives_refined)/(length(true_positives_refined)+length(false_positives))
  F1_val <- 2* (DNF*CCD)/(DNF+CCD)
  if(is.na(F1_val)){
    F1_val <- 0
  }
  return(list(F1=F1_val,ARL0=AR0,ARL1=AR1,CCD=CCD,DNF=DNF,FP=FP,Total_Positive=TotalPos))

}
