#' FDA model calculator for data without anomalies
#'
#' This function produces a weighted FDA curve for the input series.
#' The coefficients of the b-spline are calculated using poisson regression.
#' This gives the average curve, sd and a prediction interval for the curve.
#'
#' @param datamat a n*p matrix containing all previous periods of the series where n is
#' the length of the period and p is the total number of periods observed
#' @param rangex vector c(1:n), the point numbers in the period
#' @param nb number of basis functions to use
#' @param no order of basis functions
#' @param obsweights vector of weights for each column of the data
#' @param includebpoints A booliean. If TRUE then must specify breakpoints otherwise do not
#' @param bpoints Optional vector (required if includebpoints is TRUE). This details the breakpoints
#' to use for the basis functions
#'
#' @return List containing average, sd and upper and lower ranges.
#' @export
fda_mean_weighted<-function(datamat,rangex,nb,no,obsweights,includebpoints=FALSE,bpoints=NULL){
  n<-ncol(datamat)
  weight_matrix <- matrix(rep(obsweights,each=nrow(datamat)),nrow=nrow(datamat),ncol=ncol(datamat))
  weight_matrix <- weight_matrix * (1/rowSums(weight_matrix))
  data_mean <- rowSums((weight_matrix*datamat))
  if(includebpoints!=FALSE){
    basis_function <- create.bspline.basis(rangeval=c(min(rangex),max(rangex)),nbasis=nb,norder=no,breaks=bpoints)
  }
  else{
    basis_function <- create.bspline.basis(rangeval=c(min(rangex),max(rangex)),nbasis=nb,norder=no)
  }
  basis_results<-smooth.basis(rangex,data_mean,basis_function)
  resultfd<-basis_results$fd

  # Now find unbiased variance. First estimate Simga_e
  smoothdatamat<-eval.fd(rangex,resultfd)
  smoothdatamat<- smoothdatamat[,rep(1,times=n)]
  datares <- datamat - smoothdatamat
  w_sum<-sum(obsweights)
  w2_sum<-sum(obsweights^2)
  datavar<-apply((datares^2)*obsweights,1,sum)/(w_sum-(w2_sum/w_sum))

  # Now smooth the log of this and exponentiate result:
  var.fit<-smooth.basis(rangex,log(datavar),basis_function)
  var.fd<-var.fit$fd
  varvec<-exp(eval.fd(rangex,var.fd))
  SigmaE<-diag(as.vector(varvec))

  # Now calculate the actual variance of the estimated curve:
  y2cMap<-basis_results$y2cMap
  c2rMap <-eval.basis(rangex,basis_function)
  Sigmayhat <- c2rMap %*% y2cMap %*% SigmaE %*% t(y2cMap) %*% t(c2rMap)
  edge.stderr<-sqrt(diag(Sigmayhat))

  avmat<-eval.fd(rangex,resultfd)

  # Now return the average values, and the standard deviation values:
  return(list(mean_curve=avmat,st_curve=edge.stderr,weighted_average = data_mean))
}






#' FDA model calculator for data with anomalies
#'
#' This function produces a weighted FDA curve for the input series but excludes
#' anomalous points from its calculation.
#' The coefficients of the b-spline are calculated using poisson regression.
#' This gives the average curve, sd and a prediction interval for the curve.
#'
#' @param datamat a n*p matrix containing all previous periods of the series where n is
#' the length of the period and p is the total number of periods observed
#' @param rangex vector c(1:n), the point numbers in the period
#' @param nb number of basis functions to use
#' @param no order of basis functions
#' @param obsweights vector of weights for each column of the data
#' @param extreme_days a vector of values corresponds to the column of anomalous points
#' to be excluded from the model
#' @param extreme_bins a vector of values corresponding to the rows of anomalous points
#' to be excluded from the model
#' @param includebpoints A booliean. If TRUE then must specify breakpoints otherwise do not
#' @param bpoints Optional vector (required if includebpoints is TRUE). This details the breakpoints
#' to use for the basis functions
#'
#' @return List containing average, sd and upper and lower ranges.
#' @export
fda_mean_weighted_extreme<-function(datamat,rangex,nb,no,obsweights,extreme_days,extreme_bins,includebpoints=FALSE,bpoints=NULL){
  n<-ncol(datamat)
  weight_matrix <- matrix(rep(obsweights,each=nrow(datamat)),nrow=nrow(datamat),ncol=ncol(datamat))
  weight_matrix[extreme_bins,extreme_days] <- 0
  weight_matrix <- weight_matrix * (1/rowSums(weight_matrix))
  data_mean <- rowSums((weight_matrix*datamat))
  if(includebpoints!=FALSE){
    basis_function <- create.bspline.basis(rangeval=c(min(rangex),max(rangex)),nbasis=nb,norder=no,breaks=bpoints)
  }
  else{
    basis_function <- create.bspline.basis(rangeval=c(min(rangex),max(rangex)),nbasis=nb,norder=no)
  }
  basis_mat<-eval.basis(rangex,basis_function)
  basis_data<-as.data.frame(basis_mat)
  basis_data_mat<-cbind(data_mean,basis_data)
  basis_coeffs<-lm(data_mean~ . -1,data=basis_data_mat)

  avmat<-basis_coeffs$fitted.values

  # Will estimate the sd for the data using the general sd about this prediction! Nothing too fancy:
  # Non Weighted SD:
  # smoothdatamat<- as.matrix(basis_coeffs$fitted.values)[,rep(1,times=n)]
  # datares <- (datamat - smoothdatamat)^2
  # datares <- rowSums(datares)/n
  # edge.stderr<- sqrt(datares)

  # Weighted SD:
  # smoothdatamat<- as.matrix(basis_coeffs$fitted.values)[,rep(1,times=n)]
  # datares<-datamat - smoothdatamat
  # datares[weight_matrix==0]<-0
  # w_sum<-sum(obsweights)
  # w2_sum<-sum(obsweights^2)
  # datavar<-apply((datares^2)*obsweights,1,sum)/(w_sum-(w2_sum/w_sum))
  # edge.stderr<-sqrt(datavar)


  #basis_results<-smooth.basis(rangex,data_mean,basis_function)

  # resultfd<-basis_results$fd
  #
  # # Now find unbiased variance. First estimate Simga_e
  # smoothdatamat<-eval.fd(rangex,resultfd)
  # smoothdatamat<- smoothdatamat[,rep(1,times=n)]
  # datares <- datamat - smoothdatamat
  # w_sum<-sum(obsweights)
  # w2_sum<-sum(obsweights^2)
  # datavar<-apply((datares^2)*w_vec,1,sum)/(w_sum-(w2_sum/w_sum))
  #
  # # Now smooth the log of this and exponentiate result:
  # var.fit<-smooth.basis(rangex,log(datavar),basis_function)
  # var.fd<-var.fit$fd
  # varvec<-exp(eval.fd(rangex,var.fd))
  # SigmaE<-diag(as.vector(varvec))
  #
  # # Now calculate the actual variance of the estimated curve:
  # y2cMap<-basis_results$y2cMap
  # c2rMap <-eval.basis(rangex,basis_function)
  # Sigmayhat <- c2rMap %*% y2cMap %*% SigmaE %*% t(y2cMap) %*% t(c2rMap)
  # edge.stderr<-sqrt(diag(Sigmayhat))

  #avmat<-eval.fd(rangex,resultfd)

  # Now return the average values, and the standard deviation values:
  return(list(mean_curve=avmat))
}


#' FDA model calculator for data with anomalies for a single update from previous (sequential).
#'
#' This function produces a weighted FDA curve for the input series but excludes
#' anomalous points from its calculation.
#' The coefficients of the b-spline are calculated using linear regression.
#' This gives the average curve and the new sum of weights for each time within a cycle.
#'
#' @param current_data a vector detailing the data from the latest fully observed cycle
#' @param current_lambda the latest lambda estimate
#' @param sum_weights_vec the weight for each time within the cycle (these weights are
#' different for each point as anomalous points are excluded hence have weight 0)
#' @param forecast_values these are the corresponding forecast values for input current data
#' @param rangex vector c(1:n), the point numbers in the period
#' @param nb number of basis functions to use
#' @param no order of basis functions
#' @param which_anomalous details which of the current_data points are considered anomalous (boolean
#' vector with 1 if anomalous 0 if not)
#' @param includebpoints A booliean. If TRUE then must specify breakpoints otherwise do not
#' @param bpoints Optional vector (required if includebpoints is TRUE). This details the breakpoints
#' to use for the basis functions
#'
#' @return List containing average fda curve, the weighted average and the new updated sum_weights_vec which incorporates the
#' anomalous values
#' @export
fda_mean_weighted_extreme_update<-function(current_data,current_lambda,sum_weights_vec,forecast_values,rangex,nb,no,which_anomalous,includebpoints=FALSE,bpoints=NULL){
  should_include <- 1-which_anomalous
  # Calculate weighted average over historical data:
  denominator_wa <- (current_lambda*sum_weights_vec)+should_include
  numerator_wa <- (current_lambda*forecast_values*sum_weights_vec)+(should_include*current_data)
  n<-length(current_data)
  data_mean <- numerator_wa/denominator_wa

  if(includebpoints!=FALSE){
    basis_function <- create.bspline.basis(rangeval=c(min(rangex),max(rangex)),nbasis=nb,norder=no,breaks=bpoints)
  }
  else{
    basis_function <- create.bspline.basis(rangeval=c(min(rangex),max(rangex)),nbasis=nb,norder=no)
  }
  basis_mat<-eval.basis(rangex,basis_function)
  basis_data<-as.data.frame(basis_mat)

  basis_data_mat<-cbind(data_mean,basis_data)
  basis_coeffs<-lm(data_mean~ . -1,data=basis_data_mat)

  avmat<-basis_coeffs$fitted.values
  new_sum_weights_vec <- (sum_weights_vec*current_lambda) + should_include
  return(list(mean_curve = avmat,weighted_average = data_mean, new_sum_weights_vec = new_sum_weights_vec))
}



#' Stochastic gradient descent updater for forgetting factor
#'
#' This function updates the forgetting factor using stochastic
#' gradient descent.
#'
#' @param current_lambda the current value for lambda
#' @param alpha the learning rate which controls the change at each update (recommend 10^-3)
#' @param datamat a n*p matrix containing all previous periods of the series n is
#' the length of the period and p is the total number of periods observed
#' @param data_new a vector containing the values of the newly observed period
#'
#' @return new forgetting factor
#' @export
update_fda_ff <- function(current_lambda, alpha, datamat,data_new){
  # datamat bins x days
  n_val <- ncol(datamat)
  i_vals <- c(1:(n_val-1))
  weight_vals<-sapply(i_vals,FUN=ff_weight,lambda=current_lambda,n=n_val)
  n<-ncol(datamat)
  weight_matrix <- matrix(rep(weight_vals,each=nrow(datamat)),nrow=nrow(datamat),ncol=ncol(datamat))
  data_mult <- rowSums((weight_matrix*datamat))
  # Calculate the predicted x values:
  w_vec<-c()
  for(i in 1:n){
    w_vec<-c(w_vec,((current_lambda)^(n-i)))
  }
  weight_matrix <- matrix(rep(w_vec,each=nrow(datamat)),nrow=nrow(datamat),ncol=ncol(datamat))
  weight_matrix <- weight_matrix * (1/rowSums(weight_matrix))
  data_mean <- rowSums((weight_matrix*datamat))
  mult_update <- sum(sign(data_new-data_mean)*data_mult)
  return(current_lambda+((alpha/length(data_new))*mult_update))
}


#' Stochastic gradient descent updater for forgetting factor using adaptive approach
#'
#' This function updates the forgetting factor using stochastic
#' gradient descent.
#'
#' @param current_lambda the current value for lambda
#' @param previous_weights the previous weights for the data
#' @param alpha the learning rate which controls the change at each update (recommend 10^-3)
#' @param datamat a n*p matrix containing all previous periods of the series n is
#' the length of the period and p is the total number of periods observed
#' @param data_new a vector containing the values of the newly observed period
#'
#' @return new forgetting factor
#' @export
update_fda_ff_adaptive <- function(current_lambda,previous_weights, alpha, datamat,data_new){
  # datamat bins x days
  # previous_weights contains weights for days 1 to n_val - 1 as most recent has weight 1 that is independant of lambda except for denom! (already multiplied by current lambda)
  n_val <- ncol(datamat)
  prev_weights_without_lambda <- (previous_weights/current_lambda)[1:(n_val-1)]
  y_hat_prev_unweight <- rowSums(matrix(rep(prev_weights_without_lambda,each=nrow(datamat)),nrow=nrow(datamat),ncol=(ncol(datamat)-1))*datamat[,(1:n_val-1)])
  rhs_sums_numerator <- y_hat_prev_unweight - (datamat[,(n_val-1)]*sum(prev_weights_without_lambda))
  # weight_matrix_prev <- matrix(rep(previous_weights,each=nrow(datamat)),nrow=nrow(datamat),ncol=(ncol(datamat)-1))
  # rhs_sums_numerator <- rowSums(previous_weights*datamat[,(1:n_val-1)])-(datamat[,n_val]*sum(previous_weights))
  rhs_sums_denominator <- (1+sum(previous_weights))^2
  rhs_sums <- rhs_sums_numerator / rhs_sums_denominator


  weight_matrix <- matrix(rep(c(previous_weights,1),each=nrow(datamat)),nrow=nrow(datamat),ncol=(ncol(datamat)))
  data_mean <- (rowSums((weight_matrix*datamat[,(1:n_val)])))/(sum(previous_weights)+1)
  mult_update <- sum(sign(data_new-data_mean)*rhs_sums)
  return(current_lambda+((alpha/length(data_new))*mult_update))
}



#' Stochastic gradient descent updater for forgetting factor using adaptive approach
#'
#' This function updates the forgetting factor using stochastic
#' gradient descent.
#'
#' @param new_data a vector detailing the data from the latest fully observed cycle
#' @param current_data a vector detailing the data from the previous fully observed cycle
#' @param current_lambda the latest lambda estimate
#' @param sum_weights_vec the weight for each time within the cycle (these weights are
#' different for each point as anomalous points are excluded hence have weight 0). This is from the
#' previous cycle at this time.
#' @param forecast_values_new these are the corresponding forecast values for new_data
#' @param forecast_values_old these are the corresponding forecast values for current_data
#' @param alpha the learning rate which controls the change at each update (recommend 10^-3)
#' @param which_anomalous details which of the current_data points are considered anomalous (boolean
#' vector with 1 if anomalous 0 if not)
#'
#' @return new forgetting factor
#' @export
update_fda_ff_adaptive_sequential <- function(new_data,current_data,current_lambda,sum_weights_vec,forecast_values_new,forecast_values_old,alpha,which_anomalous){


  # datamat bins x days
  # previous_weights contains weights for days 1 to n_val - 1 as most recent has weight 1 that is independant of lambda except for denom! (already multiplied by current lambda)
  should_include <- 1-which_anomalous
  rhs_sums_numerator <- (forecast_values_old*sum_weights_vec)-(current_data*should_include*sum_weights_vec)
  rhs_sums_denominator <- (should_include+(current_lambda*sum_weights_vec))^2
  rhs_sums <- rhs_sums_numerator / rhs_sums_denominator

  mult_update <- sum(sign(new_data-forecast_values_new)*rhs_sums)
  return(current_lambda+((alpha/length(new_data))*mult_update))
}

#' Calculate Conformal P-Values
#'
#' This function takes a sequence of non conformity measures (in this case,
#' the absolute forecasting error is used) and produces a p-value for the
#' final value in the sequence.
#'
#' @param non_conform_measures sequence of non conformity measures
#' @return a p-value (between 0 and 1) for the final value in the input sequence
#' @export
conformal_prediction_p_val <- function(non_conform_measures){
  # Calculate p-value for last entry in that vector:
  no_greater_eq<-sum(non_conform_measures>=non_conform_measures[length(non_conform_measures)])
  p_val <- no_greater_eq/length(non_conform_measures)
  p_val
}

#' Simulate periodic data with point and contextual anomalies
#'
#' This function produces a periodic time series containing the specified number of
#' point and contextual anomalies. The data uses a sine curve for the long-term
#' component and an ARMA(1,1) process for the short-term structure. Gaussian noise
#' is added to the series
#'
#' @param n the length of the simulated sequence
#' @param burnin_period length of initial burn-in without anomalies
#' @param no_additive_long number of additive contextual anomalies
#' @param no_subtractive_long number of subtractive contextual anomalies
#' @param no_point number of point anomalies
#' @param day_freq the period length of the periodic cycle which is automatically set to 1440
#' @param ar_coef the AR parameter coefficient default = 0.05
#' @param ma_coef the MA parameter coefficient default = 0.05
#'
#' @return a list containing the series with and without anomalies added, and locations of the point and contextual anomalies
#'
#' @examples data_simulation_function(1440,100,2,2,30)
#'
#' @export
data_simulation_function <- function(n,burnin_period,no_additive_long,no_subtractive_long,no_point,day_freq=1440,ar_coef=0.05,ma_coef=0.05){
  period_investigate_long<-c(1:(day_freq*n))
  period_investigate<-c(1:(day_freq*n))
  LT_data_long<-100*(sin((((period_investigate-1)/day_freq)*(2*pi))-(pi/2))+1)

  #ARIMA_simulation_long <- arima.sim(list(order = c(1,1,1),ar=c(ar_coef),ma=c(ma_coef)), n = (day_freq*n))
  #combined_simulation_noise_long<-LT_data_long+ARIMA_simulation_long[-1]+rnorm((day_freq*n),0,sqrt(5))
  # seed 50, arima 2,1,2 0.01,0.01,0.1,0.1

  ARIMA_simulation_long <- arima.sim(list(order = c(1,0,1),ar=c(ar_coef),ma=c(ma_coef)), n = (day_freq*n))
  combined_simulation_noise_long<-LT_data_long+ARIMA_simulation_long+rnorm((day_freq*n),0,sqrt(5))



  anomalous_combined_simulation_noise_mix_long<-combined_simulation_noise_long
  #print("created non anom series")

  # Additive Long
  additive_long_x <- list()
  additive_long_y <- list()
  additive_long_y_noise <- list()
  possible_anomalous_points <- c(((burnin_period*day_freq)+1):(n*day_freq))
  for(i in c(1:no_additive_long)){
    x_anomalous_length <- sample(c(2,2.5,3,3.5,4,4.5,5),1)
    not_valid <- TRUE
    while(not_valid){
      proposed_period <- c(1:(60*x_anomalous_length))+sample(possible_anomalous_points,1)
      if(all(proposed_period %in% possible_anomalous_points)){
        not_valid <- FALSE
      }
    }
    additive_long_x[[i]] <- proposed_period
    additive_long_y[[i]] <- runif(1,50,100)*(sin((((c(1:(60*x_anomalous_length))-1)/(60*x_anomalous_length))*2*pi)-(pi/2))+1)
    additive_long_y_noise[[i]] <- additive_long_y[[i]]+rnorm(60*x_anomalous_length,0,sqrt(5))
    possible_anomalous_points <- setdiff(possible_anomalous_points,proposed_period)
    anomalous_combined_simulation_noise_mix_long[proposed_period] <- anomalous_combined_simulation_noise_mix_long[proposed_period] + additive_long_y_noise[[i]]
  }

  #print("additive long")
  # Subtractive Long
  subtractive_long_x <- list()
  subtractive_long_y <- list()
  subtractive_long_y_noise <- list()
  for(i in c(1:no_subtractive_long)){
    x_anomalous_length <- sample(c(2,2.5,3,3.5,4,4.5,5),1)
    not_valid <- TRUE
    while(not_valid){
      proposed_period <- c(1:(60*x_anomalous_length))+sample(possible_anomalous_points,1)
      if(all(proposed_period %in% possible_anomalous_points)){
        not_valid <- FALSE
      }
    }
    subtractive_long_x[[i]] <- proposed_period
    subtractive_long_y[[i]] <- -runif(1,50,100)*(sin((((c(1:(60*x_anomalous_length))-1)/(60*x_anomalous_length))*2*pi)-(pi/2))+1)
    subtractive_long_y_noise[[i]] <- subtractive_long_y[[i]]+rnorm(60*x_anomalous_length,0,sqrt(5))
    possible_anomalous_points <- setdiff(possible_anomalous_points,proposed_period)
    anomalous_combined_simulation_noise_mix_long[proposed_period] <- anomalous_combined_simulation_noise_mix_long[proposed_period] + subtractive_long_y_noise[[i]]
  }
  #print("sub long")


  # Additive Short
  x_anomalous_point <- sample(possible_anomalous_points,no_point)
  y_anomalous_point <- sample(c(1,-1),no_point,replace = T)*rnorm(no_point,200,20)
  anomalous_combined_simulation_noise_mix_long[x_anomalous_point] <- anomalous_combined_simulation_noise_mix_long[x_anomalous_point] + y_anomalous_point
  #print("point")

  anomalous_combined_simulation_noise_mix_long_adjust_value<-0
  if(min(anomalous_combined_simulation_noise_mix_long)<0){
    anomalous_combined_simulation_noise_mix_long_adjust_value<-min(anomalous_combined_simulation_noise_mix_long)
    anomalous_combined_simulation_noise_mix_long<-anomalous_combined_simulation_noise_mix_long-min(anomalous_combined_simulation_noise_mix_long)
    combined_simulation_noise_long<-combined_simulation_noise_long-anomalous_combined_simulation_noise_mix_long_adjust_value
  }


  x_anomalous<-unique(x_anomalous_point,unlist(additive_long_x),unlist(subtractive_long_x))
  x_anomalous_point<-unique(x_anomalous_point)
  x_anomalous_context<-unique(c(unlist(additive_long_x),unlist(subtractive_long_x)))

  return(list(y_no_anomalies=combined_simulation_noise_long,
              y_anomalies=anomalous_combined_simulation_noise_mix_long,
              x_anomalous=x_anomalous,
              x_anomalous_point=x_anomalous_point,
              x_anomalous_context=x_anomalous_context))
}


#' Anomaly detection performance calculator
#'
#' This function measures the performance of the anomaly detection procedure in terms of
#' F1, Recall and precision. Additionally it splits this up into correctly identified
#' point and contextual anomalies.
#'
#' @param estimated_anomalies the anomalies detected by procedure
#' @param true_anomalies_point the true point anomalies locations
#' @param true_anomalies_contextual the true contextual anomaly locations
#'
#' @return a list containing the performance of detection procedure overall and split into point and contextual anomalies
#'
#' @examples anomaly_detection_performance_split(c(1,5,4),c(1),c(5))
#'
#' @export
anomaly_detection_performance_split <- function(estimated_anomalies,true_anomalies_point,true_anomalies_contextual){

  # Point
  detect_anom_without_contextual <- setdiff(estimated_anomalies,true_anomalies_contextual)
  false_positives <- setdiff(detect_anom_without_contextual,true_anomalies_point)
  true_positives <- setdiff(detect_anom_without_contextual ,false_positives)

  CCD_point <- length(true_positives)/length(true_anomalies_point)
  DNF_point <- length(true_positives)/(length(true_positives)+length(false_positives))
  F1_val_point <- 2* (DNF_point*CCD_point)/(DNF_point+CCD_point)
  if(is.na(F1_val_point)){
    F1_val_point <- 0
  }


  # Contextual
  detect_anom_without_point <- setdiff(estimated_anomalies,true_anomalies_point)
  false_positives <- setdiff(detect_anom_without_point,true_anomalies_contextual)
  true_positives <- setdiff(detect_anom_without_point ,false_positives)

  CCD_contextual <- length(true_positives)/length(true_anomalies_contextual)
  DNF_contextual <- length(true_positives)/(length(true_positives)+length(false_positives))
  F1_val_contextual <- 2* (DNF_contextual*CCD_contextual)/(DNF_contextual+CCD_contextual)
  if(is.na(F1_val_contextual)){
    F1_val_contextual <- 0
  }


  true_anomalies <- unique(c(true_anomalies_point,true_anomalies_contextual))
  false_positives <- setdiff(estimated_anomalies,true_anomalies)
  true_positives <- setdiff(estimated_anomalies,false_positives)
  CCD <- length(true_positives)/length(true_anomalies)
  DNF <- length(true_positives)/(length(true_positives)+length(false_positives))
  F1_val <- 2* (DNF*CCD)/(DNF+CCD)
  if(is.na(F1_val)){
    F1_val <- 0
  }


  return(list(F1_val_point=F1_val_point,Recall_point=CCD_point,Precision_point=DNF_point,
              F1_val_contextual=F1_val_contextual,Recall_contextual=CCD_contextual,Precision_contextual=DNF_contextual,
              F1_val=F1_val,Recall=CCD,Precision=DNF))

}
