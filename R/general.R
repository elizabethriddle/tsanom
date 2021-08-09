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
