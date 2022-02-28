#' FDARIMA anomaly detector
#'
#' Detects anomalies by separately detecting anomalies using long-term and short-term
#' models and combining their p-values using Fishers product test statistic to give
#' final anomaly score that indicates both point and contextual anomalies.
#'
#' @param datamat a n*p matrix containing the time series where n is
#' the length of the period and p is the total number of periods observed
#' @param no_days_train number of periods of data to use as a burnin period
#' @param max_no_days this is maximum numbers of days that need to be forecasted using FDA model usually p-1
#' @param fda_forgetting_factor this is the intial forgetting factor to be used within exponential weighting (this is updated within the model)
#' @param arima_window_size number of historic data points to use to build ARIMA model
#' @param learning_rate the learning rate within stochastic gradient descent to update lambda
#' to be excluded from the model
#' @param threshold_val the threshold value used to detect anomalies
#'
#' @return returns list containing a matrix which contains the forecasts, and anomalies detected by each procedure and combined anomalies. Also
#' outputs a separate vector containing the combined anomalies only
#' @export
FDARIMA<-function(datamat,no_days_train,max_no_days,fda_forgetting_factor=0.99,arima_window_size=60,learning_rate=1e-3,threshold_val=0.005){
  datavec<-c(datamat)
  fda_predictions <- c()
  fda_variances <- c()
  arima_predictions <-c()
  strange_days <- c()
  strange_bins <- c()
  strange_data_points <- c()
  FDA_AE<-c()
  ARIMA_AE<-c()
  for(no_days in c(no_days_train:max_no_days)){
    print(no_days)
    if(no_days==no_days_train){
      # Fit poisson weighted FDA:
      lambda<-fda_forgetting_factor
      lambda_over_time<-c(rep(lambda,times=nrow(datamat)))
      no_obs<-no_days
      previous_weights <- lambda^{no_obs-c(1:(no_obs-1))}
      edge_data<-datamat[,c(1:no_obs)]
      fda_mean_weight<-fda_mean_weighted(edge_data,c(1:nrow(datamat)),102,4,c(previous_weights,1),includebpoints=FALSE)
      fda_predictions<-c(fda_predictions,fda_mean_weight$mean_curve)
      fda_variances<-c(fda_variances,fda_mean_weight$st_curve)
      FDA_AE<-c(FDA_AE,abs(fda_mean_weight$mean_curve-datamat[,(no_obs+1)]))

      #Update lambda:
      new_lambda<-update_fda_ff_adaptive(lambda,previous_weights,learning_rate,edge_data,datamat[,(no_obs+1)])
      if((0.6<=new_lambda) & (new_lambda<1)){
        lambda <- new_lambda
        lambda_over_time<-c(lambda_over_time,rep(new_lambda,times=nrow(datamat)))
      }else{
        lambda_over_time<-c(lambda_over_time,rep(lambda,times=nrow(datamat)))
      }

      # Now perform ARIMA and at each step check to see if anomaly:
      for(i in c(((no_days)*nrow(datamat)):((((no_days+1)*nrow(datamat))-1)))){
        fit<-auto.arima(datavec[setdiff(c((i-arima_window_size+1):i),strange_data_points)],max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
        fcast<-forecast(fit,h=1)$mean[1]
        arima_predictions<-c(arima_predictions,fcast)
        ARIMA_AE<-c(ARIMA_AE,abs(fcast-datavec[(i+1)]))
      }
    }
    else{
      # Fit poisson weighted FDA:
      no_obs<-no_days
      previous_weights <- lambda*c(previous_weights,1)
      edge_data<-datamat[,c(1:no_obs)]
      fda_mean_weight<-fda_mean_weighted_extreme(edge_data,c(1:nrow(datamat)),102,4,c(previous_weights,1),strange_days,strange_bins,includebpoints=FALSE)
      fda_predictions<-c(fda_predictions,fda_mean_weight$mean_curve)
      fda_variances<-c(fda_variances,fda_mean_weight$st_curve)
      FDA_AE<-c(FDA_AE,abs(fda_mean_weight$mean_curve-datamat[,(no_obs+1)]))

      if(no_days !=max_no_days){
        #Update lambda:
        new_lambda<-update_fda_ff_adaptive(lambda,previous_weights,learning_rate,edge_data,datamat[,(no_obs+1)])
        if((0.6<=new_lambda) & (new_lambda<1)){
          lambda <- new_lambda
          lambda_over_time<-c(lambda_over_time,rep(new_lambda,times=nrow(datamat)))
        }else{
          lambda_over_time<-c(lambda_over_time,rep(lambda,times=nrow(datamat)))
        }
      }
      # Now perform ARIMA and at each step check to see if anomaly:
      for(i in c(((no_days)*nrow(datamat)):((((no_days+1)*nrow(datamat))-1)))){
        arima_data_input <- datavec[setdiff(c((i-arima_window_size+1):i),strange_data_points)]
        if(length(arima_data_input)==0){
          while(length(arima_data_input)==0){
            arima_data_input<-datavec[c((i-1):i)]
          }
        }
        fit<-auto.arima(arima_data_input,max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
        fcast<-forecast(fit,h=1)$mean[1]
        arima_predictions<-c(arima_predictions,fcast)
        ARIMA_AE<-c(ARIMA_AE,abs(fcast-datavec[(i+1)]))
        # Now need to assess whether there are any anomalies
        FDA_p_value_point<-conformal_prediction_p_val(FDA_AE[1:(length(ARIMA_AE))])
        ARIMA_p_value_point<-conformal_prediction_p_val(ARIMA_AE)
        # Now combine using Fishers Product test statistic
        fisher_score <- -2*(log(FDA_p_value_point)+log(ARIMA_p_value_point))
        if(fisher_score>qchisq(1-threshold_val,4)){
          # Anomalous!
          strange_bins<-c(strange_bins,((i)%%nrow(datamat))+1)
          strange_data_points<-c(strange_data_points,i+1)
        }

      }

    }
  }
  anomalous_indicator <- rep(FALSE,length(ARIMA_AE))
  anomalous_indicator[strange_data_points-(no_days_train*nrow(datamat))] <- TRUE
  prediction_matrix<-data.frame(True = datavec[-c(1:(no_days_train*nrow(datamat)))],ARIMA = arima_predictions,ARIMA_AE = ARIMA_AE,FDA=fda_predictions,FDA_AE = FDA_AE,Anomalous=anomalous_indicator,Lambda_val = lambda_over_time)

  return(list(prediction_matrix=prediction_matrix,strange_data_points=strange_data_points))
}


#' FDARIMA anomaly detector initialiser
#'
#' Detects anomalies by separately detecting anomalies using long-term and short-term
#' models and combining their p-values using Fishers product test statistic to give
#' final anomaly score that indicates both point and contextual anomalies. This is the funciton to initiate
#' running in real time. Needed to initially start the model. Anomalies will not be detected until the conformal
#' prediction sliding window is filled (but model forecasts will still be computed).
#'
#' @param datamat a n*p matrix containing the time series where n is
#' the length of the period and p is the total number of periods observed
#' @param frequency_of_data this is the length of the period nature of the data
#' @param no_cycles_train number of periods of data to use as a burnin period
#' @param fda_forgetting_factor this is the intial forgetting factor to be used within exponential weighting (this is updated within the model)
#' @param arima_window_size number of historic data points to use to build ARIMA model
#' @param learning_rate the learning rate within stochastic gradient descent to update lambda
#' to be excluded from the model
#' @param threshold_val the threshold value used to detect anomalies
#' @param sliding_window_size the size of the sliding window for conformal prediction. This will determine whether can
#'
#' @return returns list containing a matrix which contains the forecasts, and anomalies detected by each procedure and combined anomalies. Also
#' outputs a separate vector containing the combined anomalies only
#' @export
FDARIMA_initial<-function(datamat,frequency_of_data,no_cycles_train,fda_forgetting_factor=0.99,arima_window_size=60,learning_rate=1e-3,threshold_val=0.005,sliding_window_size=14400){
  datavec<-c(datamat)
  fda_predictions <- c()
  arima_predictions <-c()
  strange_data_points <- c()
  FDA_AE<-c()
  ARIMA_AE<-c()
  max_no_days <- ncol(datamat)-1

  # Initialisation
  no_days<-no_cycles_train
  # Fit weighted FDA:
  lambda<-fda_forgetting_factor
  lambda_over_time<-c(rep(lambda,times=frequency_of_data))
  no_obs<-no_days
  previous_weights <- lambda^{no_obs-c(1:(no_obs-1))}
  sum_weights <- sum(previous_weights)+1
  edge_data<-datamat[,c(1:no_obs)]
  fda_mean_weight<-fda_mean_weighted(edge_data,c(1:frequency_of_data),102,4,c(previous_weights,1),includebpoints=FALSE)
  fda_predictions<-c(fda_predictions,fda_mean_weight$mean_curve)
  #FDA_AE<-c(FDA_AE,abs(fda_mean_weight$mean_curve-datamat[,(no_obs+1)]))



  # Now perform intial ARIMA using data in training period:
  i <- no_cycles_train*frequency_of_data
  fit <- auto.arima(datavec[c((i-arima_window_size+1):i)],max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
  fcast <- forecast(fit,h=1)$mean[1]
  arima_predictions <- c(arima_predictions,fcast)


  # Prepare for updates:
  weight_matrix <- matrix(rep(previous_weights/lambda,each=frequency_of_data),ncol=(no_cycles_train-1),nrow=frequency_of_data)
  weight_matrix <- weight_matrix * (1/rowSums(weight_matrix))
  previous_y_forecast <- rowSums((weight_matrix*datamat[,c(1:(no_cycles_train-1))]))

  current_y_forecast <- fda_mean_weight$weighted_average
  current_lambda <- lambda
  current_data <- c(datamat[,no_cycles_train])
  sum_weights_vec <- rep(sum(previous_weights/lambda),times=frequency_of_data)
  sum_weights_vec_new <- rep(sum(previous_weights)+1,times=frequency_of_data)
  new_data <-c()

  # RUN NOT ON INITIAL
  for(i in c((((no_cycles_train)*frequency_of_data)+1):((ncol(datamat)*frequency_of_data)))){
    # Assess if anomalous
    FDA_AE<-c(FDA_AE,abs(fda_predictions[i-(((no_cycles_train)*frequency_of_data))]-datavec[i]))
    ARIMA_AE<-c(ARIMA_AE,abs(fcast-datavec[(i)]))
    # Now need to assess whether there are any anomalies
    if(length(FDA_AE)>=sliding_window_size){
      FDA_p_value_point <- conformal_prediction_p_val(FDA_AE[((length(FDA_AE)-sliding_window_size+1):length(FDA_AE))])
      ARIMA_p_value_point<-conformal_prediction_p_val(ARIMA_AE[((length(ARIMA_AE)-sliding_window_size+1):length(ARIMA_AE))])
      # Now combine using Fishers Product test statistic
      fisher_score <- -2*(log(FDA_p_value_point)+log(ARIMA_p_value_point))
      if(fisher_score>qchisq(1-threshold_val,4)){
        strange_data_points<-c(strange_data_points,i)
      }
    }
    # Run ARIMA
    arima_data_input <- datavec[setdiff(c((i-arima_window_size+1):i),strange_data_points)]
    if(length(arima_data_input)==0){
      while(length(arima_data_input)==0){
        arima_data_input<-datavec[c((i-1):i)]
      }
    }
    fit<-auto.arima(arima_data_input,max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
    fcast<-forecast(fit,h=1)$mean[1]
    arima_predictions<-c(arima_predictions,fcast)

    new_data <- c(new_data,datavec[i])
    # If end of cycle update FDA:
    if(i%%frequency_of_data==0){
      # Put into vector which are anomalous
      if(length(strange_data_points)>0){
        which_anom_within_cycle <- (((strange_data_points[(i-frequency_of_data)<strange_data_points])-1)%%frequency_of_data)+1
        which_anomalous <- rep(0,times=frequency_of_data)
        which_anomalous[which_anom_within_cycle] <- 1
      }
      else{
        which_anomalous <- rep(0,times=frequency_of_data)
      }
      # Update lambda
      new_lambda <- update_fda_ff_adaptive_sequential(new_data,current_data,current_lambda,sum_weights_vec,current_y_forecast,previous_y_forecast,learning_rate,which_anomalous)
      if((0.6<=new_lambda) & (new_lambda<1)){
        current_lambda <- new_lambda
        lambda_over_time<-c(lambda_over_time,rep(new_lambda,times=frequency_of_data))
      }else{
        lambda_over_time<-c(lambda_over_time,rep(current_lambda,times=frequency_of_data))
      }
      # Run FDA
      sum_weights_vec <- sum_weights_vec_new
      fda_mean_weight<-fda_mean_weighted_extreme_update(new_data,current_lambda,sum_weights_vec,current_y_forecast,c(1:frequency_of_data),102,4,which_anomalous,includebpoints=FALSE,bpoints=NULL)
      fda_predictions<-c(fda_predictions,fda_mean_weight$mean_curve)
      previous_y_forecast <- current_y_forecast
      current_y_forecast <- fda_mean_weight$weighted_average
      sum_weights_vec_new <- fda_mean_weight$new_sum_weights_vec
      current_data <- new_data
      new_data <- c()
    }

  }



  anomalous_indicator <- rep(FALSE,length(ARIMA_AE))
  anomalous_indicator[strange_data_points-(no_cycles_train*frequency_of_data)] <- TRUE
  #print(c(length(datavec[-c(1:(no_cycles_train*frequency_of_data))]),length(arima_predictions[-c(length(arima_predictions))]),length(ARIMA_AE),length(fda_predictions),length(FDA_AE),length(anomalous_indicator),length(lambda_over_time)))
  #print(c(length(arima_predictions[-length(arima_predictions)]),length(ARIMA_AE[(length(ARIMA_AE)-length(arima_predictions)):(length(ARIMA_AE)-1)]),length(fda_predictions[(length(fda_predictions)-length(arima_predictions)):(length(fda_predictions)-1)])))
  #prediction_matrix<-data.frame(True = datavec[-c(1:(no_cycles_train*frequency_of_data))],ARIMA = arima_predictions[-length(arima_predictions)],ARIMA_AE = ARIMA_AE[(length(ARIMA_AE)-length(arima_predictions)):(length(ARIMA_AE)-1)],FDA=fda_predictions[(length(fda_predictions)-length(arima_predictions)):(length(fda_predictions)-1)],FDA_AE = FDA_AE[(length(FDA_AE)-length(arima_predictions)):(length(FDA_AE)-1)],Anomalous=anomalous_indicator[(length(anomalous_indicator)-length(arima_predictions)):(length(anomalous_indicator)-1)],Lambda_val = lambda_over_time[(length(lambda_over_time)-length(arima_predictions)):(length(lambda_over_time)-1)])
  prediction_matrix<-data.frame(True = datavec[-c(1:(no_cycles_train*frequency_of_data))],ARIMA = arima_predictions[-c(length(arima_predictions))],ARIMA_AE = ARIMA_AE,FDA=fda_predictions[1:length(ARIMA_AE)],FDA_AE = FDA_AE[1:length(ARIMA_AE)],Anomalous=anomalous_indicator,Lambda_val = lambda_over_time[1:length(ARIMA_AE)])

  if(length(FDA_AE)>sliding_window_size){
  return(list(prediction_matrix=prediction_matrix,strange_data_points=strange_data_points,
              new_data = new_data,
              current_data = current_data,
              fda_prediction_update = fda_predictions[c((length(fda_predictions)-frequency_of_data+1):length(fda_predictions))],
              ARIMA_input_data = datavec[c((length(datavec)-arima_window_size+1):length(datavec))],
              previous_FDA_wm = previous_y_forecast,
              current_FDA_wm = current_y_forecast,
              current_lambda = current_lambda,
              ARIMA_fcast = fcast,
              sum_weights_vec = sum_weights_vec,
              sum_weights_vec_new = sum_weights_vec_new,
              FDA_AE = FDA_AE[c(((length(FDA_AE)-sliding_window_size+1)):length(FDA_AE))],
              ARIMA_AE = ARIMA_AE[c(((length(ARIMA_AE)-sliding_window_size+1)):length(ARIMA_AE))]))
  }
  else{
  return(list(prediction_matrix=prediction_matrix,strange_data_points=strange_data_points,
              new_data = new_data,
              current_data = current_data,
              fda_prediction_update = fda_predictions[c((length(fda_predictions)-frequency_of_data+1):length(fda_predictions))],
              ARIMA_input_data = datavec[c((length(datavec)-arima_window_size+1):length(datavec))],
              previous_FDA_wm = previous_y_forecast,
              current_FDA_wm = current_y_forecast,
              current_lambda = current_lambda,
              ARIMA_fcast = fcast,
              sum_weights_vec = sum_weights_vec,
              sum_weights_vec_new = sum_weights_vec_new,
              FDA_AE = FDA_AE,
              ARIMA_AE = ARIMA_AE))
  }
}


#' FDARIMA anomaly detector updater
#'
#' Detects anomalies by separately detecting anomalies using long-term and short-term
#' models and combining their p-values using Fishers product test statistic to give
#' final anomaly score that indicates both point and contextual anomalies. This is the funciton to update
#' existing calculations as new data becomes available.
#'
#' @param i the current data point number
#' @param new_data_point the next value in the series that is used to update the models and asssess if anomalous
#' @param new_data this is a vector of all data points from the current cycle up till time i-1
#' @param current_data this is the value at the previous time point in the series
#' @param fda_predictions this is a vector detailing the current FDA forecasts for the current cycle
#' @param ARIMA_input_data this is a sliding window of data required to run the ARIMA model. It excludes new_data_point
#' @param previous_FDA_wm the weighted average from the previous update (at $i-2$)
#' @param current_FDA_wm the weighted average from update at time $i$
#' @param ARIMA_fcast this is the ARIMA forecast for time $i$
#' @param sum_weights_vec this is a vector detailing the weights for each time within the cycle(for the FDA model) for time i-2. This accounts for anomalous data.
#' @param sum_weights_vec_new this is a vector detailing the weights for each time within the cycle(for the FDA model) for time i-1. This accounts for anomalous data.
#' @param current_lambda this is the current estimate for the forgetting factor within FDA
#' @param FDA_AE the absolute error for the FDA mdoel. This is collected to run the conformal predictor and is kept at the lenght of the sliding window.
#' @param ARIMA_AE the absolute error for the ARIMA model. This is collected to run the conformal predictor and is kept at the lenght of the sliding window.
#' @param strange_data_points this outlines anomalous points and is a vector of times of anomalous points identified (only require  points in the past cycle).
#' @param frequency_of_data this is the length of the period nature of the data
#' @param fda_forgetting_factor this is the intial forgetting factor to be used within exponential weighting (this is updated within the model)
#' @param arima_window_size number of historic data points to use to build ARIMA model
#' @param learning_rate the learning rate within stochastic gradient descent to update lambda
#' to be excluded from the model
#' @param threshold_val the threshold value used to detect anomalies
#' @param sliding_window_size the size of the sliding window for conformal prediction. This will determine whether can
#'
#' @return returns list containing a matrix which contains the forecasts, and anomalies detected by each procedure and combined anomalies. Also
#' outputs a separate vector containing the combined anomalies only
#' @export
FDARIMA_updater<-function(i,new_data_point,new_data,current_data,fda_predictions,ARIMA_input_data,previous_FDA_wm,current_FDA_wm,ARIMA_fcast,sum_weights_vec,sum_weights_vec_new,current_lambda,FDA_AE,ARIMA_AE,strange_data_points,frequency_of_data,fda_forgetting_factor=0.99,arima_window_size=60,learning_rate=1e-3,threshold_val=0.005,sliding_window_size=14400){

  # Assess if anomalous
  FDA_AE<-c(FDA_AE,abs(fda_predictions[((i-1)%%frequency_of_data)+1]-new_data_point))
  FDA_fcast <- fda_predictions[((i-1)%%frequency_of_data)+1]
  ARIMA_AE<-c(ARIMA_AE,abs(ARIMA_fcast-new_data_point))

  # Now need to assess whether there are any anomalies
  if(length(FDA_AE)>sliding_window_size){
    FDA_p_value_point <- conformal_prediction_p_val(FDA_AE[((length(FDA_AE)-sliding_window_size+1):length(FDA_AE))])
    ARIMA_p_value_point<-conformal_prediction_p_val(ARIMA_AE[((length(ARIMA_AE)-sliding_window_size+1):length(ARIMA_AE))])

    # Now combine using Fishers Product test statistic
    fisher_score <- -2*(log(FDA_p_value_point)+log(ARIMA_p_value_point))
    if(fisher_score>qchisq(1-threshold_val,4)){
      strange_data_points<-c(strange_data_points,i)
    }
    FDA_AE <- FDA_AE[-c(1)]
    ARIMA_AE <- ARIMA_AE[-c(1)]
  }
  # Run ARIMA
  ARIMA_input_data_new <- c(ARIMA_input_data[-c(1)],new_data_point)
  arima_data_input <- ARIMA_input_data_new[setdiff(c((i-arima_window_size+1):i),strange_data_points)-(i-arima_window_size)]
  if(length(arima_data_input)==0){
    while(length(arima_data_input)==0){
      arima_data_input<-ARIMA_input_data_new[c(arima_window_size-1,arima_window_size)]
    }
  }
  fit<-auto.arima(arima_data_input,max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
  fcast<-forecast(fit,h=1)$mean[1]
  new_data <- c(new_data,new_data_point)
  # If end of cycle update FDA:
  if(i%%frequency_of_data==0){
    # Put into vector which are anomalous
    if(length(strange_data_points)>0){
      which_anom_within_cycle <- (((strange_data_points[(i-frequency_of_data)<strange_data_points])-1)%%frequency_of_data)+1
      which_anomalous <- rep(0,times=frequency_of_data)
      which_anomalous[which_anom_within_cycle] <- 1
    }
    else{
      which_anomalous <- rep(0,times=frequency_of_data)
    }
    # Update lambda
    new_lambda <- update_fda_ff_adaptive_sequential(new_data,current_data,current_lambda,sum_weights_vec,current_FDA_wm,previous_FDA_wm,learning_rate,which_anomalous)
    if((0.6<=new_lambda) & (new_lambda<1)){
      current_lambda <- new_lambda
    }
    # Run FDA
    sum_weights_vec <- sum_weights_vec_new
    fda_mean_weight<-fda_mean_weighted_extreme_update(new_data,current_lambda,sum_weights_vec,current_FDA_wm,c(1:frequency_of_data),102,4,which_anomalous,includebpoints=FALSE,bpoints=NULL)
    fda_predictions<-c(fda_predictions,fda_mean_weight$mean_curve)
    previous_FDA_wm <- current_FDA_wm
    current_FDA_wm <- fda_mean_weight$weighted_average
    sum_weights_vec_new <- fda_mean_weight$new_sum_weights_vec
    current_data <- new_data
    new_data <- c()
  }
  if(length(FDA_AE)>=sliding_window_size){
    return(list(new_data = new_data,
                current_data = current_data,
                strange_data_points = strange_data_points,
                fda_prediction_update = fda_predictions[c((length(fda_predictions)-frequency_of_data+1):length(fda_predictions))],
                ARIMA_input_data = ARIMA_input_data_new,
                previous_FDA_wm = previous_FDA_wm,
                current_FDA_wm = current_FDA_wm,
                current_lambda = current_lambda,
                ARIMA_fcast = fcast,
                FDA_fcast = FDA_fcast,
                sum_weights_vec = sum_weights_vec,
                sum_weights_vec_new = sum_weights_vec_new,
                FDA_AE = FDA_AE[c(((length(FDA_AE)-sliding_window_size+1)):length(FDA_AE))],
                ARIMA_AE = ARIMA_AE[c(((length(ARIMA_AE)-sliding_window_size+1)):length(ARIMA_AE))]))
  }
  else{
    return(list(new_data = new_data,
                current_data = current_data,
                strange_data_points = strange_data_points,
                fda_prediction_update = fda_predictions[c((length(fda_predictions)-frequency_of_data+1):length(fda_predictions))],
                ARIMA_input_data = ARIMA_input_data_new,
                previous_FDA_wm = previous_FDA_wm,
                current_FDA_wm = current_FDA_wm,
                current_lambda = current_lambda,
                ARIMA_fcast = fcast,
                FDA_fcast = FDA_fcast,
                sum_weights_vec = sum_weights_vec,
                sum_weights_vec_new = sum_weights_vec_new,
                FDA_AE = FDA_AE,
                ARIMA_AE = ARIMA_AE))
  }
}



#' Uses Shermann-Morrison theorem to update A matrix
#'
#' Updates A matrix for sequential updates of weighted linear regression
#'
#' @param A A matrix
#' @param v new values
#' @param lambda forgetting factor
#'
#' @return new A matrix
#' @export
newestAff<-function(A,v,lambda){
  (1/lambda)*A-(((1/(lambda^2))*A%*%v%*%t(v)%*%A)*as.numeric(solve(1+(1/lambda)*t(v)%*%A%*%v)))
}

#' Uses Shermann-Morrison theorem to update paramter vector
#'
#' Updates parameter vector matrix sequential updates of weighted linear regression
#'
#' @param beta parameter vector
#' @param A A matrix
#' @param v new values
#' @param y new response
#' @param lambda forgetting factor
#'
#' @return new parameter vector beta
#' @export
newestbetaff<-function(beta,A,v,y,lambda){
  beta+(((1/lambda)*A%*%v%*%(y-(t(v)%*%beta)))%*%as.numeric(solve(1+(1/lambda)*t(v)%*%A%*%v)))
}

#' Forecast combiner using weighted linear regression
#'
#' Combines short-term and long term forecasts using fixed exponentially weighted linear regression.
#'
#' @param True_values the real time series
#' @param FDA_values the forecasted FDA values for the input series
#' @param ARIMA_values the forecasted ARIMA values for the input series
#' @param start_val the number of points to initialise with
#' @param end_val when to stop the estimation
#' @param lambda the fixed forgetting factor to use
#'
#' @return combined forecast values
#' @export
regression_combination_forecast<-function(True_values,FDA_values,ARIMA_values,start_val,end_val,lambda=0.98){
  prediction_matrix <- data.frame(True = True_values,ARIMA=ARIMA_values,FDA=FDA_values)
  ARIMA_coef_forgetting<-c()
  FDA_coef_forgetting<-c()
  forgetting_predicted_value<-c()

  i<-start_val

  w_vec<-c()
  no_obs<-i
  for(j in 1:no_obs){
    w_vec<-c(w_vec,((lambda)^(no_obs-j)))
  }
  prediction_lm_forgetting<-lm(True~ARIMA+FDA+0,data=prediction_matrix[(1:i),],weights=w_vec)
  ARIMA_coef_forgetting<-c(ARIMA_coef_forgetting,prediction_lm_forgetting$coefficients[1])
  FDA_coef_forgetting<-c(FDA_coef_forgetting,prediction_lm_forgetting$coefficients[2])
  forgetting_predicted_value<-c(forgetting_predicted_value,predict(prediction_lm_forgetting,newdata=prediction_matrix[(i+1),]))

  A_train<-solve(t(as.matrix(prediction_matrix[(1:i),c(2,3)]))%*%as.matrix(prediction_matrix[(1:i),c(2,3)]),tol = 1e-19)
  A_new<-A_train

  beta_seq<-data.frame(ARIMA=prediction_lm_forgetting$coefficients[1],FDA=prediction_lm_forgetting$coefficients[2])
  beta_seq<-t(beta_seq)

  for(i in c((start_val+1):end_val)){

    v<-prediction_matrix[i,c(2,3)]
    y<-prediction_matrix[i,1]

    beta_seq<-newestbetaff(beta_seq,A_new,t(v),y,lambda)
    A_new<-newestAff(A_new,t(v),lambda)

    ARIMA_coef_forgetting<-c(ARIMA_coef_forgetting,beta_seq[1])
    FDA_coef_forgetting<-c(FDA_coef_forgetting,beta_seq[2])
    v_next <- prediction_matrix[(i+1),c(2,3)]
    forgetting_predicted_value<-c(forgetting_predicted_value,sum(v_next*beta_seq))
  }
  return(list(ARIMA_coef_forgetting=ARIMA_coef_forgetting,FDA_coef_forgetting=FDA_coef_forgetting,forgetting_predicted_value=forgetting_predicted_value))
}
