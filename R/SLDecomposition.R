#' SL Decomposition anomaly detector
#'
#' Detects anomalies by separately detecting anomalies using long-term and combined forecast
#' from decomposition and combining their p-values using Fishers product test statistic to give
#' final anomaly score that indicates both point and contextual anomalies.
#'
#' @param datamat a n*p matrix containing the time series where n is
#' the length of the period and p is the total number of periods observed
#' @param arima_window_size number of historic data points to use to build ARIMA model
#' @param no_days_train number of periods of data to use as a burnin period
#' @param starting_no_obs this is the number of days to start collecting residuals in long-term process to
#' forecast the short-term process using ARIMA
#' @param fda_forgetting_factor this is the intial forgetting factor to be used within exponential weighting (this is updated within the model)
#' @param learning_rate the learning rate within stochastic gradient descent to update lambda
#' to be excluded from the model
#' @param threshold_val the threshold value used to detect anomalies
#'
#' @return returns list containing a matrix which contains the forecasts, and anomalies detected by each procedure and combined anomalies. Also
#' outputs a separate vector containing the combined anomalies only
#' @export
SLDecomposition <- function(datamat,frequency_of_data=1440,arima_window =60,no_days_train=10,starting_no_obs=5,fda_forgetting_factor=0.99,learning_rate=1e-3,threshold_val = 0.005){
  datavec <- c(datamat)
  SLD_fda_anomaly<-c()
  SLD_fda_error_anomaly<-c()
  SLD_arima_anomaly<-c()
  # now sig level set to 0.01

  lambda<-fda_forgetting_factor
  lambda_over_time<-c(rep(lambda,times=frequency_of_data))
  no_obs<-starting_no_obs
  previous_weights <- lambda^{no_obs-c(1:(no_obs-1))}
  edge_data<-datamat[,c(1:no_obs)]
  fda_mean_weight<-fda_mean_weighted(edge_data,c(1:nrow(datamat)),102,4,c(previous_weights,1),includebpoints=FALSE)
  SLD_fda_error_anomaly<-c(SLD_fda_error_anomaly,datamat[,(no_obs+1)]-fda_mean_weight$mean_curve)
  SLD_arima_anomaly<-c()
  SLD_combined_prediction_anomaly<-c()
  SLD_combined_prediction_error_anomaly<-c()
  strange_bins_sld_anomaly<-c()
  strange_days_sld_anomaly<-c()
  strange_bins_sld_FDA_anomaly<-c()
  strange_days_sld_FDA_anomaly<-c()
  strange_points_sld_anomaly<-c()

  strange_points_sld_FDA_anomaly<-c()
  strange_points_sld_combined_anomaly<-c()


  # Now Update Lambda:
  #Update lambda:
  new_lambda<-update_fda_ff_adaptive(lambda,previous_weights,learning_rate,edge_data,datamat[,(no_obs+1)])
  if((0.6<=new_lambda) & (new_lambda<1)){
    lambda <- new_lambda
  }
  print(lambda)
  lambda_over_time<-c(lambda_over_time,rep(lambda,times=frequency_of_data))
  for(no_days in c((starting_no_obs+1):(no_days_train-1))){
    print(no_days)
    no_obs<-no_days
    previous_weights <- lambda*c(previous_weights,1)
    edge_data<-datamat[,c(1:no_obs)]
    fda_mean_weight<-fda_mean_weighted(edge_data,c(1:nrow(datamat)),102,4,c(previous_weights,1),includebpoints=FALSE)
    SLD_fda_anomaly<-c(SLD_fda_anomaly,fda_mean_weight$mean_curve)
    SLD_fda_error_anomaly<-c(SLD_fda_error_anomaly,datamat[,(no_obs+1)]-fda_mean_weight$mean_curve)

    # Now Update Lambda:
    #Update lambda:
    new_lambda<-update_fda_ff_adaptive(lambda,previous_weights,learning_rate,edge_data,datamat[,(no_obs+1)])
    if((0.6<=new_lambda) & (new_lambda<1)){
      lambda <- new_lambda
    }
    lambda_over_time<-c(lambda_over_time,rep(lambda,times=frequency_of_data))
    for(i in c(((no_obs-starting_no_obs)*nrow(datamat)):((((no_obs+1-starting_no_obs)*nrow(datamat))-1)))){
      fit<-auto.arima(SLD_fda_error_anomaly[((i-(arima_window-1)):i)],max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
      fcast<-forecast(fit,h=1)
      SLD_arima_anomaly<-c(SLD_arima_anomaly,fcast$mean[1])
      SLD_combined_prediction_anomaly<-c(SLD_combined_prediction_anomaly,SLD_fda_anomaly[(i+1)-nrow(datamat)]+fcast$mean[1])
      SLD_combined_prediction_error_anomaly<-c(SLD_combined_prediction_error_anomaly,SLD_fda_anomaly[(i+1-(nrow(datamat)))]+fcast$mean[1]-datavec[(i+1+(nrow(datamat)*starting_no_obs))])
    }
  }




  for(no_days in c(no_days_train:(ncol(datamat)-1))){
    # Fit weighted FDA:
    print(no_days)
    no_obs<-no_days
    previous_weights <- lambda*c(previous_weights,1)
    edge_data<-datamat[,c(1:no_obs)]
    fda_mean_weight<-fda_mean_weighted_extreme(edge_data,c(1:nrow(datamat)),102,4,c(previous_weights,1),strange_days_sld_anomaly,strange_bins_sld_anomaly,includebpoints=FALSE)
    SLD_fda_anomaly<-c(SLD_fda_anomaly,fda_mean_weight$mean_curve)
    SLD_fda_error_anomaly<-c(SLD_fda_error_anomaly,(datamat[,(no_obs+1)]-fda_mean_weight$mean_curve))

    # Now Update Lambda:
    #Update lambda:
    new_lambda<-update_fda_ff_adaptive(lambda,previous_weights,learning_rate,edge_data,datamat[,(no_obs+1)])
    if((0.6<=new_lambda) & (new_lambda<1)){
      lambda <- new_lambda
    }
    lambda_over_time<-c(lambda_over_time,rep(lambda,times=frequency_of_data))
    # Now perform ARIMA and at each step check to see if anomaly:
    for(i in c(((no_days-starting_no_obs)*nrow(datamat)):((((no_days+1-starting_no_obs)*nrow(datamat))-1)))){
      if(length(setdiff(c((i-(arima_window-1)):i),strange_points_sld_combined_anomaly-(nrow(datamat)*starting_no_obs)))==0){
        fit<-auto.arima(SLD_fda_error_anomaly[c(i-1,i)],max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
      }
      else{
        fit<-auto.arima(SLD_fda_error_anomaly[setdiff(c((i-(arima_window-1)):i),strange_points_sld_combined_anomaly-(nrow(datamat)*starting_no_obs))],max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
      }
      fcast<-forecast(fit,h=1)
      SLD_arima_anomaly<-c(SLD_arima_anomaly,fcast$mean[1])
      # Check if anomalous point i+1 that forecasting here:
      FDA_p_value_point<-conformal_prediction_p_val(abs(SLD_fda_error_anomaly[1:(i+1)]))
      SLD_combined_prediction_anomaly<-c(SLD_combined_prediction_anomaly,SLD_fda_anomaly[(i+1-(nrow(datamat)))]+fcast$mean[1])
      SLD_combined_prediction_error_anomaly<-c(SLD_combined_prediction_error_anomaly,abs(SLD_fda_anomaly[(i+1-(nrow(datamat)))]+fcast$mean[1]-datavec[(i+1+(nrow(datamat)*starting_no_obs))]))
      Combined_p_value_point<-conformal_prediction_p_val(SLD_combined_prediction_error_anomaly)
      if(FDA_p_value_point<threshold_val){
        strange_points_sld_FDA_anomaly<-c(strange_points_sld_FDA_anomaly,(i+1+(nrow(datamat)*starting_no_obs)))
        strange_bins_sld_FDA_anomaly<-c(strange_bins_sld_FDA_anomaly,((i)%%nrow(datamat))+1)
        strange_days_sld_FDA_anomaly<-c(strange_days_sld_FDA_anomaly,no_days+1)
      }
      if(Combined_p_value_point<threshold_val){
        strange_points_sld_combined_anomaly<-c(strange_points_sld_combined_anomaly,(i+1+(nrow(datamat)*starting_no_obs)))
      }
      fisher_score <- -2*(log(FDA_p_value_point)+log(Combined_p_value_point))
      if(fisher_score>qchisq(1-threshold_val,4)){
        # Anomalous!
        strange_bins_sld_anomaly<-c(strange_bins_sld_anomaly,((i)%%nrow(datamat))+1)
        strange_days_sld_anomaly<-c(strange_days_sld_anomaly,no_days+1)
        strange_points_sld_anomaly<-c(strange_points_sld_anomaly,(i+1+(nrow(datamat)*starting_no_obs)))
      }
    }
  }

  prediction_matrix_SLD_anomaly<-data.frame(True = c(datamat[,c((no_days_train+1):ncol(datamat))]),Combined = SLD_combined_prediction_anomaly[-c(1:(nrow(datamat)*(starting_no_obs-1)))],FDA=SLD_fda_anomaly[-c(1:(nrow(datamat)*(starting_no_obs-1)))],Day=rep(c((no_days_train+1):ncol(datamat)),each=nrow(datamat)),anomalous_fisher=0,anomalous_fda=0,anomalous_combined=0)

  prediction_matrix_SLD_anomaly$anomalous_fisher[strange_points_sld_anomaly-(nrow(datamat)*no_days_train)]<-1
  prediction_matrix_SLD_anomaly$anomalous_fda[strange_points_sld_FDA_anomaly-(nrow(datamat)*no_days_train)]<-1
  prediction_matrix_SLD_anomaly$anomalous_combined[strange_points_sld_combined_anomaly-(nrow(datamat)*no_days_train)]<-1
  SLD_anomaly_times_005<-c(((nrow(datamat)*no_days_train+1)):(nrow(datamat)*ncol(datamat)))[prediction_matrix_SLD_anomaly$anomalous_fisher==1]

  return(list(prediction_matrix=prediction_matrix_SLD_anomaly,SLD_anomalous_times = SLD_anomaly_times_005,lambda_over_time=lambda_over_time))
}





#' SLD Decomposition anomaly detector initialiser
#'
#' Detects anomalies by separately detecting anomalies using long-term and combined forecast
#' from decomposition and combining their p-values using Fishers product test statistic to give
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
SLD_initial<-function(datamat,frequency_of_data,no_cycles_train,fda_forgetting_factor=0.99,arima_window_size=60,learning_rate=1e-3,threshold_val=0.005,sliding_window_size=14400){
  datavec<-c(datamat)
  fda_predictions <- c()
  arima_predictions <-c()
  strange_data_points <- c()
  FDA_AE<-c()
  ARIMA_AE<-c()

  SLD_arima_anomaly<-c()
  SLD_fda_error_anomaly <- c()
  SLD_combined_prediction_anomaly<-c()
  SLD_combined_prediction_error_anomaly<-c()
  strange_bins_sld_anomaly<-c()
  strange_days_sld_anomaly<-c()
  strange_bins_sld_FDA_anomaly<-c()
  strange_days_sld_FDA_anomaly<-c()
  strange_points_sld_anomaly<-c()
  SLD_fda_anomaly<-c()

  strange_points_sld_FDA_anomaly<-c()
  strange_points_sld_combined_anomaly<-c()

  max_no_days <- ncol(datamat)-1

  # Initialisation
  # Need to get FDA predictions for no_cycles_train cycle for ARIMA:
  no_days<-no_cycles_train-1
  # Fit weighted FDA:
  lambda<-fda_forgetting_factor
  lambda_over_time<-c(rep(lambda,times=frequency_of_data))
  no_obs<-no_days
  previous_weights <- lambda^{no_obs-c(1:(no_obs-1))}
  sum_weights <- sum(previous_weights)+1
  edge_data<-datamat[,c(1:no_obs)]
  fda_mean_weight<-fda_mean_weighted(edge_data,c(1:frequency_of_data),102,4,c(previous_weights,1),includebpoints=FALSE)
  SLD_fda_anomaly<-c(SLD_fda_anomaly,fda_mean_weight$mean_curve)

  # Now run for data up till no_cycles_train:
  no_days<-no_cycles_train
  # Fit weighted FDA:
  lambda<-fda_forgetting_factor
  lambda_over_time<-c(rep(lambda,times=frequency_of_data))
  no_obs<-no_days
  previous_weights <- lambda^{no_obs-c(1:(no_obs-1))}
  sum_weights <- sum(previous_weights)+1
  edge_data<-datamat[,c(1:no_obs)]
  fda_mean_weight<-fda_mean_weighted(edge_data,c(1:frequency_of_data),102,4,c(previous_weights,1),includebpoints=FALSE)
  SLD_fda_anomaly<-c(SLD_fda_anomaly,fda_mean_weight$mean_curve)
  #FDA_AE<-c(FDA_AE,abs(fda_mean_weight$mean_curve-datamat[,(no_obs+1)]))

  # Calculate FDA errors over cycle no_cycles_train:
  SLD_fda_error_anomaly <- c(SLD_fda_error_anomaly,abs(SLD_fda_anomaly[1:frequency_of_data]-datamat[,no_cycles_train]))


  # Now perform intial ARIMA using fda errors to get first prediction for next cycle:
  i <- (no_cycles_train-no_cycles_train+1)*frequency_of_data
  fit <- auto.arima(SLD_fda_error_anomaly[c((i-arima_window_size+1):i)],max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
  fcast <- forecast(fit,h=1)$mean[1]
  SLD_arima_anomaly <- c(SLD_arima_anomaly,fcast)
  SLD_combined_prediction_anomaly<-c(SLD_combined_prediction_anomaly,SLD_fda_anomaly[frequency_of_data+1]+fcast)

  # Prepare for updates:
  # weight_matrix <- matrix(rep(previous_weights/lambda,each=no_cycles_train),nrow=no_cycles_train,ncol=frequency_of_data)
  # weight_matrix <- weight_matrix * (1/rowSums(weight_matrix))
  # previous_FDA_wm <- rowSums((weight_matrix*datamat[,c(1:(no_cycles_train-1))]))

  weight_matrix <- matrix(rep(previous_weights/lambda,each=frequency_of_data),ncol=(no_cycles_train-1),nrow=frequency_of_data)
  weight_matrix <- weight_matrix * (1/rowSums(weight_matrix))
  previous_FDA_wm <- rowSums((weight_matrix*datamat[,c(1:(no_cycles_train-1))]))

  current_FDA_wm <- fda_mean_weight$weighted_average
  current_lambda <- lambda
  current_data <- c(datamat[,no_cycles_train])
  sum_weights_vec <- rep(sum(previous_weights/lambda),times=frequency_of_data)
  sum_weights_vec_new <- rep(sum(previous_weights)+1,times=frequency_of_data)
  new_data <-c()
  # RUN NOT ON INITIAL
  for(i in c((((no_cycles_train)*frequency_of_data)+1):((ncol(datamat)*frequency_of_data)))){
    # Assess if anomalous
    SLD_combined_prediction_anomaly<-c(SLD_combined_prediction_anomaly,SLD_fda_anomaly[i-(((no_cycles_train-1)*frequency_of_data))]+SLD_arima_anomaly[length(SLD_arima_anomaly)])
    SLD_fda_error_anomaly<-c(SLD_fda_error_anomaly,abs(SLD_fda_anomaly[i-(((no_cycles_train-1)*frequency_of_data))]-datavec[i]))
    SLD_combined_prediction_error_anomaly<-c(SLD_combined_prediction_error_anomaly,abs(SLD_combined_prediction_anomaly[i-(((no_cycles_train)*frequency_of_data))]-datavec[i]))

    # Now need to assess whether there are any anomalies
    if(length(SLD_combined_prediction_error_anomaly)>=sliding_window_size){
      FDA_p_value_point<-conformal_prediction_p_val(abs(SLD_fda_error_anomaly)[((length(SLD_fda_error_anomaly)-sliding_window_size+1):length(SLD_fda_error_anomaly))])
      Combined_p_value_point<-conformal_prediction_p_val(SLD_combined_prediction_error_anomaly[((length(SLD_combined_prediction_error_anomaly)-sliding_window_size+1):length(SLD_combined_prediction_error_anomaly))])
      if(FDA_p_value_point<threshold_val){
        strange_points_sld_FDA_anomaly<-c(strange_points_sld_FDA_anomaly,i)
      }
      if(Combined_p_value_point<threshold_val){
        strange_points_sld_combined_anomaly<-c(strange_points_sld_combined_anomaly,i)
      }
      fisher_score <- -2*(log(FDA_p_value_point)+log(Combined_p_value_point))
      if(fisher_score>qchisq(1-threshold_val,4)){
        # Anomalous!
        strange_points_sld_anomaly<-c(strange_points_sld_anomaly,i)
      }
    }


    # Run ARIMA
    arima_data_input <- SLD_fda_error_anomaly[setdiff(c((i-arima_window_size+1):i),strange_points_sld_anomaly)-(((no_cycles_train-1)*frequency_of_data))]
    if(length(arima_data_input)==0){
      while(length(arima_data_input)==0){
        arima_data_input<-SLD_fda_error_anomaly[c((i-1):(i))-(((no_cycles_train-1)*frequency_of_data))]
      }
    }
    fit<-auto.arima(arima_data_input,max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
    fcast<-forecast(fit,h=1)$mean[1]
    SLD_arima_anomaly<-c(SLD_arima_anomaly,fcast)

    new_data <- c(new_data,datavec[i])
    # If end of cycle update FDA:
    if(i%%frequency_of_data==0){
      # Put into vector which are anomalous
      if(length(strange_points_sld_anomaly)>0){
        which_anom_within_cycle <- (((strange_points_sld_anomaly[(i-frequency_of_data)<strange_points_sld_anomaly])-1)%%frequency_of_data)+1
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
      SLD_fda_anomaly<-c(SLD_fda_anomaly,fda_mean_weight$mean_curve)
      previous_FDA_wm <- current_FDA_wm
      current_FDA_wm <- fda_mean_weight$weighted_average
      sum_weights_vec_new <- fda_mean_weight$new_sum_weights_vec
      current_data <- new_data
      new_data <- c()
    }

  }


  anomalous_indicator <- rep(FALSE,length(SLD_combined_prediction_anomaly))
  anomalous_indicator[strange_data_points-(no_cycles_train*frequency_of_data)] <- TRUE

  print(c(length(datavec[-c(1:(no_cycles_train*frequency_of_data))]),length(SLD_combined_prediction_anomaly[-1]),length(SLD_fda_anomaly[-c(1:(frequency_of_data))][c(1:length(SLD_combined_prediction_anomaly))])))
  prediction_matrix_SLD_anomaly <- data.frame(True = datavec[-c(1:(no_cycles_train*frequency_of_data))],Combined = SLD_combined_prediction_anomaly[-1],FDA=SLD_fda_anomaly[-c(1:(frequency_of_data))][c(1:length(SLD_combined_prediction_anomaly)-1)],anomalous_fisher=0,anomalous_fda=0,anomalous_combined=0)

  prediction_matrix_SLD_anomaly$anomalous_fisher[strange_points_sld_anomaly-(frequency_of_data*no_cycles_train)] <- 1
  prediction_matrix_SLD_anomaly$anomalous_fda[strange_points_sld_FDA_anomaly-(frequency_of_data*no_cycles_train)] <- 1
  prediction_matrix_SLD_anomaly$anomalous_combined[strange_points_sld_combined_anomaly-(frequency_of_data*no_cycles_train)] <- 1
  strange_data_points <- strange_points_sld_anomaly

  if(length(SLD_combined_prediction_error_anomaly)>=sliding_window_size){
  return(list(prediction_matrix=prediction_matrix_SLD_anomaly,
              strange_data_points=strange_data_points,
              new_data = new_data,
              current_data = current_data,
              SLD_fda_anomaly = SLD_fda_anomaly[c((length(SLD_fda_anomaly)-frequency_of_data+1):length(SLD_fda_anomaly))],
              previous_FDA_wm = previous_FDA_wm,
              current_FDA_wm = current_FDA_wm,
              current_lambda = current_lambda,
              ARIMA_fcast = fcast,
              sum_weights_vec = sum_weights_vec,
              sum_weights_vec_new = sum_weights_vec_new,
              SLD_fda_error_anomaly = SLD_fda_error_anomaly[c(((length(SLD_fda_error_anomaly)-sliding_window_size+1)):length(SLD_fda_error_anomaly))],
              SLD_combined_prediction_error_anomaly = SLD_combined_prediction_error_anomaly[c(((length(SLD_combined_prediction_error_anomaly)-sliding_window_size+1)):length(SLD_combined_prediction_error_anomaly))]))
  }
  else{
    return(list(prediction_matrix=prediction_matrix_SLD_anomaly,
                strange_data_points=strange_data_points,
                new_data = new_data,
                current_data = current_data,
                SLD_fda_anomaly = SLD_fda_anomaly[c((length(SLD_fda_anomaly)-frequency_of_data+1):length(SLD_fda_anomaly))],
                previous_FDA_wm = previous_FDA_wm,
                current_FDA_wm = current_FDA_wm,
                current_lambda = current_lambda,
                ARIMA_fcast = fcast,
                sum_weights_vec = sum_weights_vec,
                sum_weights_vec_new = sum_weights_vec_new,
                SLD_fda_error_anomaly = SLD_fda_error_anomaly[-c(1:frequency_of_data)],
                SLD_combined_prediction_error_anomaly = SLD_combined_prediction_error_anomaly))
  }
}


#' SL Decomposition anomaly detector updater
#'
#' Detects anomalies by separately detecting anomalies using long-term and combined forecast
#' from decomposition and combining their p-values using Fishers product test statistic to give
#' final anomaly score that indicates both point and contextual anomalies. This is the funciton to update
#' existing calculations as new data becomes available.
#'
#' @param i the current data point number
#' @param new_data_point the next value in the series that is used to update the models and asssess if anomalous
#' @param current_data this is the value at the previous time point in the series
#' @param SLD_fda_anomaly this is a vector detailing the current FDA forecasts for the current cycle
#' @param previous_FDA_wm the weighted average from the previous update (at $i-2$)
#' @param current_FDA_wm the weighted average from update at time $i$
#' @param ARIMA_fcast this is the ARIMA forecast for time $i$
#' @param sum_weights_vec this is a vector detailing the weights for each time within the cycle(for the FDA model) for time i-2. This accounts for anomalous data.
#' @param sum_weights_vec_new this is a vector detailing the weights for each time within the cycle(for the FDA model) for time i-1. This accounts for anomalous data.
#' @param SLD_fda_error_anomaly the absolute error for the FDA mdoel. This is collected to run the conformal predictor and is kept at the lenght of the sliding window.
#' @param SLD_combined_prediction_error_anomaly the absolute error for the combined forecast for SLD. This is collected to run the conformal predictor and is kept at the lenght of the sliding window.
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
SLD_updater<-function(i,new_data_point,new_data,current_data,SLD_fda_anomaly,previous_FDA_wm,current_FDA_wm,ARIMA_fcast,sum_weights_vec,sum_weights_vec_new,current_lambda,SLD_fda_error_anomaly,SLD_combined_prediction_error_anomaly,strange_data_points,frequency_of_data,fda_forgetting_factor=0.99,arima_window_size=60,learning_rate=1e-3,threshold_val=0.005,sliding_window_size=14400){

  # Assess if anomalous
  FDA_fcast <- SLD_fda_anomaly[((i-1)%%frequency_of_data)+1]
  SLD_fda_error_anomaly<-c(SLD_fda_error_anomaly,abs(SLD_fda_anomaly[((i-1)%%frequency_of_data)+1]-new_data_point))
  SLD_combined_prediction_error_anomaly<-c(SLD_combined_prediction_error_anomaly,abs(SLD_fda_anomaly[((i-1)%%frequency_of_data)+1]+ARIMA_fcast-new_data_point))
  Combined_fcast <- SLD_fda_anomaly[((i-1)%%frequency_of_data)+1]+ARIMA_fcast
  # Now need to assess whether there are any anomalies
  if(length(SLD_fda_error_anomaly)>=sliding_window_size){
    FDA_p_value_point<-conformal_prediction_p_val(abs(SLD_fda_error_anomaly)[((length(SLD_fda_error_anomaly)-sliding_window_size+1):length(SLD_fda_error_anomaly))])
    Combined_p_value_point<-conformal_prediction_p_val(SLD_combined_prediction_error_anomaly[((length(SLD_combined_prediction_error_anomaly)-sliding_window_size+1):length(SLD_combined_prediction_error_anomaly))])
    # if(FDA_p_value_point<threshold_val){
    #   strange_points_sld_FDA_anomaly<-c(strange_points_sld_FDA_anomaly,i)
    # }
    # if(Combined_p_value_point<threshold_val){
    #   strange_points_sld_combined_anomaly<-c(strange_points_sld_combined_anomaly,i)
    # }
    fisher_score <- -2*(log(FDA_p_value_point)+log(Combined_p_value_point))
    if(fisher_score>qchisq(1-threshold_val,4)){
      # Anomalous!
      strange_data_points<-c(strange_data_points,i)
    }
  }

  # Run ARIMA
  arima_data_input <- SLD_fda_error_anomaly[setdiff(c((i-arima_window_size+1):i),strange_data_points)-(i-arima_window_size)]
  if(length(arima_data_input)==0){
    while(length(arima_data_input)==0){
      arima_data_input<-SLD_fda_error_anomaly[c(length(SLD_fda_error_anomaly)-1,length(SLD_fda_error_anomaly))]
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
    SLD_fda_anomaly<-c(SLD_fda_anomaly,fda_mean_weight$mean_curve)
    previous_FDA_wm <- current_FDA_wm
    current_FDA_wm <- fda_mean_weight$weighted_average
    sum_weights_vec_new <- fda_mean_weight$new_sum_weights_vec
    current_data <- new_data
    new_data <- c()
  }
  if(length(SLD_combined_prediction_error_anomaly)>=sliding_window_size){
  return(list(new_data = new_data,
              current_data = current_data,
              strange_data_points = strange_data_points,
              SLD_fda_anomaly = SLD_fda_anomaly[c((length(SLD_fda_anomaly)-frequency_of_data+1):length(SLD_fda_anomaly))],
              previous_FDA_wm = previous_FDA_wm,
              current_FDA_wm = current_FDA_wm,
              current_lambda = current_lambda,
              ARIMA_fcast = fcast,
              FDA_fcast = FDA_fcast,
              Combined_fcast = Combined_fcast,
              sum_weights_vec = sum_weights_vec,
              sum_weights_vec_new = sum_weights_vec_new,
              SLD_fda_error_anomaly = SLD_fda_error_anomaly[c(((length(SLD_fda_error_anomaly)-sliding_window_size+1)):length(SLD_fda_error_anomaly))],
              SLD_combined_prediction_error_anomaly = SLD_combined_prediction_error_anomaly[c(((length(SLD_combined_prediction_error_anomaly)-sliding_window_size+1)):length(SLD_combined_prediction_error_anomaly))]))

  }
  else{
    return(list(new_data = new_data,
                current_data = current_data,
                strange_data_points = strange_data_points,
                SLD_fda_anomaly = SLD_fda_anomaly[c((length(SLD_fda_anomaly)-frequency_of_data+1):length(SLD_fda_anomaly))],
                previous_FDA_wm = previous_FDA_wm,
                current_FDA_wm = current_FDA_wm,
                current_lambda = current_lambda,
                ARIMA_fcast = fcast,
                FDA_fcast = FDA_fcast,
                Combined_fcast = Combined_fcast,
                sum_weights_vec = sum_weights_vec,
                sum_weights_vec_new = sum_weights_vec_new,
                SLD_fda_error_anomaly = SLD_fda_error_anomaly,
                SLD_combined_prediction_error_anomaly = SLD_combined_prediction_error_anomaly))

  }
}
