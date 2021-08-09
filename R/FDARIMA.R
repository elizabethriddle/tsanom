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
#' @param p the significance level to use for the prediction interval
#'
#' @return List containing average, sd and upper and lower ranges.
#' @export
fda_mean_weighted_poisson<-function(datamat,rangex,nb,no,obsweights,includebpoints=FALSE,bpoints=NULL,p=0.999){
  n<-ncol(datamat)
  data_mean<-rowWeightedMeans(datamat,obsweights)
  if(includebpoints!=FALSE){
    basis_function <- create.bspline.basis(rangeval=c(min(rangex),max(rangex)),nbasis=nb,norder=no,breaks=bpoints)
  }
  else{
    basis_function <- create.bspline.basis(rangeval=c(min(rangex),max(rangex)),nbasis=nb,norder=no)
  }
  basis_mat<-eval.basis(rangex,basis_function)
  basis_data<-as.data.frame(basis_mat)
  basis_data_mat<-cbind(data_mean,basis_data)
  basis_coeffs<-glm(data_mean~ . -1,family="poisson",data=basis_data_mat)

  avmat<-basis_coeffs$fitted.values

  # Weighted SD:
  smoothdatamat<- as.matrix(basis_coeffs$fitted.values)[,rep(1,times=n)]
  datares<-datamat - smoothdatamat
  w_sum<-sum(obsweights)
  w2_sum<-sum(obsweights^2)
  datavar<-apply((datares^2)*obsweights,1,sum)/(w_sum-(w2_sum/w_sum))
  edge.stderr<-sqrt(datavar)

  # Now find the lower and upper and upper limits at
  q<-0.5+(p/2)
  number_points <- rep(ncol(datamat),times=nrow(datamat))
  mult_val <- qt(q,number_points-1)

  lower_val <- avmat - (mult_val*edge.stderr*sqrt(1+(1/number_points)))
  upper_val <- avmat + (mult_val*edge.stderr*sqrt(1+(1/number_points)))


  # Now return the average values, and the standard deviation values:
  return(list(mean_curve=avmat,st_curve=edge.stderr,upper_val=upper_val,lower_val=lower_val))
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
#' @param p the significance level to use for the prediction interval
#'
#' @return List containing average, sd and upper and lower ranges.
#' @export
fda_mean_weighted_poisson_extreme<-function(datamat,rangex,nb,no,obsweights,extreme_days,extreme_bins,includebpoints=FALSE,bpoints=NULL,p=0.999){
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
  basis_coeffs<-glm(data_mean~ . -1,family="poisson",data=basis_data_mat)

  avmat<-basis_coeffs$fitted.values


  # Weighted SD:
  smoothdatamat<- as.matrix(basis_coeffs$fitted.values)[,rep(1,times=n)]
  datares<-datamat - smoothdatamat
  datares[weight_matrix==0]<-0
  w_sum<-sum(obsweights)
  w2_sum<-sum(obsweights^2)
  datavar<-apply((datares^2)*obsweights,1,sum)/(w_sum-(w2_sum/w_sum))
  edge.stderr<-sqrt(datavar)

  # Now find the lower and upper and upper limits at
  q<-0.5+(p/2)
  number_points <- rowSums(weight_matrix>0)
  mult_val <- qt(q,number_points-1)

  lower_val <- avmat - (mult_val*edge.stderr*sqrt(1+(1/number_points)))
  upper_val <- avmat + (mult_val*edge.stderr*sqrt(1+(1/number_points)))

  # Now return the average values, and the standard deviation values:
  return(list(mean_curve=avmat,st_curve=edge.stderr,upper_val=upper_val,lower_val=lower_val))
}

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
FDARIMA<-function(datamat,no_days_train,max_no_days,fda_forgetting_factor=0.86,arima_window_size=60,learning_rate=1e-3,threshold_val=0.005){
  datavec<-c(datamat)
  poisson_fda_predictions <- c()
  poisson_fda_variances <- c()
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
      lambda_over_time<-c(rep(lambda,times=1440))
      w_vec<-c()
      no_obs<-no_days
      for(i in 1:no_obs){
        w_vec<-c(w_vec,((lambda)^(no_obs-i)))
      }
      edge_data<-datamat[,c(1:no_obs)]
      fda_mean_weight_poisson<-fda_mean_weighted_poisson(edge_data,c(1:1440),102,4,w_vec,includebpoints=FALSE)
      poisson_fda_predictions<-c(poisson_fda_predictions,fda_mean_weight_poisson$mean_curve)
      poisson_fda_variances<-c(poisson_fda_variances,fda_mean_weight_poisson$st_curve)
      FDA_AE<-c(FDA_AE,abs(fda_mean_weight_poisson$mean_curve-datamat[,(no_obs+1)]))

      #Update lambda:
      new_lambda<-update_fda_ff(lambda,learning_rate,edge_data,datamat[,(no_obs+1)])
      if((0.6<=new_lambda) & (new_lambda<1)){
        lambda <- new_lambda
        lambda_over_time<-c(lambda_over_time,rep(new_lambda,times=1440))
      }else{
        lambda_over_time<-c(lambda_over_time,rep(lambda,times=1440))
      }

      # Now perform ARIMA and at each step check to see if anomaly:
      for(i in c(((no_days)*1440):((((no_days+1)*1440)-1)))){
        fit<-auto.arima(datavec[setdiff(c((i-arima_window_size+1):i),strange_data_points)],max.p = 10,max.q = 10, max.d = 5,seasonal = FALSE)
        fcast<-forecast(fit,h=1)$mean[1]
        arima_predictions<-c(arima_predictions,fcast)
        ARIMA_AE<-c(ARIMA_AE,abs(fcast-datavec[(i+1)]))
      }
    }
    else{
      # Fit poisson weighted FDA:
      w_vec<-c()
      no_obs<-no_days
      for(i in 1:no_obs){
        w_vec<-c(w_vec,((lambda)^(no_obs-i)))
      }
      edge_data<-datamat[,c(1:no_obs)]
      fda_mean_weight_poisson<-fda_mean_weighted_poisson_extreme(edge_data,c(1:1440),102,4,w_vec,strange_days,strange_bins,includebpoints=FALSE)
      poisson_fda_predictions<-c(poisson_fda_predictions,fda_mean_weight_poisson$mean_curve)
      poisson_fda_variances<-c(poisson_fda_variances,fda_mean_weight_poisson$st_curve)
      FDA_AE<-c(FDA_AE,abs(fda_mean_weight_poisson$mean_curve-datamat[,(no_obs+1)]))

      if(no_days !=max_no_days){
        #Update lambda:
        new_lambda<-update_fda_ff(lambda,learning_rate,edge_data,datamat[,(no_obs+1)])
        if((0.6<=new_lambda) & (new_lambda<1)){
          lambda <- new_lambda
          lambda_over_time<-c(lambda_over_time,rep(new_lambda,times=1440))
        }else{
          lambda_over_time<-c(lambda_over_time,rep(lambda,times=1440))
        }
      }
      # Now perform ARIMA and at each step check to see if anomaly:
      for(i in c(((no_days)*1440):((((no_days+1)*1440)-1)))){
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
          strange_bins<-c(strange_bins,((i)%%1440)+1)
          if(strange_bins==1){
            strange_days<-c(strange_days,no_days+1)
          }
          else{
            strange_days<-c(strange_days,no_days)
          }
          strange_data_points<-c(strange_data_points,i+1)
        }

      }

    }
  }
  anomalous_indicator <- rep(FALSE,length(ARIMA_AE))
  anomalous_indicator[strange_data_points-(no_days_train*1440)] <- TRUE
  prediction_matrix<-data.frame(True = datavec[-c(1:(no_days_train*1440))],ARIMA = arima_predictions,ARIMA_AE = ARIMA_AE,FDA=poisson_fda_predictions,FDA_AE = FDA_AE,Anomalous=anomalous_indicator,Lambda_val = lambda_over_time)

  return(list(prediction_matrix=prediction_matrix,strange_data_points=strange_data_points))
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
