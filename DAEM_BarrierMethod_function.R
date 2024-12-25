library(reshape2)
library(dplyr)
library(ggplot2)
library(patchwork)

if (!requireNamespace("rootSolve", quietly = TRUE)) {
  install.packages("rootSolve")
}
library(rootSolve)

theta_df_full = NULL
####### Denote function ####

Qfunc = function(beta = 1 ,lambda, latentZ_mat, j=1){
  sum(latentZ_mat[,j]*event_vec*log(lambda)) +
    sum(latentZ_mat[,j]*event_vec*log(beta)) +
    sum(latentZ_mat[,j]*event_vec*(beta-1)*log(time_vec)) -
    sum(latentZ_mat[,j]*lambda*time_vec^beta) +
    sum(latentZ_mat[,j]*log(pi_vec[j]))
}

Qfunc_onlyB = function(beta = 1 ,latentZ_mat, j=1){
  sum(latentZ_mat[,j]*event_vec)*log(sum(latentZ_mat[,j]*event_vec)/sum(latentZ_mat[,j]*(time_vec^beta)))+
    sum(latentZ_mat[,j]*event_vec*log(beta))+
    sum(latentZ_mat[,j]*event_vec*(beta-1)*log(time_vec))-
    sum(latentZ_mat[,j]*event_vec)
}

sumQfunc = function(beta_vec,lambda_vec,latentZ_mat ){
  Qfunc(beta=beta_vec[1],lambda_vec[1],latentZ_mat,j=1)+Qfunc(beta=beta_vec[2],lambda_vec[2],latentZ_mat,j=2)+Qfunc(beta=beta_vec[3],lambda_vec[3],latentZ_mat,j=3)
}


# hazardrate = function(t,beta,lambda){
#   beta*lambda *t^(beta-1)
# }

weibull_func = function(t,beta,lambda){
  lambda*beta*t^(beta-1)*exp(-lambda*(t^beta))
}




diffB_onlyB = function(beta,latentZ_mat,j){
  sum(latentZ_mat[,j]*event_vec)/beta + 
    sum(latentZ_mat[,j]*event_vec*log(time_vec))-
    sum(latentZ_mat[,j]*event_vec)*sum(latentZ_mat[,j]*(time_vec^beta)*log(time_vec))/sum(latentZ_mat[,j]*(time_vec^beta))
}



DecisionBoundary = function(t,beta_vec,lambda_vec,j=1){
  (pi_vec[j]/pi_vec[2])*(beta_vec[j]/beta_vec[2])*(lambda_vec[j]/lambda_vec[2])*t^(beta_vec[j]-1)*exp(-lambda_vec[j]*t^{beta_vec[j]}+lambda_vec[2]*t)
}


Estep_result = function(beta,lambda,pi_vec,alpha=1){
  wpdf1=weibull_func(time_vec,beta[1],lambda[1])*pi_vec[1]
  wpdf2=weibull_func(time_vec,beta[2],lambda[2])*pi_vec[2]
  wpdf3=weibull_func(time_vec,beta[3],lambda[3])*pi_vec[3]
  wpdf1 = wpdf1^alpha
  wpdf2 = wpdf2^alpha
  wpdf3 = wpdf3^alpha
  sumWeibull = wpdf1+wpdf2+wpdf3
  data.frame(V1=wpdf1/sumWeibull,
             V2=wpdf2/sumWeibull,
             V3=wpdf3/sumWeibull)
}


barrierFunc_1 = function(beta,latentZ_mat,bp){
  result =  diffB_onlyB(beta, latentZ_mat, j=1)+(1/beta -1/(1-beta))*(1/bp)
  return(result)
}

barrierFunc_3 = function(beta,latentZ_mat,bp){
  result =  diffB_onlyB(beta, latentZ_mat, j=3)+(1/(beta-1)+1/(beta-600))*(1/bp)
  return(result)
}

barrier_beta1 = function(beta,latentZ_mat,bp){
  result <- uniroot(function(beta) barrierFunc_1(beta,latentZ_mat,bp),
  interval = c(0,1),tol=1e-10)
  return(result$root)
}

barrier_beta3 = function(beta,latentZ_mat,bp){
  maxRange=beta
  while(!is.na(diffB_onlyB(maxRange, latentZ_mat, j=3))){
  maxRange=maxRange+1
  }
  maxRange = maxRange - 1
  result = uniroot(function(beta) barrierFunc_3(beta,latentZ_mat,bp),
  interval = c(1, maxRange),tol=1e-10)
  return(result$root)
}




initial_lambda_calc = function(time_vec,event_vec,beta_vec,censored1,censored3){
  library(survival)
  surv_obj <- Surv(time_vec, event_vec)
  km_fit <- survfit(surv_obj ~ 1)
  cum_hazard <- cumsum(km_fit$n.event / km_fit$n.risk)
  hazard_rate <- diff(cum_hazard) / diff(km_fit$time)
  # plot(hazard_rate)

  # censored1 = floor((length(time_vec))/10)+1
  # censored3 = floor((length(time_vec))*0.6)
  # censored1까지의 시간 벡터 및 위험률 추출
  time_censored1 <-unique(time_vec)[1:censored1]
  # beta_vec[1]을 사용하여 시간 벡터를 제곱
  time_transformed1 <- beta_vec[1]*time_censored1^(beta_vec[1]-1)
  # 위험률 벡터의 검열된 값 추출
  hazard_censored1 <- hazard_rate[1:censored1]
  # 선형 회귀 실행
  lm_fit1 <- lm(hazard_censored1 ~ time_transformed1)

  # censored3 시간 벡터 및 위험률 추출
  time_censored3 <-unique(time_vec)[censored1:length(cum_hazard)]
  time_transformed3 <- beta_vec[3]*time_censored3^(beta_vec[3]-1)
  # 위험률 벡터의 검열된 값 추출
  hazard_censored3 <- hazard_rate[censored1:length(cum_hazard)]
  # 선형 회귀 실행
  lm_fit3 <- lm(hazard_censored3 ~ time_transformed3)



  ## initial lambda
  lambda_vec = c(abs(lm_fit1$coefficients[2]),abs(mean(hazard_rate[censored1:censored3],na.rm=T)),abs(lm_fit3$coefficients[2]) )%>% as.vector()
  return(lambda_vec)
}























