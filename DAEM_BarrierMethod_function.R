
library(reshape2)
library(dplyr)
library(ggplot2)
library(patchwork)


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


barrierFunc_1 = function(beta,latentZ_mat,bp){
  result =  diffB_onlyB(beta, latentZ_mat, j=1)+(1/beta -1/(1-beta))*(1/bp)
  return(result)
}

barrierFunc_3 = function(beta,latentZ_mat,bp){
  result =  diffB_onlyB(beta, latentZ_mat, j=3)+log(beta-1))/bp
  return(result)
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

find_zero_new_beta1 <- function(latentZ_mat, j = 1, interval = c(0, 1)) {
  result <- uniroot(function(new_beta1) diffB_onlyB(new_beta1, latentZ_mat, j),
  interval = interval,tol=1e-10)
  return(result$root)
}

find_zero_new_beta3 <- function(latentZ_mat, j = 3, interval = c(1, 200)) {
  result <- uniroot(function(new_beta3) diffB_onlyB(new_beta3, latentZ_mat, j),
  interval = interval,tol=1e-10)
  return(result$root)
}


barrier_beta1_NR = function(beta,latentZ_mat,j,bp){
  stepsize = 1e-3
  if( diffB_onlyB(1,latentZ_mat,j=1)<0){
    beta_candi = find_zero_new_beta1(latentZ_mat, j=1)
  }else{
    while(TRUE){
      gradient = diffB_onlyBeta1_with_Barrier(beta,latentZ_mat,j,bp)
      beta_candi = beta + stepSize_func(gradient,beta,stepsize)
      if(abs(diffB_onlyB(beta_candi,latentZ_mat,j=1))<1e-3){break}else{
        stepsize = 0.1*stepsize
      }
      if(beta_candi>1){
        # print("!!!!!!!!! Beta1 Over 1 !!!!!!!!!!")
        beta_candi = (beta+0.999)/2}
      if(beta_candi<0){beta_candi = beta/2}
      if(abs(beta_candi-beta)<tot){break}
      beta = beta_candi
    }
  }
  return(beta_candi)
}


barrier_beta3_NR = function(beta,latentZ_mat,j,bp){
  stepsize = 1e-3
  if( diffB_onlyB(1,latentZ_mat,j=3)>0){
    beta_candi = find_zero_new_beta3(latentZ_mat, j=3)
  }else{
    while(TRUE){
      gradient = diffB_onlyBeta3_with_Barrier(beta,latentZ_mat,j,bp)
      beta_candi = beta + stepSize_func(gradient,beta,stepsize)
      
      if(abs(diffB_onlyB(beta_candi,latentZ_mat,j=3))<1e-3){break}else{
        stepsize = 0.1*stepsize
      }
      
      if(beta_candi<1){beta_candi = (1+beta)/2}
      if(abs(beta_candi-beta)<tot){break}
      beta = beta_candi
    }
  }
  return(beta_candi)
}

stepSize_func = function(grad,beta,step_size){
  grad_size = abs(grad)
  if(grad_size<step_size){
    grad_step = grad
  }else{
    while(grad_size>step_size){
      grad_size = grad_size/10
    }
    grad_step = sign(grad)*grad_size
  }
  return(grad_step)
}

barrier_beta1 = function(beta,latentZ_mat,bp){
  result <- uniroot(function(beta) barrierFunc_1(beta,latentZ_mat,bp),
  interval = c(0,1),tol=1e-10)
  return(result$root)
}


barrier_beta3 = function(beta,latentZ_mat,bp){
  result <- uniroot(function(beta) barrierFunc_3(beta,latentZ_mat,bp),
  interval = c(1, beta*2),tol=1e-10)
  return(result$root)
}
















