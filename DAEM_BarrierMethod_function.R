library(reshape2)
library(dplyr)
library(ggplot2)
library(patchwork)
library(survival)
library(ggbreak)
library(cowplot)
library(gridExtra)
library(grid)
library(scales)
library(foreach)
library(doParallel)
library(rootSolve)



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


hazardrate = function(t,beta,lambda){
  beta*lambda *t^(beta-1)
}

weibull_func = function(t,beta,lambda){
  lambda*beta*t^(beta-1)*exp(-lambda*(t^beta))
}




diffB_onlyB = function(beta,latentZ_mat,j){
  sum(latentZ_mat[,j]*event_vec)/beta + 
    sum(latentZ_mat[,j]*event_vec*log(time_vec))-
    sum(latentZ_mat[,j]*event_vec)*sum(latentZ_mat[,j]*(time_vec^beta)*log(time_vec))/sum(latentZ_mat[,j]*(time_vec^beta))
}

diffL = function(beta,lambda,latentZ_mat,j){
  sum(latentZ_mat[,j]*(event_vec/lambda - time_vec^beta))
}

DecisionBoundary = function(t,beta_vec,lambda_vec,j=1){
  (pi_vec[j]/pi_vec[2])*(beta_vec[j]/beta_vec[2])*(lambda_vec[j]/lambda_vec[2])*t^(beta_vec[j]-1)*exp(-lambda_vec[j]*t^{beta_vec[j]}+lambda_vec[2]*t)
}




newton_raphson <- function(beta_init,lambda,latentZ_mat, j, tol = 1e-6, max_iter = 100) {
  beta <- beta_init
  for (i in 1:max_iter) {
    # 함수값 f(beta)
    f_beta <- Qfunc(beta ,lambda, latentZ_mat,j)
    df_beta <- -diffB_onlyB(beta,latentZ_mat,j)
    
    print(df_beta)
    # 뉴턴-랩슨 업데이트
    beta_new <- beta - f_beta / df_beta
    
    print(f_beta / df_beta)
    # 수렴 확인
    if (abs(beta_new - beta) < tol) {
      return(beta_new)
    }
    
    beta <- beta_new
  }
  
  warning("Maximum iterations reached without convergence")
  return(beta)
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

# Estep_result = function(beta,lambda,pi_vec,alpha=1){
#   # 각 값의 개수 계산
#   value_counts <-as.vector(table(time_vec))
#   # # 원본 벡터와 같은 순서로 각 값의 개수 적용
# 
#   # beta = beta_vec
#   # lambda = lambda_vec
#   # alpha=1
#   wpdf1=(weibull_func(unique(time_vec),beta[1],lambda[1])*pi_vec[1])^value_counts
#   wpdf2=(weibull_func(unique(time_vec),beta[2],lambda[2])*pi_vec[2])^value_counts
#   wpdf3=(weibull_func(unique(time_vec),beta[3],lambda[3])*pi_vec[3])^value_counts
#   wpdf1 = wpdf1^alpha
#   wpdf2 = wpdf2^alpha
#   wpdf3 = wpdf3^alpha
#   sumWeibull = wpdf1+wpdf2+wpdf3
#   
# 
#   # 사후확률 계산 (고유한 시간 값 기준)
#   posterior_df = data.frame(
#     V1 = wpdf1 / sumWeibull,
#     V2 = wpdf2 / sumWeibull,
#     V3 = wpdf3 / sumWeibull
#   )
#   
#   # 중복된 시간값만큼 행을 반복하여 확장 (각 행을 value_counts 만큼 반복)
#   expanded_df <- posterior_df[rep(1:nrow(posterior_df), times = value_counts), ]
#   
# }

# Estep_result2 = function(beta,lambda,pi_vec,alpha=1){
#   wpdf1=weibull_func(time_vec,beta[1],lambda[1])*pi_vec[1]
#   wpdf2=weibull_func(time_vec,beta[2],lambda[2])*pi_vec[2]
#   # wpdf3=weibull_func(time_vec,beta[3],lambda[3])*pi_vec[3]
#   wpdf1 = wpdf1^alpha
#   wpdf2 = wpdf2^alpha
#   # wpdf3 = wpdf3^alpha
#   sumWeibull = wpdf1+wpdf2
#   data.frame(V1=wpdf1/sumWeibull,
#              V2=wpdf2/sumWeibull)
#   # V3=wpdf3/sumWeibull)
# }



Estep_result2 = function(beta,lambda,pi_vec,alpha=1){
  # unq_time_vec = unique(time_vec)  
  # 각 값의 개수 계산
  value_counts <- table(time_vec)
  # 원본 벡터와 같은 순서로 각 값의 개수 적용
  count_vec <-   as.numeric(sapply(time_vec, function(x) value_counts[as.character(x)]))
  wpdf1=weibull_func(time_vec,beta[1],lambda[1])*pi_vec[1]
  wpdf2=weibull_func(time_vec,beta[2],lambda[2])*pi_vec[2]
  # wpdf3=weibull_func(time_vec,beta[3],lambda[3])*pi_vec[3]
  wpdf1 = (wpdf1^alpha)/count_vec
  wpdf2 = (wpdf2^alpha)/count_vec
  # wpdf3 = wpdf3^alpha
  sumWeibull = (wpdf1+wpdf2)
  data.frame(V1=(wpdf1/sumWeibull),
             V2=(wpdf2/sumWeibull)  )
  # V3=wpdf3/sumWeibull)
}

barrierFunc_1 = function(beta,latentZ_mat,bp,barrierExist=1){
  result =  diffB_onlyB(beta, latentZ_mat, j=1)+barrierExist*(1/beta -1/(1-beta))*(1/bp)
  return(result)
}

barrierFunc_3 = function(beta,latentZ_mat,bp,barrierExist=1){
  result =  diffB_onlyB(beta, latentZ_mat, j=3)+barrierExist*(1/bp)*(1/(beta-1))
  # result =  diffB_onlyB(beta, latentZ_mat, j=3)+barrierExist*(1/bp)*(1/(beta-1)-1/(500-beta))
  return(result)
}

barrier_beta1 = function(beta,latentZ_mat,bp){
  result <- uniroot(function(beta) barrierFunc_1(beta,latentZ_mat,bp),
  interval = c(0,1),tol=1e-10)
  return(result$root)
}

barrier_beta3 = function(beta,latentZ_mat,bp){
  maxRange=beta
  # while(!is.na(diffB_onlyB(maxRange, latentZ_mat, j=3))&& !is.na(barrierFunc_3(maxRange,latentZ_mat,bp))){
  while(!is.na(barrierFunc_3(maxRange,latentZ_mat,bp))){
  maxRange=maxRange+1
  }
  # print(maxRange)
  # print(bp)
  # print(barrierFunc_3(maxRange,latentZ_mat,bp))
  # while (barrierFunc_3(maxRange,latentZ_mat,bp)<0 && !is.na(barrierFunc_3(maxRange,latentZ_mat,bp)) ) {
  #   bp = bp*1.1
  #   print(bp)
  # }
  maxRange = maxRange - 1
  result = uniroot(function(beta) barrierFunc_3(beta,latentZ_mat,bp),
  interval = c(1, maxRange),tol=1e-10)
  return(result$root)
}

# nobarrier_beta1 = function(beta,latentZ_mat,bp){
#   result <- uniroot(function(beta) barrierFunc_1(beta,latentZ_mat,bp,barrierExist=0),
#                     interval = c(0.000001,1),tol=1e-10)
#   return(result$root)
# }

nobarrier_beta1 = function(beta,latentZ_mat,bp){
  result <- uniroot(function(beta) diffB_onlyB(beta, latentZ_mat, j=1),
                    interval = c(0,1),tol=1e-10)
  return(result$root)
}

nobarrier_beta3 = function(beta,latentZ_mat,bp){
  maxRange=beta
  while(!is.na(diffB_onlyB(maxRange, latentZ_mat, j=3))){
    maxRange=maxRange+1
  }
  maxRange = maxRange - 1
  result = uniroot(function(beta) barrierFunc_3(beta,latentZ_mat,bp,barrierExist=0),
                   interval = c(1, maxRange),tol=1e-10)
  return(result$root)
}



findBeta1 = function(){
  beta <- tryCatch({
    # nobarrier_beta1 실행
    nobarrier_beta1(candi_before_vec[1], latentZ_mat, bp = bp1)
  }, error = function(e) {
    # 오류 발생 시 barrier_beta1 실행
    barrier_beta1(candi_before_vec[1], latentZ_mat, bp = bp1)
  })
  
  diff <- tryCatch({
    # nobarrier_beta1 실행
    betaIn = nobarrier_beta1(candi_before_vec[1], latentZ_mat, bp = bp1)
    diffB_onlyB(beta, latentZ_mat, j=1)
  }, error = function(e) {
    # 오류 발생 시 barrier_beta1 실행
    betaIn = barrier_beta1(candi_before_vec[1], latentZ_mat, bp = bp1)
    barrierFunc_1(beta,latentZ_mat,bp1)
  })
  
  result = list(beta=beta ,diff=diff)
  return(result)
}

findBeta3 = function(){
  beta <- tryCatch({
    # nobarrier_beta1 실행
    nobarrier_beta3(candi_before_vec[3], latentZ_mat, bp = bp3)
  }, error = function(e) {
    # 오류 발생 시 barrier_beta1 실행
    barrier_beta3(candi_before_vec[3], latentZ_mat, bp = bp3)
  })
  
  diff <- tryCatch({
    # nobarrier_beta1 실행
    betaIn = nobarrier_beta3(candi_before_vec[3], latentZ_mat, bp = bp3)
    diffB_onlyB(beta, latentZ_mat, j=3)
  }, error = function(e) {
    # 오류 발생 시 barrier_beta1 실행
    betaIn = barrier_beta3(candi_before_vec[3], latentZ_mat, bp = bp3)
    barrierFunc_3(beta,latentZ_mat,bp3)
  })

  result = list(beta=beta ,diff=diff)
  return(result)
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


initial_lambda_calc2 = function(time_vec,event_vec,beta_vec,censored1,censored3){
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
  
  # # censored3 시간 벡터 및 위험률 추출
  # time_censored3 <-unique(time_vec)[censored1:length(cum_hazard)]
  # time_transformed3 <- beta_vec[3]*time_censored3^(beta_vec[3]-1)
  # # 위험률 벡터의 검열된 값 추출
  # hazard_censored3 <- hazard_rate[censored1:length(cum_hazard)]
  # # 선형 회귀 실행
  # lm_fit3 <- lm(hazard_censored3 ~ time_transformed3)
  
  
  
  ## initial lambda
  # lambda_vec = c(abs(lm_fit1$coefficients[2]),abs(mean(hazard_rate[censored1:censored3],na.rm=T)),abs(lm_fit3$coefficients[2]) )%>% as.vector()
  lambda_vec = c(abs(lm_fit1$coefficients[2]),abs(mean(hazard_rate[censored1:censored3],na.rm=T)) )%>% as.vector()
  return(lambda_vec)
}



printResult = function(){
  print("#####################################################################################################################")
  print( paste0( "EM iteration : " ,ITerAnneal ," / " , iter ," / "," sumQ :",sumQfunc(beta_vec,lambda_vec,latentZ_mat)))
  print( paste0( "Annealing : " ,annealingPara))
  print(paste0(c("pi_vec : " , sapply(pi_vec , function(i) round(i,2))),collapse = " / "))
  print(paste0(c("Lambda : " , sapply(lambda_vec , function(i) round(i,4))),collapse = " / "))
  print(paste0(c("Beta :",sapply(beta_vec , function(i) round(i,4))) ,collapse = " / "))
  print(paste0("Beta diff : " , diffB_onlyB(beta_vec[1],latentZ_mat,j=1) %>% abs+ diffB_onlyB(beta_vec[3],latentZ_mat,j=3) %>% abs))
  print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))
  print(paste0(c("Beta1 at 3 is positive :",diffB_onlyB(1,latentZ_mat,j=3)>0,diffB_onlyB(1,latentZ_mat,j=3)) ,collapse = " / "))
  # print(paste0(" data save : ", nrow(theta_df)))
  print(paste0(" bpBase : ", bpBase ))
  print(paste0(" Init Beta : ",initial_beta))
  print(paste0(" Init Pi : ", initial_pi ))
  print(paste0(" bp1 , bp3 : ", bp1 , "/", bp3 ))
  # print(paste0(" betaTot : ", betaTot )) 
  print("#####################################################################################################################")
}

printResult2 = function(){
  print("#####################################################################################################################")
  print( paste0( "EM iteration : " ,ITerAnneal ," / " , iter ))
  print( paste0( "Annealing : " ,annealingPara))
  print(paste0(c("pi_vec : " , sapply(pi_vec , function(i) round(i,2))),collapse = " / "))
  print(paste0(c("Lambda : " , sapply(lambda_vec , function(i) round(i,4))),collapse = " / "))
  print(paste0(c("Beta :",sapply(beta_vec , function(i) round(i,4))) ,collapse = " / "))
  print(paste0("Beta diff : " , diffB_onlyB(beta_vec[1],latentZ_mat,j=1) %>% abs))
  print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))
  # print(paste0(c("Beta1 at 3 is positive :",diffB_onlyB(1,latentZ_mat,j=3)>0,diffB_onlyB(1,latentZ_mat,j=3)) ,collapse = " / "))
  # print(paste0(" data save : ", nrow(theta_df)))
  print(paste0(" bpBase : ", bpBase ))
  print(paste0(" Init Beta : ",initial_beta))
  print(paste0(" Init Pi : ", initial_pi ))
  # print(paste0(" bp1 , bp3 : ", bp1 , "/", bp3 ))
  # print(paste0(" betaTot : ", betaTot )) 
  print("#####################################################################################################################")
}

# generate_latentZ_mat = function(n) {
#   if (n < 2) stop("n은 2 이상이어야 합니다.")
#   
#   # 첫 번째 열: 0.5에서 0.1로 선형 감소
#   col1 = seq(10, 0.01, length.out = n)
#   
#   # 두 번째 열: 0.25에서 시작, 중간에서 0.9로 증가했다가 다시 0.25로 감소
#   mid_point = ceiling(n / 2)  # 중간 지점
#   col2_first_half = seq(1, 10, length.out = mid_point)
#   col2_second_half = seq(10, 1, length.out = n - mid_point)
#   col2 = c(col2_first_half, col2_second_half)
#   
#   # 세 번째 열: 1 - col1 - col2
#   col3 = 10 - col1 - col2
#   sumCol = abs(col1)+abs(col2)+abs(col3)
#   col1 = abs(col1)/sumCol
#   col2 = abs(col2)/sumCol
#   col3 = abs(col3)/sumCol
#   # 행렬 생성
#   latentZ_mat = cbind(col1, col2, col3)
#   
#   # 유효성 검사: 각 행의 합이 1인지 확인
#   if (!all(abs(rowSums(latentZ_mat) - 1) < 1e-10)) {
#     stop("행렬의 각 행의 합이 1이 아닙니다.")
#   }
#   
#   return(latentZ_mat)
# }




generate_latentZ_mat = function(n, lambda = 1, scale_output = TRUE) {
  if (n < 2) stop("n은 2 이상이어야 합니다.")
  
  # col1: exp 함수로 감소
  col1 = exp(-lambda * seq(0, 1, length.out = n))
  
  # col2: 중간 행의 값을 일정하게 유지
  mid_index = ceiling(n / 2)
  col2 = rep(col1[mid_index], n)
  
  # col3: col1의 역순
  col3 = rev(col1)
  
  # 행렬 생성
  latentZ_mat = cbind(col1, col2, col3)
  
  # scale로 각 행의 합을 1로 조정
  if (scale_output) {
    row_sums = rowSums(latentZ_mat)
    latentZ_mat = latentZ_mat / row_sums
  }
  
  return(latentZ_mat)
}













