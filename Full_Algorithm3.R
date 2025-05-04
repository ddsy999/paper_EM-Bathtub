
source("DAEM_BarrierMethod_function.R")
library(muhaz)
library(ggplot2)
library(dplyr)

# 1. 혼합 비율 및 분포 설정
pi_exp <- runif(1, 0.4, 0.6)
pi_weib <- 1 - pi_exp
N <- 100
n_exp <- round(N * pi_exp)
n_weib <- N - n_exp





target_mode <- 320
target_var <- 3500

loss_fn <- function(par) {
  shape <- par[1]
  scale <- par[2]
  
  if (shape <= 1 || scale <= 0) return(1e6)  # mode가 존재하려면 shape > 1
  
  mode <- scale * ((shape - 1) / shape)^(1 / shape)
  var <- scale^2 * (gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2)
  
  (mode - target_mode)^2 + (var - target_var)^2
}

# 초기값
res <- optim(c(2, 200), loss_fn, method = "L-BFGS-B", lower = c(1.01, 0.1))

shape_w <- res$par[1]
scale_w <- res$par[2]

# 확인
mode_weib <- scale_w * ((shape_w - 1) / shape_w)^(1 / shape_w)
var_weib <- scale_w^2 * (gamma(1 + 2 / shape_w) - gamma(1 + 1 / shape_w)^2)

cat("shape_w =", shape_w, "\nscale_w =", scale_w, "\n")
cat("mode =", mode_weib, "\nvariance =", var_weib, "\n")

















# Weibull 고장시간 생성
# shape_w <- 35.0267 #30*(1+0.1*runif(1,0,1))
# scale_w <- 228.9428  #5*(1+0.3*runif(1,0,1))
lambda_w = scale_w^(-shape_w)
t_weib <- rweibull(n_weib, shape = shape_w, scale = scale_w)
mean_weib <- scale_w * gamma(1 + 1 / shape_w)
mode_weib = scale_w * ((shape_w - 1) / shape_w)^(1 / shape_w)
var_weib <- scale_w^2 * (gamma(1 + 2 / shape_w) - (gamma(1 + 1 / shape_w))^2)


# Exponential 고장시간 생성
lambda_exp <- (1/(600))*(1+0.05*runif(1,0,1))
t_exp <- rexp(n_exp, rate = lambda_exp)
mean_exp = 1/lambda_exp


# 결합된 고장시간 데이터
failure_times <- c(t_exp, t_weib)
min_time <- min(failure_times)
max_time <- max(failure_times)


DBbound3 = DecisionBoundary(seq(min_time, max_time, length.out = 500),beta_vec = c(1,1,shape_w),
                            lambda_vec = c(1,lambda_exp,lambda_w),pi_vec=c(1,pi_exp,pi_weib),j=3)
changePoint3 = seq(min_time, max_time, length.out = 500)[which.min(DBbound3<1)]



df = data.frame(time = c(t_exp,t_weib) , model=rep(c("Exponential", "Weibull"), times = c(length(t_exp),length(t_weib))))
df = df %>% filter(time<300)

trueDF = data.frame(beta2=1,beta3=shape_w ,
                    lambda1=lambda_exp, lambda3=lambda_w,
                    pi2=pi_exp,pi3=pi_weib
)



###############
df = df %>% mutate(event = ifelse(time<mode_weib,1,0)) 
surv_obj <- Surv(time = df$time, event = df$event)
fit <- survfit(surv_obj ~ 1)

# 시간과 누적 생존율
times <- fit$time
surv_probs <- fit$surv

# 누적 hazard (대략적 추정)
cumhaz <- -log(surv_probs)

# 시간 구간별 변화량
delta_time <- diff(c(0, times))
delta_hazard <- diff(c(0, cumhaz))
hazard_rate <- delta_hazard / delta_time

df$hazard_rate= hazard_rate 
###############


# 2. 시간축: 실제 데이터 범위 기반
t_grid <- seq(min_time, max_time, length.out = 500)

# 3. 이론적 hazard 계산
# Exponential
f_exp <- lambda_exp * exp(-lambda_exp * t_grid)
S_exp <- exp(-lambda_exp * t_grid)
h_exp <- rep(lambda_exp, length(t_grid))

# Weibull
f_weib <- (shape_w / scale_w) * (t_grid / scale_w)^(shape_w - 1) * exp(-(t_grid / scale_w)^shape_w)
S_weib <- exp(-(t_grid / scale_w)^shape_w)
h_weib <- (shape_w / scale_w) * (t_grid / scale_w)^(shape_w - 1)

df_t_grid = data.frame(time=t_grid ,h_exp, h_weib) 

################################


# fdata = df %>% mutate(event = ifelse(time>max_time*0.9,0,1)) %>% mutate(time=ifelse(event==1,time,max_time*0.9)) %>% arrange(time)

fdata = df %>% arrange(time)

dataName = "simul"
# Data preprocessing
N = nrow(fdata)
k=2 
event_vec = as.numeric(fdata$event)
time_vec = as.numeric(fdata$time)

###### TTT 

# 고장 데이터만 정렬
failures <- sort(fdata$time[fdata$event == 1])
n <- length(failures)
T_total <- sum(failures)

# TTT 좌표 계산
x_vals <- (1:n) / n
y_vals <- sapply(1:n, function(i) {
  (sum(failures[1:i]) + (n - i) * failures[i]) / T_total
})

ttt_df <- data.frame(x = x_vals, y = y_vals)

# TTT plot
TTT2=ggplot(ttt_df, aes(x = x, y = y)) +
  geom_step(direction = "hv") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Device G",
    x = "Proportion of Failures (i/n)",
    y = "TTT Transform"
  ) +
  theme_minimal()

TTT2+df %>% ggplot(aes(x=time,y=hazard_rate,color=model))+geom_point()+df_t_grid %>% ggplot(aes(x=t_grid,h_exp))+geom_line()+geom_line(aes(x=t_grid,y=h_weib),color="red")

c(mean_exp,mean_weib)
mode_weib
var_weib
changePoint3






#######################

tot=1e-9
maxEMIter=1e+5
maxIterAnnealing = 100

annealingSchedule = seq(0.1,0.999999999,length.out=maxIterAnnealing) 
bpBaseSchedule =exp(seq(log(1e+1), log(1e+7), length.out = maxIterAnnealing))

## result 기록 
theta_df_full = NULL
result_latentZ_mat = NULL
## initial Parameter : beta , lambda , pi
initial_beta = c(1,30)
initial_pi_set = c(1,1)
initial_pi = initial_pi_set / sum(initial_pi_set)

initial_lambda = initial_lambda_func3(time_vec,event_vec,
                                      initial_beta,ratio3=0.5)

# initial_lambda = c(0.0101530363319033,5.10733273080284e-11)
# 


# initial_lambda = c(1e-34,1e-40)
# initial_lambda = c(2.096286 ,9.332562e-12)


## setting parameter
beta_vec = initial_beta
pi_vec = initial_pi
lambda_vec = initial_lambda
latentZ_mat = Estep_result2(beta_vec,lambda_vec,pi_vec,alpha=0.1)

latentZ_mat %>% round(5)


print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))
print(paste0(c("Beta3 at 1 is positive :",diffB_onlyB(1,latentZ_mat,j=2)>0,diffB_onlyB(1,latentZ_mat,j=2)) ,collapse = " / "))




for( ITerAnneal in 1:maxIterAnnealing){
  annealingPara = annealingSchedule[ITerAnneal]
  bpBase  = bpBaseSchedule[ITerAnneal]
  
  
  for( iter in 1:maxEMIter){ 
    
    #### M-Step ####
    #### Pi : M-Step :  #### 
    new_pi = colSums(latentZ_mat)/N
    # new_pi = colSums(latentZ_mat %>% distinct())/sum(latentZ_mat %>% distinct())
    #### Beta : M-Step : bijection method #### 
    candi_before_vec = beta_vec
    candi_lambda_vec = lambda_vec
    
    # bp1 = bpBase
    bp3 = bpBase
    
    # new_beta1 = barrier_beta1(candi_before_vec[1],latentZ_mat,bp=bpBase)
    new_beta3 = barrier_beta3(candi_before_vec[2],latentZ_mat=cbind(0,latentZ_mat),bp=bpBase)
    # print(new_beta3)
    #### Update Parameter ####
    new_beta = c(1,new_beta3)
    new_lambda = sapply(1:k , function(i)  sum(latentZ_mat[,i]*event_vec)/sum(latentZ_mat[,i]*(time_vec^new_beta[i])))
    
    
    parameter_diff = sqrt(sum((beta_vec-new_beta)^2+(pi_vec-new_pi)^2))
    
    ### Update Parameter record ###
    beta_vec = new_beta
    pi_vec = new_pi
    lambda_vec = new_lambda
    alpha_temper = annealingPara
    
    #### Stopping rule ####
    if(parameter_diff<tot || iter==maxEMIter){
      printResult3()
      #print(c("betaTot : ",betaTot))
      print("!!!!!!!!!!!!!!!!!!!! parameter diff Break !!!!!!!!!!!!!!!")
      theta_df_full = rbind(theta_df_full,
                            data.frame(data_Name = dataName,
                                       Iter = ITerAnneal,
                                       # beta1=beta_vec[1],
                                       beta3=beta_vec[2],
                                       lambda1=lambda_vec[1],
                                       lambda3=lambda_vec[2],
                                       # lambda3=lambda_vec[3],
                                       diffBeta1=diffB_onlyB(beta_vec[1],latentZ_mat,j=1),
                                       diffBeta3=diffB_onlyB(beta_vec[2],latentZ_mat,j=2),
                                       # diffBeta3=diffB_onlyB(beta_vec[3],latentZ_mat,j=3),
                                       diffLambda1 = diffL(beta_vec[1],lambda_vec[1],latentZ_mat,j=1),
                                       diffLambda3 = diffL(beta_vec[2],lambda_vec[2],latentZ_mat,j=2),
                                       pi1=pi_vec[1],
                                       pi3=pi_vec[2],
                                       # pi3=pi_vec[3],
                                       init_beta1 = initial_beta[1],
                                       init_beta3 = initial_beta[2],
                                       # init_beta3 = initial_beta[3],
                                       init_pi1   = initial_pi[1],
                                       init_pi3   = initial_pi[2],
                                       # init_pi3   = initial_pi[3],
                                       EM_endIter = iter,
                                       tempering = alpha_temper,
                                       # Qlike = sumQfunc(beta_vec,lambda_vec,latentZ_mat),
                                       Q1=Qfunc_onlyB(beta=beta_vec[1],latentZ_mat,j=1),
                                       Q3=Qfunc_onlyB(beta=beta_vec[2],latentZ_mat,j=2),
                                       # Q3=Qfunc_onlyB(beta=beta_vec[3],latentZ_mat,j=3),
                                       bpBase = bpBase,
                                       # bp1=bp1
                                       bp3=bp3
                            ))
      
      result_latentZ_mat[[ITerAnneal]]=latentZ_mat
      break}
    
    #### E-Step ####
    latentZ_mat = Estep_result2(beta_vec,lambda_vec,pi_vec,alpha=alpha_temper)
  }
  
}


plot(theta_df_full$diffBeta3 %>% abs)

optData = theta_df_full[which.min(theta_df_full$diffBeta3 %>% abs),]
optData
trueDF


df$hz2 = hazardrate(df$time,1,optData$lambda1)
df$hz3 = hazardrate(df$time,optData$beta3,optData$lambda3)

df %>% ggplot(aes(x=time,y=hazard_rate))+geom_point()+
  geom_line(aes(y=hz2),color="blue")+geom_line(aes(y=hz3),color="red")

