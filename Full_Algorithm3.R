

library(muhaz)
library(ggplot2)
library(dplyr)

# 1. 혼합 비율 및 분포 설정
pi_exp <- runif(1, 0.6, 0.8)
pi_weib <- 1 - pi_exp
N <- 100
n_exp <- round(N * pi_exp)
n_weib <- N - n_exp

# Exponential 고장시간 생성
lambda_exp <- 1
t_exp <- rexp(n_exp, rate = lambda_exp)

# Weibull 고장시간 생성
shape_w <- 4
scale_w <-0.5
t_weib <- rweibull(n_weib, shape = shape_w, scale = scale_w)

# 결합된 고장시간 데이터
failure_times <- c(t_exp, t_weib)
min_time <- min(failure_times)
max_time <- max(failure_times)


df = data.frame(time = c(t_exp,t_weib) , model=rep(c("Exponential", "Weibull"), times = c(length(t_exp),length(t_weib))))

# 2. muhaz로 추정
haz_exp <- muhaz(df$time[df$model == "Exponential"])
haz_weib <- muhaz(df$time[df$model == "Weibull"])
haz_all <- muhaz(df$time)

# 3. data.frame으로 변환
df_exp <- data.frame(time = haz_exp$est.grid, hazard = haz_exp$haz.est, model = "Exponential")
df_weib <- data.frame(time = haz_weib$est.grid, hazard = haz_weib$haz.est, model = "Weibull")
df_all <- data.frame(time = haz_all$est.grid, hazard = haz_all$haz.est, model = "All")

# 4. 결합
df_hazards <- bind_rows(df_exp, df_weib, df_all)

# 5. ggplot 시각화
ggplot(df_hazards, aes(x = time, y = hazard, color = model)) +
  geom_line(size = 1.2) +
  labs(title = "Empirical Hazard Functions (muhaz)",
       x = "Time", y = "Estimated Hazard Rate") +
  theme_minimal()




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
df_t_grid %>% ggplot(aes(x=t_grid,h_exp))+geom_line()+geom_line(aes(x=t_grid,y=h_weib))

time_criterion = t_grid[which((h_exp/h_weib<1))[1]]<max_time*0.8 & t_grid[which((h_exp/h_weib<1))[1]]>max_time*0.3
isTRUE(time_criterion)

truePara = c(lambda_exp,shape_w,scale_w^(-shape_w),pi_exp,pi_weib,t_grid[which((h_exp/h_weib<1))[1]])
names(truePara) = c("lambda_exp","shape_w","scale_w^(-shape_w)","pi_exp","pi_weib","changePoint")

truePara %>% round(3)


################################
source("DAEM_BarrierMethod_function.R")

# fdata = df %>% mutate(event = ifelse(time>max_time*0.9,0,1)) %>% mutate(time=ifelse(event==1,time,max_time*0.9)) %>% arrange(time)
fdata = df %>% mutate(event = 1) %>% arrange(time)

dataName = "simul"
# Data preprocessing
N = nrow(fdata)
k=2 
event_vec = as.numeric(fdata$event)
time_vec = as.numeric(fdata$time)

tot=1e-9
maxEMIter=1e+5
maxIterAnnealing = 100

annealingSchedule = seq(0.1,0.999999999,length.out=maxIterAnnealing) 
bpBaseSchedule =exp(seq(log(1e+2), log(1e+4), length.out = maxIterAnnealing))

## result 기록 
theta_df_full = NULL
result_latentZ_mat = NULL
## initial Parameter : beta , lambda , pi
initial_beta = c(1,10)
initial_pi_set = c(1,1)
initial_pi = initial_pi_set / sum(initial_pi_set)

initial_lambda = initial_lambda_func3(time_vec,event_vec,
                                      initial_beta,ratio3=0.7)


## setting parameter
beta_vec = initial_beta
pi_vec = initial_pi
lambda_vec = initial_lambda
latentZ_mat = Estep_result2(beta_vec,lambda_vec,pi_vec,alpha=1)


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
    new_beta3 = barrier_beta3(candi_before_vec[2],cbind(0,latentZ_mat),bp=bpBase)
    
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


theta_df_full$diffBeta1+theta_df_full$diffBeta2
truePara

