
source("DAEM_BarrierMethod_function.R")

library(fitdistrplus)
library(dplyr)
n_total = 50
n_const_prop = 0.5
n_wear_prop = 0.2
n_const <- n_total*n_const_prop
n_wear <- n_total*n_wear_prop
shape_w = 78
lambda_w = 2.18e-151
scale_w = lambda_w^(-1/shape_w)
# lambda_w = scale_w^(-shape_w)
lambda_e = 0.0252
scale_e = 1/lambda_e
time <- c(rweibull(n_const, shape = 1, scale = scale_e),
          rweibull(n_wear, shape = shape_w, scale = scale_w))
event <- rep(1, length(time))
df <- data.frame(time = time, event = event)
df$distn = c(rep("exp",n_const),rep("wear",n_wear))
wearMax =df %>% filter(distn=="wear") %>% pull(time) %>% max
df = df %>% filter(time<wearMax)
df = df %>% arrange(time)

df$hazard = compute_empirical_hazard_surv(df)



# hazard point plot
ggplot(df, aes(x = time, y = hazard ,color=distn )) +
  geom_point() +
  geom_line(color = "blue") +
  labs(title = "Empirical Hazard Rate (전체 데이터)",
       x = "Time",
       y = "Hazard Rate") +
  theme_minimal()

# ggplot(df, aes(x = time, y = weibull_fit )) +
#   geom_point() +
#   geom_line(color = "blue")



fdata = df %>% arrange(time)

dataName = "simul"
# Data preprocessing
N = nrow(fdata)
k=2 
event_vec = as.numeric(fdata$event)
time_vec = as.numeric(fdata$time)


########################## 



#######################

tot=1e-9
maxEMIter=1e+5
maxIterAnnealing = 100

annealingSchedule = seq(0.1,0.999999999,length.out=maxIterAnnealing) 
bpBaseSchedule =exp(seq(log(1e-1), log(1e+7), length.out = maxIterAnnealing))

## result 기록 
theta_df_full = NULL
result_latentZ_mat = NULL
## initial Parameter : beta , lambda , pi
# initial_beta = c(1,35)
initial_beta = c(1,30)
initial_pi_set = c(1,1)
initial_pi = initial_pi_set / sum(initial_pi_set)

initial_lambda = initial_lambda_func3(time_vec,event_vec,
                                      initial_beta,ratio3=0.9)

# initial_lambda = c(0.01287726 ,3.066099e-69 )
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

trueDF = data.frame(beta2=1,lambda2=lambda_e ,beta3=shape_w,lambda3=lambda_w 
                    ,pi2=n_const/(n_const+n_wear) , pi3= n_wear/(n_const+n_wear))
trueDF

df$hz2 = hazardrate(df$time,1,optData$lambda1)
df$hz3 = hazardrate(df$time,optData$beta3,optData$lambda3)

df %>% ggplot(aes(x=time,y=hazard ))+geom_point()+
  geom_line(aes(x=time,y=hz2),color="blue")+geom_line(aes(x=time,y=hz3),color="red")






