
source("DAEM_BarrierMethod_function.R")


# Reading Data

file_name = 'LEDLife_114C_30A.txt'
fdata = read.table(file_name,header = T)
dataName = tools::file_path_sans_ext(file_name)

# Data preprocessing
N = nrow(fdata)
k=2 
event_vec = as.numeric(fdata[,2])
time_vec = as.numeric(fdata[,1])



########################## 



#######################

tot=1e-9
maxEMIter=1e+5
maxIterAnnealing = 100

annealingSchedule = seq(0.1,0.999999999,length.out=maxIterAnnealing) 
bpBaseSchedule =exp(seq(log(1e-1), log(1e+5), length.out = maxIterAnnealing))

## result 기록 
theta_df_full = NULL
result_latentZ_mat = NULL
## initial Parameter : beta , lambda , pi
# initial_beta = c(1,35)
initial_beta = c(1,10)
initial_pi_set = c(1,1)
initial_pi = initial_pi_set / sum(initial_pi_set)

initial_lambda = initial_lambda_func3(time_vec,event_vec,
                                      initial_beta,ratio3=0.9)



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
plot(theta_df_full$diffBeta1 %>% abs)

plot(theta_df_full$diffBeta1 %>% abs+theta_df_full$diffBeta3 %>% abs)

plot(theta_df_full$pi1 %>% abs)

optData = theta_df_full[which.min(theta_df_full$diffBeta3 %>% abs),]
optData


hazard_em = compute_empirical_hazard_surv(fdata)

df_with_hazard <- fdata %>%
  left_join(hazard_em %>% select(time, hazard) %>% rename(hazard_em = hazard),
            by = "time")

df_with_hazard$hz2 = hazardrate(df_with_hazard$time,1,optData$lambda1)
df_with_hazard$hz3 = hazardrate(df_with_hazard$time,optData$beta3,optData$lambda3)

head(df_with_hazard)
df_with_hazard %>% ggplot(aes(x=time,y=hazard_em ))+
  geom_point()+
  geom_point(aes(x=time,y=hz2),color="red")+
  geom_point(aes(x=time,y=hz3),color="blue")



#########
time_DBlist = seq(min(time_vec),max(time_vec),by=0.1)

# DBbound1 = DecisionBoundary(time_DBlist,beta_vec = c(optBeta1,1,optBeta3),lambda_vec = c(optLambda1,optLambda2,optLambda3),pi_vec=c(optPi1,optPi2,optPi3),j=1)
DBbound3 = DecisionBoundary(time_DBlist,beta_vec = c(1,1,optData$beta3),lambda_vec = c(1,optData$lambda1,optData$lambda3),pi_vec=c(1,optData$pi1,optData$pi3),j=3)

# changePoint1 = time_DBlist[which.min(DBbound1>1)]
changePoint3 = time_DBlist[which.min(DBbound3<1)]
# DBdata = data.frame(time = time_DBlist, postProbRatio1 =DBbound1, postProbRatio3 = DBbound3 )
DBdata = data.frame(time = time_DBlist, postProbRatio3 = DBbound3 )
# DB1 = DBdata %>% ggplot(aes(x=time,y=postProbRatio1))+geom_line()+  theme_minimal()+
#   geom_hline(yintercept = 1,lty=2)+
#   annotate("text", x = changePoint1+10, y = 2,label=paste0("Time :",changePoint1))+
#   labs(
#     title = "Posterior Ratio (infant vs constant)",
#     x = "",
#     y = ""
#   )+ scale_y_continuous(
#     breaks = c(1,seq(2, 7, by = 5)) # y축 눈금을 5 단위로 설정
#   )+
  # geom_point(data = data.frame(x = c(changePoint1), y = c(1)), aes(x = x, y = y), color = "blue", size = 3)
DB3 = DBdata %>% ggplot(aes(x=time,y=postProbRatio3))+geom_line()+ 
  labs(
    title = "Posterior Ratio (wear-out vs constant)",
    x = "",
    y = ""
  ) +
  annotate("text", x = changePoint3, y = 1,label=paste0("Time :",changePoint3))+
  geom_point(data = data.frame(x = c(changePoint3), y = c(1)), aes(x = x, y = y), color = "blue", size = 3)+
  theme_minimal()+geom_hline(yintercept = 1,lty=2)


DB3

