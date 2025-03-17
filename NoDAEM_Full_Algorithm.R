source("DAEM_BarrierMethod_function.R")

# Reading Data
# file_name = 'Aarest_data.txt'
# file_name = 'FRT_censord.txt'
# file_name = 'RearDump.txt'
file_name = 'SerumReversal.txt'

fdata = read.table(file_name,header = T)
dataName = tools::file_path_sans_ext(file_name)

# Data preprocessing
N = nrow(fdata)
k=3 
event_vec = as.numeric(fdata[,2])
time_vec = as.numeric(fdata[,1])

tot=1e-9
maxEMIter=1e+5
maxIterAnnealing = 100
learningRateBp = 2

# annealingSchedule = seq(0.6,0.999999999,length.out=maxIterAnnealing) 
annealingSchedule = rep(1,maxIterAnnealing)
# annealingSchedule = 1 - exp(seq(log(1 - 0.7), log(1 - 0.99999), length.out = maxIterAnnealing))
# annealingSchedule = seq(0.7,0.99999999,length.out=maxIterAnnealing) # FRT
bpBaseSchedule =exp(seq(log(1e+2), log(1e+7), length.out = maxIterAnnealing))
# RearDump.txt 의 경우 bP를 아래와 같이 해야한다. 
# bpBaseSchedule =exp(seq(log(1e+1), log(1e+5), length.out = maxIterAnnealing)) # 'RearDump.txt'

## result 기록 
theta_df_full = NULL
result_latentZ_mat = NULL
## initial Parameter : beta , lambda , pi
initial_beta = c(0.5,1,2)
initial_pi_set = c(1,1,1)
initial_pi = initial_pi_set / sum(initial_pi_set)
initial_lambda = initial_lambda_calc(time_vec,event_vec,
                                     initial_beta,censored1 = ceiling(N*0.1),censored3 = ceiling(N*0.9))

# initial_lambda = initial_lambda_calc(time_vec,event_vec,
#                                      initial_beta,censored1 = ceiling(N*0.1),censored3 = ceiling(N*0.9))
## setting parameter
beta_vec = initial_beta
pi_vec = initial_pi
lambda_vec = initial_lambda
latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=1)



for( ITerAnneal in 1:maxIterAnnealing){
  annealingPara = annealingSchedule[ITerAnneal]
  bpBase  = bpBaseSchedule[ITerAnneal]
  
  
  for( iter in 1:maxEMIter){ 
    
    #### M-Step ####
    #### Pi : M-Step :  #### 
    new_pi = colSums(latentZ_mat)/N
    #### Beta : M-Step : bijection method #### 
    candi_before_vec = beta_vec
    candi_lambda_vec = lambda_vec
    
    bp1 = bpBase
    bp3 = bpBase
    
    new_beta1 = barrier_beta1(candi_before_vec[1],latentZ_mat,bp=bpBase)
    new_beta3 = barrier_beta3(candi_before_vec[3],latentZ_mat,bp=bpBase)
    
    #### Update Parameter ####
    new_beta = c(new_beta1,1,new_beta3)
    new_lambda = sapply(1:k , function(i)  sum(latentZ_mat[,i]*event_vec)/sum(latentZ_mat[,i]*(time_vec^new_beta[i])))
    
    
    parameter_diff = sqrt(sum((beta_vec-new_beta)^2+(pi_vec-new_pi)^2))
    
    ### Update Parameter record ###
    beta_vec = new_beta
    pi_vec = new_pi
    lambda_vec = new_lambda
    alpha_temper = annealingPara
    
    #### Stopping rule ####
    if(parameter_diff<tot || iter==maxEMIter){
      printResult()
      #print(c("betaTot : ",betaTot))
      print("!!!!!!!!!!!!!!!!!!!! parameter diff Break !!!!!!!!!!!!!!!")
      theta_df_full = rbind(theta_df_full,
                            data.frame(data_Name = dataName,
                                       Iter = ITerAnneal,
                                       beta1=beta_vec[1],beta3=beta_vec[3],
                                       lambda1=lambda_vec[1],lambda2=lambda_vec[2],lambda3=lambda_vec[3],
                                       diffBeta1=diffB_onlyB(beta_vec[1],latentZ_mat,j=1),
                                       diffBeta3=diffB_onlyB(beta_vec[3],latentZ_mat,j=3),
                                       diffLambda1 = diffL(beta_vec[1],lambda_vec[1],latentZ_mat,j=1),
                                       diffLambda3 = diffL(beta_vec[3],lambda_vec[3],latentZ_mat,j=3),
                                       pi1=pi_vec[1],
                                       pi2=pi_vec[2],
                                       pi3=pi_vec[3],
                                       init_beta1 = initial_beta[1],
                                       init_beta2 = initial_beta[2],
                                       init_beta3 = initial_beta[3],
                                       init_pi1   = initial_pi[1],
                                       init_pi2   = initial_pi[2],
                                       init_pi3   = initial_pi[3],
                                       EM_endIter = iter,
                                       tempering = alpha_temper,
                                       Qlike = sumQfunc(beta_vec,lambda_vec,latentZ_mat),
                                       Q1=Qfunc_onlyB(beta=beta_vec[1],latentZ_mat,j=1),
                                       Q2=Qfunc_onlyB(beta=beta_vec[2],latentZ_mat,j=2),
                                       Q3=Qfunc_onlyB(beta=beta_vec[3],latentZ_mat,j=3),
                                       bpBase = bpBase,
                                       bp1=bp1,
                                       bp3=bp3
                            ))
      
      result_latentZ_mat[[ITerAnneal]]=latentZ_mat
      break}
    
    #### E-Step ####
    latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=alpha_temper)
  }
  
}



result_Name = paste0("Proposed_Result_NoDAEM","_",dataName,".txt")
write.table(theta_df_full , file = result_Name)














