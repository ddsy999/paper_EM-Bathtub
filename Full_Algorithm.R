source("DAEM_BarrierMethod_function.R")

# Reading Data
# file_name = 'Aarest_data.txt'
# file_name = 'RearDump.txt'
file_name = 'SerumReversal.txt'
# file_name = 'FRT_censord.txt'
fdata = read.table(file_name,header = T)

# Data preprocessing
N = nrow(fdata)
k=3 
event_vec = fdata[,2] %>% as.numeric()
time_vec = fdata[,1]%>% as.numeric()
# time_vec = time_vec+1
# time_vec = time_vec/(max(time_vec)*1.1)
tot=1e-7
maxIter=20000
# maxIterAnnealing = 1000000
maxIterAnnealing = 100
learningRate = 1.1
learningRateBp = 2
bpBase  = 0.1
initialbpBase = 1e+1
# betaTot = 1e-4
initalAnnealingPara = 0.6
# 0.991597
# annealingSchedule = 1-1/(seq(2,12,length.out=maxIterAnnealing)^2)
# annealingSchedule = 1-1/(seq(50,100,length.out=maxIterAnnealing))
annealingSchedule = seq(0.7,0.9999,length.out=maxIterAnnealing)
bpBaseSchedule = seq(initialbpBase,1e+4,length.out=maxIterAnnealing)
betaTotSchedule = seq(1e-2,1e-8,length.out=maxIterAnnealing)



# 11개의 열을 가진 빈 result data.frame 생성
column_names <- c("beta1", "lambda1", "beta2", "lambda2", "beta3", "lambda3", "sumQfunc","diffB_beta1","diffB_beta3","bp","iter","Beta1 at 1","Q1","Q2","Q3","pi1","pi2","pi3")
theta_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(theta_df) <- column_names
result_latentZ_mat = list()

## initial Parameter : beta , lambda , pi
# initial_beta = c(0.5,1,5)
initial_beta = c(0.5,1,2)
initial_pi_set = c(1,1,1)
initial_pi = initial_pi_set / sum(initial_pi_set)


initial_lambda = initial_lambda_calc(time_vec,event_vec,
                                     initial_beta,censored1 = 20,censored3 = 25)
initial_lambda
beta_vec = initial_beta
pi_vec = initial_pi

# initial lambda를 아래처럼 한다면?
# init_latentmat = generate_latentZ_mat(N,lambda =0.01)
# initial_lambda = sapply(1:k , function(i)  sum(init_latentmat[,i]*event_vec)/sum(init_latentmat[,i]*(time_vec^initial_beta[i])))
lambda_vec = initial_lambda
latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=1)

# 
# print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))
# print(paste0(c("Beta1 at 3 is positive :",diffB_onlyB(1,latentZ_mat,j=3)>0,diffB_onlyB(1,latentZ_mat,j=3)) ,collapse = " / "))




for( ITerAnneal in 1:maxIterAnnealing){
  

  annealingPara = annealingSchedule[ITerAnneal]
  bpBase  = bpBaseSchedule[ITerAnneal]
  betaTot = betaTotSchedule[ITerAnneal]


  # print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))
  # print(paste0(c("Beta1 at 3 is positive :",diffB_onlyB(1,latentZ_mat,j=3)>0,diffB_onlyB(1,latentZ_mat,j=3)) ,collapse = " / "))
  # 
  
  bp1=bpBase
  bp3=bpBase
  
for( iter in 1:maxIter){ 

  #### M-Step ####
  #### Pi : M-Step :  #### 
  new_pi = colSums(latentZ_mat)/N
  #### Beta : M-Step : bijection method #### 
  candi_before_vec = beta_vec
  candi_lambda_vec = lambda_vec
  
  
  ###############
  # print(" ### Beta1 Update ###")
  # bp1=bpBase
  # for( i in 1:10000){
  #   beta1Result = findBeta1()
  #   new_beta1 = beta1Result$beta
  #   if(new_beta1==1){
  #     print("########## Beta 1 ==1 Break")
  #     break}
  #   if(abs(beta1Result$diff)<betaTot){break}
  #   bp1=bp1*learningRate
  # }
  # 
  # bp3=bpBase
  # for( i in 1:10000){
  #   # beta3Result = findBeta3()
  #   new_beta3 = beta3Result$beta
  #   if(abs(beta3Result$diff)<betaTot){break}
  #   bp3=bp3*learningRate
  # }
  

  # bp1=bpBase
  # for( i in 1:10000){
  #   bp1=bp1*learningRate
  #   new_beta1 = barrier_beta1(candi_before_vec[1],latentZ_mat,bp=bp1 )
  #   if(abs(barrierFunc_1(new_beta1,latentZ_mat,bp1))<betaTot){break}
  #   # if(abs(diffB_onlyB(new_beta1,latentZ_mat,j=1))<betaTot){break}
  #   if(bp1>1e+200){break}
  # }
  # bp3=bpBase
  # for( i in 1:10000){
  #   bp3=bp3*learningRate
  #   new_beta3 = barrier_beta3(candi_before_vec[3],latentZ_mat,bp=bp3 )
  #   if(abs(barrierFunc_3(new_beta3,latentZ_mat,bp3))<betaTot){break}
  #   # if(abs(diffB_onlyB(new_beta3,latentZ_mat,j=3))<betaTot){break}
  # 
  #   if(bp3>1e+200){break}
  # }
  
  
  bp1=bp1*learningRate
  new_beta1 = barrier_beta1(candi_before_vec[1],latentZ_mat,bp=bp1 )
  bp3=bp3*learningRate
  new_beta3 = barrier_beta3(candi_before_vec[3],latentZ_mat,bp=bp3 )

  
  #### Update Parameter ####
  new_beta = c(new_beta1,1,new_beta3)
  new_lambda = sapply(1:k , function(i)  sum(latentZ_mat[,i]*event_vec)/sum(latentZ_mat[,i]*(time_vec^new_beta[i])))
  
  
  # print(paste0("differential : ", diffB_onlyB(new_beta[1],latentZ_mat,j=1)," ",diffB_onlyB(new_beta[3],latentZ_mat,j=3)))
  
  
  parameter_diff = sqrt(sum((beta_vec-new_beta)^2+(pi_vec-new_pi)^2))
  
  
  ### Update Parameter record ###
  beta_vec = new_beta
  pi_vec = new_pi
  lambda_vec = new_lambda
  alpha_temper = annealingPara
  
  #### Stopping rule ####
  # if(parameter_diff<1e-10 || iter==maxIter){
  if(parameter_diff<tot || iter==maxIter){
    printResult()
    print("!!!!!!!!!!!!!!!!!!!! parameter diff Break !!!!!!!!!!!!!!!")
    theta_df_full = rbind(theta_df_full,
                          data.frame(beta1=beta_vec[1],beta3=beta_vec[3],
                                     lambda1=lambda_vec[1],lambda3=lambda_vec[3],
                                     diffBeta1=diffB_onlyB(beta_vec[1],latentZ_mat,j=1),
                                     diffBeta3=diffB_onlyB(beta_vec[3],latentZ_mat,j=3),
                                     pi1=pi_vec[1],
                                     pi2=pi_vec[2],
                                     pi3=pi_vec[3],
                                     init_beta1 = initial_beta[1],
                                     init_beta2 = initial_beta[2],
                                     init_beta3 = initial_beta[3],
                                     init_pi1   = initial_pi[1],
                                     init_pi2   = initial_pi[2],
                                     init_pi3   = initial_pi[3],
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
  # if(diffB_onlyB(1,latentZ_mat,j=1)>0){
  #   print( c("diffB_onlyB(1,latentZ_mat,j=1)>0" , diffB_onlyB(1,latentZ_mat,j=1)))
  #   break
  # }
}


