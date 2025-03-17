source("DAEM_BarrierMethod_function.R")

# Reading Data
# file_name = 'Aarest_data.txt'
file_name = 'FRT_censord.txt'
# file_name = 'RearDump.txt'
# file_name = 'SerumReversal.txt'

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

annealingSchedule = seq(0.6,0.999999999,length.out=maxIterAnnealing) 
# annealingSchedule = 1 - exp(seq(log(1 - 0.7), log(1 - 0.99999), length.out = maxIterAnnealing))
# annealingSchedule = seq(0.7,0.99999999,length.out=maxIterAnnealing) # FRT
bpBaseSchedule =exp(seq(log(1e+1), log(1e+6), length.out = maxIterAnnealing))


## result 기록 
theta_df = NULL
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



for( ITerAnneal in 1:1){
  annealingPara = 1
  # bpBase  = bpBaseSchedule[ITerAnneal]
  bpBase  = 1e+5
  
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
    
    
    theta_df = rbind(theta_df,
                          data.frame(data_Name = dataName,
                                     EM_iter = iter,
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
    result_latentZ_mat[[iter]]=latentZ_mat
    #### Stopping rule ####
    if(parameter_diff<tot || iter==maxEMIter){
      printResult()
      #print(c("betaTot : ",betaTot))
      print("!!!!!!!!!!!!!!!!!!!! parameter diff Break !!!!!!!!!!!!!!!")
      
      break}
    
    #### E-Step ####
    latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=alpha_temper)
  }
}

# colnames(theta_df) <- column_names
# result_latentZ_mat



theta_df_08 = theta_df
result_latentZ_mat_08 = result_latentZ_mat
theta_df_1 = theta_df
result_latentZ_mat_1 = result_latentZ_mat


plot(as.numeric(theta_df_08$beta1))
plot(as.numeric(theta_df_1$beta1))
# plot(as.numeric(theta_df_1$beta3))

result_latentZ_mat_08
result_latentZ_mat_1
theta_df_08

ratio1_Z1 = numeric()
for( i in 1:length(result_latentZ_mat_1)){
  ratio1_Z1[i] = result_latentZ_mat_1[[i]][1,1]/result_latentZ_mat_1[[i]][1,3]
}

result1 = cbind(theta_df_1,ratio1_Z1)
# result1 = result1[1:100,]
result1$iter =  as.numeric(result1$EM_iter )
result1$beta1 =  as.numeric(result1$beta1)

R1=ggplot(result1, aes(x = iter, y = log(ratio1_Z1, 10), color = iter)) +
  geom_point() + geom_line() +
  scale_color_gradient(low = "darkgrey", high = "red") +  # 작은 값 -> 회색, 큰 값 -> 빨간색
  labs(y = expression(log[10](Z[11]/Z[13])), x = "") +
  theme_minimal()+
  theme(legend.position="none",axis.title.y = element_text(size = 16))

R2=ggplot(result1,aes(x=iter,y=beta1, color = iter))+geom_point()+
  scale_color_gradient(low = "darkgrey", high = "red") +  # 작은 값 -> 회색, 큰 값 -> 빨간색
  labs(y =expression(beta[1]),
       x = "")+geom_line()+theme_minimal()+
  theme(legend.position="none",axis.title.y = element_text(size = 16))


result1_1 = bind_rows(result_latentZ_mat_1, .id = "iter")
result1_1$x <- rep(1:30, times = length(result_latentZ_mat_1))
result1_1$iter = as.numeric(result1_1$iter)
R3=ggplot(result1_1, aes(x = x, y = V1, group = iter, color = iter)) +
  geom_line(alpha = 0.3) +  # 선의 투명도 조절
  theme_minimal() +
  geom_hline(yintercept = 1e-5,lty=2)+
  theme(legend.position = "none") +  # legend 삭제
  labs(title = expression("Latent variable "*Z[i1]),
       x = "Index i= 1~30 ",
       y = expression(Z[i1]))+
  scale_color_gradient(
    low = "grey", high = "red",  # 그룹 숫자 작을 때 회색, 클 때 빨간색
    name = "Iteration"  # 🔥 Legend 제목 추가
  ) +
  annotate(
    "text", 
    x = 25, y = 0.02, # 특정 좌표 지정
    label = 1.8e-05, # 텍스트 내용
    size = 4, # 텍스트 크기
    color = "black", # 텍스트 색상
    # fontface = "bold" # 텍스트 스타일
  )+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.key.height = unit(1.5, "cm") 
        )

A10=(R1/R2|R3)

A10

ratio08_Z1 = numeric()
for( i in 1:length(result_latentZ_mat_08)){
  ratio08_Z1[i] = result_latentZ_mat_08[[i]][1,1]/result_latentZ_mat_08[[i]][1,3]
}


result08 = cbind(theta_df_08,ratio08_Z1)
result08$iter =  as.numeric(result08$EM_iter)
result08$beta1 =  as.numeric(result08$beta1)
result08 = result08
R1=ggplot(result08,aes(x=iter,y=log(ratio08_Z1,10), color = iter))+geom_point()+geom_line()+
  scale_color_gradient(low = "darkgrey", high = "red") +  # 작은 값 -> 회색, 큰 값 -> 빨간색
  labs(y = expression(log[10](Z[11]/Z[13])), x = "") +
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.key.height = unit(1.5, "cm"),
        axis.title.y = element_text(size = 16)
  )

# R2=ggplot(result08,aes(x=iter,y=beta1))+geom_point()+
#   labs(y =expression(beta[1]),
#        x = "Iteration")+geom_line()+theme_minimal()

R2=ggplot(result08,aes(x=iter,y=beta1, color = iter))+geom_point()+
  scale_color_gradient(low = "darkgrey", high = "red") +  # 작은 값 -> 회색, 큰 값 -> 빨간색
  labs(y =expression(beta[1]),
       x = "")+geom_line()+theme_minimal()+
  theme(legend.position="none",axis.title.y = element_text(size = 16))


result08_1 = bind_rows(result_latentZ_mat_08, .id = "iter")
result08_1$x <- rep(1:30, times = length(result_latentZ_mat_08))
result08_1$iter = as.numeric(result08_1$iter)
R3=ggplot(result08_1, aes(x = x, y = V1, group = iter, color = iter)) +
  geom_line(alpha = 0.3) +  # 선의 투명도 조절
  theme_minimal() +
  geom_hline(yintercept = 0.127,lty=2)+
  theme(legend.position = "none") +  # legend 삭제
  labs(title = expression("Latent variable "*Z[i1]),
       x = "Index i= 1~30 ",
       y = expression(Z[i1]))+
  scale_color_gradient(
    low = "dimgray", high = "red",  # 그룹 숫자 작을 때 회색, 클 때 빨간색
    name = "Iteration"  # 🔥 Legend 제목 추가
  ) +
  annotate(
    "text", 
    x = 25, y = 0.15, # 특정 좌표 지정
    label = 0.127, # 텍스트 내용
    size = 4, # 텍스트 크기
    color = "black", # 텍스트 색상
    # fontface = "bold" # 텍스트 스타일
  )+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.height = unit(1.5, "cm") 
  )


A08=(R1/R2|R3)


A08

A10+A08


result08_1 %>% tail
