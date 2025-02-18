source("DAEM_BarrierMethod_function.R")

# Reading Data
# file_name = 'Aarest_data.txt'
file_name = 'FRT_censord.txt'
fdata = read.table(file_name,header = T)

# Data preprocessing
N = nrow(fdata)
k=3 
event_vec = fdata[,2] %>% as.numeric()
time_vec = fdata[,1]%>% as.numeric()
time_vec = time_vec+1
# time_vec = time_vec/(max(time_vec)*1.1)
tot=1e-8
maxIter=3000
learningRate = 1.1
learningRateBp = 2

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
                                     initial_beta,censored1 = ceiling(N*0.3),censored3 = ceiling(N*0.7))

bpBase  = 1e+2
betaTot = 1e-4

beta_vec = initial_beta
pi_vec = initial_pi
lambda_vec = initial_lambda
## initial latent variable 
latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=1)
alpha_temper = 1

diffB_onlyB(1,latentZ_mat ,j=1)

for( iter in 1:maxIter){ 
  print("#####################################################################################################################")
  print( paste0( "EM iteration : " , iter ," sumQ :",sumQfunc(beta_vec,lambda_vec,latentZ_mat)))
  print(paste0(c("pi_vec : " , sapply(pi_vec , function(i) round(i,2))),collapse = " / "))
  print(paste0(c("Lambda : " , sapply(lambda_vec , function(i) round(i,2))),collapse = " / "))
  print(paste0(c("Beta :",sapply(beta_vec , function(i) round(i,4))) ,collapse = " / "))
  print(paste0("Beta diff : " , diffB_onlyB(beta_vec[1],latentZ_mat,j=1) %>% abs+ diffB_onlyB(beta_vec[3],latentZ_mat,j=3) %>% abs))
  print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))
  print(paste0(c("Beta1 at 3 is positive :",diffB_onlyB(1,latentZ_mat,j=3)>0,diffB_onlyB(1,latentZ_mat,j=3)) ,collapse = " / "))
  print(paste0(" data save : ", nrow(theta_df)))
  print("#####################################################################################################################")
  
  #### M-Step ####
  #### Pi : M-Step :  #### 
  new_pi = colSums(latentZ_mat)/N
  #### Beta : M-Step : bijection method #### 
  candi_before_vec = beta_vec
  candi_lambda_vec = lambda_vec
  
  
  bp1=bpBase
  for( i in 1:10000){
    new_beta1 = barrier_beta1(candi_before_vec[1],latentZ_mat,bp=bp1 )
    bp1=bp1*learningRateBp
    if(bp1>1e+8){break}
  }
  # barrier_beta1
  # barrierFunc_1
  bp3=bpBase
  for( i in 1:10000){
    new_beta3 = barrier_beta3(candi_before_vec[3],latentZ_mat,bp=bp3 )
    bp3=bp3*learningRateBp
    if(bp3>1e+8){break}
  }
  
  # diffB_onlyB(candi_before_vec[3], latentZ_mat, j=3)+(1/(candi_before_vec[3]-1)+1/(candi_before_vec[3]-500))*(1/bp3)
  # diffB_onlyB(candi_before_vec[3], latentZ_mat, j=3)+(1/(candi_before_vec[3]-1)+1/(candi_before_vec[3]-500))*(1/bp3)
  # diffB_onlyB(candi_before_vec[3], latentZ_mat, j=3)
  # diffB_onlyB(candi_before_vec[1], latentZ_mat, j=1)
  
  ################
  
  #### Update Parameter ####
  new_beta = c(new_beta1,1,new_beta3)
  new_lambda = sapply(1:k , function(i)  sum(latentZ_mat[,i]*event_vec)/sum(latentZ_mat[,i]*(time_vec^new_beta[i])))
  
  
  print(paste0("differential : ", diffB_onlyB(new_beta[1],latentZ_mat,j=1)," ",diffB_onlyB(new_beta[3],latentZ_mat,j=3)))
  
  
  parameter_diff = sqrt(sum((beta_vec-new_beta)^2+(pi_vec-new_pi)^2))
  
  
  ### Update Parameter record ###
  beta_vec = new_beta
  pi_vec = new_pi
  lambda_vec = new_lambda
  
  theta_df = rbind(theta_df,
                   c(beta_vec[1],lambda_vec[1],beta_vec[2],lambda_vec[2],beta_vec[3],lambda_vec[3],sumQfunc(beta_vec,lambda_vec,latentZ_mat)
                     ,diffB_onlyB(beta_vec[1],latentZ_mat,j=1) ,
                     diffB_onlyB(beta_vec[3],latentZ_mat,j=3) ,
                     "bp",iter,
                     diffB_onlyB(1,latentZ_mat,j=1),
                     Qfunc_onlyB(beta=beta_vec[1],latentZ_mat,j=1),
                     Qfunc_onlyB(beta=beta_vec[2],latentZ_mat,j=2),
                     Qfunc_onlyB(beta=beta_vec[3],latentZ_mat,j=3),
                     pi_vec[1],
                     pi_vec[2],
                     pi_vec[3]
                   ))
  
  
  result_latentZ_mat[[iter]]=latentZ_mat
  
  #### Stopping rule ####
  
  if(parameter_diff<1e-6){
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
                                     tempering = alpha_temper
                          ))
    
    result_latentZ_mat[[iter]]=latentZ_mat
    break}
  
  #### E-Step ####
  latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=alpha_temper)
  
}

colnames(theta_df) <- column_names
result_latentZ_mat



theta_df_08 = theta_df
result_latentZ_mat_08 = result_latentZ_mat
theta_df_1 = theta_df
result_latentZ_mat_1 = result_latentZ_mat


plot(as.numeric(theta_df_08$beta1))
plot(as.numeric(theta_df_1$beta1))
plot(as.numeric(theta_df_1$beta3))

result_latentZ_mat_08
result_latentZ_mat_1
theta_df_08

ratio1_Z1 = numeric()
for( i in 1:length(result_latentZ_mat_1)){
  ratio1_Z1[i] = result_latentZ_mat_1[[i]][1,1]/result_latentZ_mat_1[[i]][1,3]
}

result1 = cbind(theta_df_1,ratio1_Z1)
# result1 = result1[1:100,]
result1$iter =  as.numeric(result1$iter)
result1$beta1 =  as.numeric(result1$beta1)
R1=ggplot(result1,aes(x=iter,y=log(ratio1_Z1,10)))+geom_point()+geom_line()+
  labs(y = expression(log[10](Z[11]/ Z[13])),
       x = "Iteration")+theme_minimal()
R2=ggplot(result1,aes(x=iter,y=beta1))+geom_point()+
  labs(y =expression(beta[1]),
       x = "Iteration")+geom_line()+theme_minimal()


result1_1 = bind_rows(result_latentZ_mat_1, .id = "iter")
result1_1$x <- rep(1:30, times = length(result_latentZ_mat_1))
result1_1$iter = as.numeric(result1_1$iter)
R3=ggplot(result1_1, aes(x = x, y = V1, group = iter, color = iter)) +
  geom_line(alpha = 0.3) +  # 선의 투명도 조절
  scale_color_gradient(low = "grey", high = "red") +  # 🔥 그룹 숫자 작을 때 회색, 클 때 빨간색
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
  theme(legend.position = "right")

A10=(R1/R2|R3)+plot_annotation(title = expression("Annealing parameter " * gamma == 1), 
                           theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))


A10

ratio08_Z1 = numeric()
for( i in 1:length(result_latentZ_mat_08)){
  ratio08_Z1[i] = result_latentZ_mat_08[[i]][1,1]/result_latentZ_mat_08[[i]][1,3]
}


result08 = cbind(theta_df_08,ratio08_Z1)
result08$iter =  as.numeric(result08$iter)
result08$beta1 =  as.numeric(result08$beta1)
result08 = result08
R1=ggplot(result08,aes(x=iter,y=log(ratio08_Z1,10)))+geom_point()+geom_line()+
  labs(y = expression(log[10](Z[11]/ Z[13])))+theme_minimal()
R2=ggplot(result08,aes(x=iter,y=beta1))+geom_point()+
  labs(y =expression(beta[1]),
       x = "Iteration")+geom_line()+theme_minimal()


result08_1 = bind_rows(result_latentZ_mat_08, .id = "iter")
result08_1$x <- rep(1:30, times = length(result_latentZ_mat_08))
result08_1$iter = as.numeric(result08_1$iter)
R3=ggplot(result08_1, aes(x = x, y = V1, group = iter, color = iter)) +
  geom_line(alpha = 0.3) +  # 선의 투명도 조절
  scale_color_gradient(low = "grey", high = "red") +  # 🔥 그룹 숫자 작을 때 회색, 클 때 빨간색
  theme_minimal() +
  geom_hline(yintercept = 0.127,lty=2)+
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
    x = 25, y = 0.11, # 특정 좌표 지정
    label = 0.127, # 텍스트 내용
    size = 4, # 텍스트 크기
    color = "black", # 텍스트 색상
    # fontface = "bold" # 텍스트 스타일
  )+
  theme(legend.position = "right")

A08=(R1/R2|R3)+plot_annotation(title = expression("Annealing parameter " * gamma == 0.8), 
                           theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

A10
A08

A10/A08


result08_1 %>% tail
