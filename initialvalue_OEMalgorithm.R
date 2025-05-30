source("DAEM_BarrierMethod_function.R")


# Reading Data
# file_name = 'Aarest_data.txt'
# file_name = 'RearDump.txt'
# file_name = 'SerumReversal.txt'
file_name = 'FRT_censord.txt'
fdata = read.table(file_name,header = T)
theta_df_EMfull = data.frame()
# Data preprocessing
N = nrow(fdata)
k=3 
event_vec = fdata[,2] %>% as.numeric()
time_vec = fdata[,1]%>% as.numeric()
# time_vec = time_vec/(max(time_vec)*1.1)
tot=1e-8
maxIter = 2000

# 11개의 열을 가진 빈 result data.frame 생성
column_names <- c("beta1", "lambda1", "beta2", "lambda2", "beta3", "lambda3", "sumQfunc","diffB_beta1","diffB_beta3","bp","iter","Beta1 at 1","Q1","Q2","Q3","pi1","pi2","pi3")
theta_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(theta_df) <- column_names
result_latentZ_mat = list()



# 초기값 설정
initial_beta_list <- list(
  c(0.5, 1, 2),
  c(0.5, 1, 25),
  # c(0.5, 1, 50),
  c(0.2, 1, 5),
  c(0.5, 1, 5),
  c(0.8, 1, 5)
)

initial_pi_list <- list(
  c(1, 1, 1),
  c(1, 2, 1),
  c(1, 10, 1)
)

# 모든 조합 생성
combination_list <- expand.grid(
  beta = seq_along(initial_beta_list),  # initial_beta의 인덱스
  pi = seq_along(initial_pi_list)       # initial_pi_set의 인덱스
)

# 각 조합을 출력
for (i in 1:nrow(combination_list)) {
  initial_beta <- initial_beta_list[[combination_list$beta[i]]]
  initial_pi_set <- initial_pi_list[[combination_list$pi[i]]]
  
  
  
initial_pi = initial_pi_set / sum(initial_pi_set)

initial_lambda = initial_lambda_calc(time_vec,event_vec,
                                     beta_vec=initial_beta,censored1 =5,censored3 = 30)


beta_vec = initial_beta
pi_vec = initial_pi
lambda_vec = initial_lambda
## initial latent variable 
latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=1)


objective_function1 <- function(new_beta1) {
  diffB_onlyB(new_beta1, latentZ_mat , j = 1)
}

objective_function3 <- function(new_beta1) {
  diffB_onlyB(new_beta3, latentZ_mat , j = 3)
}


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
  
  
  ###############
  print(" ### Beta1 Update ###")
  result1 = multiroot(
    f = function(new_beta1) diffB_onlyB(new_beta1, latentZ_mat , j = 1), 
    start = candi_before_vec[1])
  new_beta1 = result1$root
  
  print(" ### Beta3 Update ###")
  result3 = multiroot(
    f = function(new_beta3) diffB_onlyB(new_beta3, latentZ_mat , j = 3), 
    start = candi_before_vec[3])
  new_beta3 = result3$root
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
  alpha_temper = 1
  if(parameter_diff<1e-4){
    print("!!!!!!!!!!!!!!!!!!!! parameter diff Break !!!!!!!!!!!!!!!")
    theta_df_EMfull = rbind(theta_df_EMfull,
                          data.frame(
                            filename = paste0(unlist(strsplit(file_name, "\\."))[1]),
                            beta1=beta_vec[1],beta3=beta_vec[3],
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
                            Q3=Qfunc_onlyB(beta=beta_vec[3],latentZ_mat,j=3)
                          ))
    
    result_latentZ_mat[[iter]]=latentZ_mat
    break}
  
  #### E-Step ####
  latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=alpha_temper)
  
  
  # if(beta_vec[1]<0 || beta_vec[1]>2){break}
}

}
colnames(theta_df) <- column_names

theta_df %>% tail()

theta_df_EMfull %>% mutate( inital_Beta = paste0(init_beta1 ,"_",init_beta2,"_",init_beta3,"_",init_pi1,"_",init_pi2,"_",init_pi3)) %>% 
  group_by(inital_Beta) %>%       # Qlike 최대값 필터링
  slice(1) %>%     ungroup()   %>% as.data.frame()       


write.table(theta_df_EMfull,paste0(unlist(strsplit(file_name, "\\."))[1],"_EMInitResult.txt"),row.names = F)



tail(theta_df_EMfull)
plot(theta_df_EMfull$beta1)
plot(theta_df_EMfull$beta3)



LL = numeric()
for( i in 1:length(result_latentZ_mat)){
  LL[i] = result_latentZ_mat[[i]][1,1]/result_latentZ_mat[[i]][1,3]
}
result_latentZ_mat[[1]][1,1]/result_latentZ_mat[[1]][1,3]
LL







