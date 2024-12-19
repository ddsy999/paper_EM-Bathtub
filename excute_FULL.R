source("DAEM_BarrierMethod_function.R")






# Reading Data
fdata = read.table('Aarest_data.txt',header = T)

# Data preprocessing
N = nrow(fdata)
k=3 
event_vec = fdata[,2] %>% as.numeric()
time_vec = fdata[,1]%>% as.numeric()
# time_vec = time_vec/(max(time_vec)*1.1)
tot=1e-8
maxBp = 1e+5
learning_rate = 1.01

# 11개의 열을 가진 빈 result data.frame 생성
column_names <- c("beta1", "lambda1", "beta2", "lambda2", "beta3", "lambda3", "sumQfunc","diffB_beta1","diffB_beta3","bp","iter","Beta1 at 1","Q1","Q2","Q3","pi1","pi2","pi3")
theta_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(theta_df) <- column_names

latent_1 = data.frame(V1 = numeric(), V2 = numeric(), V3 = numeric(), V4=numeric())
latent_2 = data.frame(V1 = numeric(), V2 = numeric(), V3 = numeric(), V4=numeric())
latent_3 = data.frame(V1 = numeric(), V2 = numeric(), V3 = numeric(), V4=numeric())
result_latentZ_mat = list()
# TODO : delete latent_1 ,latent_2,latent_3



## initial Parameter : beta , lambda , pi
initial_beta = c(0.5,1,5)
initial_pi_set = c(1,1,1)
initial_pi = initial_pi_set / sum(initial_pi_set)
initial_lambda = c(3,2,1.5)
beta_vec = initial_beta
pi_vec = initial_pi
lambda_vec = initial_lambda
## initial latent variable 
latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec)
latentZ_mat



for( iter in 1:100){ 
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
  #### Beta : M-Step : G-D method #### 
  candi_before_vec = beta_vec
  candi_lambda_vec = lambda_vec
  

barrier_beta1(candi_before_vec[1],latentZ_mat,bp=1e1 )

abs(diffB_onlyB(1,latentZ_mat,j=1))

  ###############
  print(" ### Beta1 Update ###")
  bp1=1e+50*(iter+1)
  new_beta1 = barrier_beta1(candi_before_vec[1],latentZ_mat,bp=bp1 )
  
  
  print(" ### Beta3 Update ###")
  bp3=1e+50*(iter+1)
  for( i in 1:10000){
    new_beta3=barrier_beta3_NR(candi_before_vec[3],latentZ_mat,j=3,bp3) 
    if(abs(diffB_onlyB(new_beta3,latentZ_mat,j=3))<tot){
      print(paste0("@@@ Beta3 diff limit complete ",diffB_onlyB(new_beta3,latentZ_mat,j=3)))
      break}
    bp3=bp3*1.01
    candi_before_vec[3]=new_beta3
  }
  
  ################
  
  
  new_beta = c(new_beta1,1,new_beta3)
  new_lambda = sapply(1:k , function(i)  sum(latentZ_mat[,i]*event_vec)/sum(latentZ_mat[,i]*(time_vec^new_beta[i])))
  
  
  print(c("differential : ", diffB_onlyB(new_beta[1],latentZ_mat,j=1)," ",diffB_onlyB(new_beta[3],latentZ_mat,j=3)))
  
  
  
  
  beta1Up = 0 
  beta3Up = 0
  ### Beta Update ### 
  if(Qfunc_onlyB(beta=beta_vec[1],latentZ_mat,j=1)<Qfunc_onlyB(beta=new_beta[1],latentZ_mat,j=1)){
    print(" Beta 1 Update")
    beta1Up=1
  }
  if(Qfunc_onlyB(beta=beta_vec[3],latentZ_mat,j=3)<Qfunc_onlyB(beta=new_beta[3],latentZ_mat,j=3)){
    print(" Beta 3 Update")
    beta3Up=1
  }
  
  
  
  parameter_diff = sqrt(sum((beta_vec-new_beta)^2+(pi_vec-new_pi)^2))
  
  ### Pi Update ###
  pi_vec = new_pi
  
  ### Beta Lambda Update ### 
  if(beta1Up>0){ 
    beta_vec[1]=new_beta[1]
    lambda_vec[1]=new_lambda[1]
  }
  if(beta3Up>0){ 
    beta_vec[3]=new_beta[3]
    lambda_vec[3]=new_lambda[3]
  }
  beta_vec[2]=1
  lambda_vec[2] = new_lambda[2]
  
  Update_start=1
  if(Update_start == 1){
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
  }
  
  candi_latent_1 = data.frame(V1=sum(latentZ_mat[,1]*event_vec)/beta_vec[1],
                              V2=sum(latentZ_mat[,1]*event_vec*log(time_vec)),
                              V3=-sum(latentZ_mat[,1]*event_vec)*sum(latentZ_mat[,1]*(time_vec^beta_vec[1])*log(time_vec))/sum(latentZ_mat[,1]*(time_vec^beta_vec[1]))
                              )
  latent_1 = rbind(latent_1,candi_latent_1)
  
  candi_latent_3 = data.frame(V1=sum(latentZ_mat[,3]*event_vec)/beta_vec[3],
                              V2=sum(latentZ_mat[,3]*event_vec*log(time_vec)),
                              V3=-sum(latentZ_mat[,3]*event_vec)*sum(latentZ_mat[,3]*(time_vec^beta_vec[3])*log(time_vec))/sum(latentZ_mat[,3]*(time_vec^beta_vec[3]))
  )
  latent_3 = rbind(latent_3,candi_latent_3)
  
  
  result_latentZ_mat[[iter]]=latentZ_mat
  
  #### Stopping rule ####
  alpha_temper = 0.8
  if(parameter_diff<1e-7){
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
    
    
    break}
  
  #### E-Step ####
  latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=alpha_temper)
  
}


write.csv(theta_df_full,"init_para.csv",row.names = F)









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
)



colnames(theta_df) <- column_names
theta_df$beta1 = as.numeric(theta_df$beta1)
# plot(theta_df[,'sumQfunc'])
tail(theta_df)
# theta_df
p1=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=sumQfunc %>% as.numeric()))+geom_point()+geom_line()+ggtitle("sumQfunc")
p2=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=beta1 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Beta1")
p3=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=beta3 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Beta3")

p4=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=diffB_beta1 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("diffB_beta1")
p5=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=diffB_beta3 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("diffB_beta3")

p6=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=Q1 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Q1")
p7=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=Q2 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Q2")
p8=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=Q3 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Q3")
# 두 개의 그림을 결합하여 하나로 표시
combined_plot <- p1+p2+p3+p4+p5+p6+p7+p8

# 결합된 그림을 출력
print(combined_plot)


aa=latent_1 %>% mutate(total=V1+V2+V3)
bb=latent_3 %>% mutate(total=V1+V2+V3)
# plot(theta_df[,"sumQfunc"])

plot(aa[1:27,"total"])
plot(aa[,"total"])
# theta_df[26,]

plot(aa[,"V1"])
plot(aa[,"V2"])
plot(aa[,"V3"])
plot(aa[,"total"])





k1=30
A1=weibull_func(time_vec,theta_df$beta1[k1] %>% as.numeric(),theta_df$lambda1[k1] %>% as.numeric()) 
A2=weibull_func(time_vec,theta_df$beta2[k1] %>% as.numeric(),theta_df$lambda2[k1] %>% as.numeric()) 
A3=weibull_func(time_vec,theta_df$beta3[k1] %>% as.numeric(),theta_df$lambda3[k1] %>% as.numeric()) 

k2=15
B1=weibull_func(time_vec,theta_df$beta1[k2] %>% as.numeric(),theta_df$lambda1[k2] %>% as.numeric()) 
B2=weibull_func(time_vec,theta_df$beta2[k2] %>% as.numeric(),theta_df$lambda2[k2] %>% as.numeric()) 
B3=weibull_func(time_vec,theta_df$beta3[k2] %>% as.numeric(),theta_df$lambda3[k2] %>% as.numeric()) 


ppd=data.frame(time = 1:length(A1),A1,A2,A3,B1,B2,B3)
ppd %>% ggplot(aes(x=time,y=A1 %>% log))+geom_line(color="red")+
  geom_line(data=ppd,aes(x=time,y=B1 %>% log),color="blue")

ppd %>% ggplot(aes(x=time,y=A1 ))+geom_line(color="red")+
  geom_line(data=ppd,aes(x=time,y=B1 ),color="blue")

ppd %>% ggplot(aes(x=time,y=A2 %>% log))+geom_line(color="red")+
  geom_line(data=ppd,aes(x=time,y=B2 %>% log),color="blue")

ppd %>% ggplot(aes(x=time,y=A3 %>% log))+geom_line(color="red")+
  geom_line(data=ppd,aes(x=time,y=B3 %>% log),color="blue")