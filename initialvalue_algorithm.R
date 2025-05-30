source("DAEM_BarrierMethod_function.R")

# Reading Data
# file_name = 'Aarest_data.txt'
# file_name = 'RearDump.txt'
# file_name = 'SerumReversal.txt'
file_name = 'FRT_censord.txt'
fdata = read.table(file_name,header = T)

# Data preprocessing
N = nrow(fdata)
k=3 
event_vec = fdata[,2] %>% as.numeric()
time_vec = fdata[,1]%>% as.numeric()

tot=1e-9
maxEMIter=1e+5
maxIterAnnealing = 1e+2
learningRateBp = 2
initialbpBase = 1e+2
# betaTot = 1e-4
initalAnnealingPara = 0.6
annealingSchedule = seq(0.66,0.9999,length.out=maxIterAnnealing)
bpBaseSchedule = seq(initialbpBase,1e+4,length.out=maxIterAnnealing)
betaTotSchedule = seq(1e+1,1e-8,length.out=maxIterAnnealing)



## result 기록 
result_latentZ_mat = list()
theta_df_full = NULL
theta_df = NULL


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
                                     initial_beta,censored1 = ceiling(N*0.1),censored3 = ceiling(N*0.4))
# initial_lambda = initial_lambda_calc(time_vec,event_vec,
#                                      beta_vec=initial_beta,censored1 =5,censored3 = 30)
## setting parameter
beta_vec = initial_beta
pi_vec = initial_pi
lambda_vec = initial_lambda
latentZ_mat = Estep_result(beta_vec,lambda_vec,pi_vec,alpha=1)


# latentZ_mat %>% round(4)
# print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))
# print(paste0(c("Beta1 at 3 is positive :",diffB_onlyB(1,latentZ_mat,j=3)>0,diffB_onlyB(1,latentZ_mat,j=3)) ,collapse = " / "))



for( ITerAnneal in 1:maxIterAnnealing){
  
  annealingPara = annealingSchedule[ITerAnneal]
  bpBase  = bpBaseSchedule[ITerAnneal]
  betaTot = betaTotSchedule[ITerAnneal]
  
  bp1=bpBase
  bp3=bpBase
  
  for( iter in 1:maxEMIter){ 
    
    #### M-Step ####
    #### Pi : M-Step :  #### 
    new_pi = colSums(latentZ_mat)/N
    #### Beta : M-Step : bijection method #### 
    candi_before_vec = beta_vec
    candi_lambda_vec = lambda_vec
    
    
    bp1=bpBase
    for( i in 1:10000){
      new_beta1 = barrier_beta1(candi_before_vec[1],latentZ_mat,bp=bp1 )
      # if(abs(barrierFunc_1(new_beta1,latentZ_mat,bp1))<tot){break}
      # if(abs(barrierFunc_1(new_beta1,latentZ_mat,bp1))<betaTot){break}
      # if(abs(diffB_onlyB(new_beta1,latentZ_mat,j=1))<betaTot){break}
      # if(abs(diffB_onlyB(new_beta1,latentZ_mat,j=1))<tot){break}
      bp1=bp1*learningRateBp
      if(bp1>1e+6){break}
    }
    bp3=bpBase
    for( i in 1:10000){
      new_beta3 = barrier_beta3(candi_before_vec[3],latentZ_mat,bp=bp3 )
      # if(abs(barrierFunc_3(new_beta3,latentZ_mat,bp3))<tot){break}
      # if(abs(barrierFunc_3(new_beta3,latentZ_mat,bp3))<betaTot){break}
      # if(abs(diffB_onlyB(new_beta3,latentZ_mat,j=3))<betaTot){break}
      # if(abs(diffB_onlyB(new_beta3,latentZ_mat,j=3))<tot){break}
      bp3=bp3*learningRateBp
      if(bp3>1e+6){break}
    }
    
    
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
    # if(parameter_diff<1e-10 || iter==maxEMIter){
    if(parameter_diff<tot || iter==maxEMIter){
      printResult()
      print(c("betaTot : ",betaTot))
      print("!!!!!!!!!!!!!!!!!!!! parameter diff Break !!!!!!!!!!!!!!!")
      theta_df_full = rbind(theta_df_full,
                            data.frame(filename = paste0(unlist(strsplit(file_name, "\\."))[1]),
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
    
    # print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))
    # print(paste0(c("Beta1 at 3 is positive :",diffB_onlyB(1,latentZ_mat,j=3)>0,diffB_onlyB(1,latentZ_mat,j=3)) ,collapse = " / "))
    
  }
}

}




# theta_df_full %>% tail
# plot(theta_df_full$Q1)
# plot(theta_df_full$Q2)
# plot(theta_df_full$Q3)
plot(theta_df_full$Qlike)
# plot(theta_df_full$Qlike[80:100])
# plot(theta_df_full$Q1+theta_df_full$Q3+theta_df_full$Q2)
plot(theta_df_full$beta1)
plot(theta_df_full$beta3)


theta_df_full %>% tail(35)
appAnnealLimit = theta_df_full$tempering[min(which(theta_df_full$Qlike > -589.7))]
theta_df_full %>% filter(tempering==appAnnealLimit)
qa = theta_df_full %>% ggplot(aes(x=tempering,y=Qlike))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(A.1) Likelihood Q("* beta*")"),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  annotate("text", x = appAnnealLimit, y = -690, label = appAnnealLimit, color = "red",
           hjust = 1,size=3)+
  theme_minimal()

q1 = theta_df_full %>% ggplot(aes(x=tempering,y=Q1))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(A.2) Likelihood "*Q[1]*"("* beta[1]*")"),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  # annotate("text", x = appAnnealLimit, y = -100, label = appAnnealLimit, color = "red", hjust = 1.0,size=3)+
  theme_minimal()

q2 = theta_df_full %>% ggplot(aes(x=tempering,y=Q2))+geom_line(lty=2)+geom_point()+
  labs(
    title =expression("(A.3) Likelihood "*Q[2]*"("* beta[2]*")"),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  # annotate("text", x = appAnnealLimit, y = -146, label = appAnnealLimit, color = "red", hjust = 1.0,size=3)+
  theme_minimal()

q3 = theta_df_full %>% ggplot(aes(x=tempering,y=Q3))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(A.4) Likelihood "*Q[3]*"("* beta[3]*")"),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  # annotate("text", x = appAnnealLimit, y = -21.6, label = appAnnealLimit, color = "red", hjust = 1.0,size=3)+
  theme_minimal()

beta1_gg = theta_df_full %>% ggplot(aes(x=tempering,y=beta1))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(B.1) Learning " * beta[1]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()

beta3_gg = theta_df_full %>% ggplot(aes(x=tempering,y=beta3))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(B.3) Learning " * beta[3]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()


beta3_ggDiff = theta_df_full %>% ggplot(aes(x=tempering,y=diffBeta3))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(B.4) Derivative of " * beta[3]),
    x = "Annealing parameter",
    y = ""
  ) +
  # geom_hline(yintercept = 0)+
  # coord_cartesian(ylim=c(-1e-5,0))+
  # scale_y_break(c(-1e-4, -1e-6), space = 0.3 , scales="free") +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()

beta1_ggDiff = theta_df_full %>% ggplot(aes(x=tempering,y=diffBeta1))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(B.2) Derivative of " * beta[1]),
    x = "Annealing parameter",
    y = ""
  ) +
  # scale_y_break(c(1e-4, 0.015), space = 0.3 , scales="free") +
  # coord_cartesian(ylim=c(0,1e-4))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()



Lambda1_gg = theta_df_full %>% ggplot(aes(x=tempering,y=lambda1))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(C.1) Learning " * lambda[1]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()

Lambda3_gg = theta_df_full %>% ggplot(aes(x=tempering,y=lambda3))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(C.2) Learning " * lambda[3]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()


pi_gg = theta_df_full %>% ggplot(aes(x=tempering,y=pi3))+
  geom_line(lty=2,size=1,color="black")+
  geom_line(aes(x=tempering,y=pi1),lty=3,size=1,color="black")+
  geom_line(aes(x=tempering,y=pi2),lty=4,size=1,color="black")+
  labs(
    title = expression("(C.3) Learning " * pi["k"]),
    x = "Annealing parameter",
    y = ""
  ) +
  # coord_cartesian(ylim=c(0.25,0.45))+
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  annotate("text", x = 0.95, y = 0.90, label = expression(pi[1]),size=5)+
  annotate("text", x = 0.95, y = 0.05, label = expression(pi[2]),size=5)+
  annotate("text", x = 0.95, y = 0.3, label = expression(pi[3]),size=5)+
  theme_minimal()
# pi_gg
# 

title_grob <- textGrob(
  paste0(unlist(strsplit(file_name, "\\."))[1]," Data"),
  gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
)

grid.arrange(
  qa,q1,q2,q3,
  beta1_gg,beta1_ggDiff,beta3_gg,beta3_ggDiff,
  Lambda1_gg,Lambda3_gg,pi_gg,
  ncol = 4,            # 열의 개수
  # heights = c(1, 1)     # 행 높이 비율
  top = title_grob
)

optimalData = theta_df_full %>% filter(tempering==appAnnealLimit) 
optimalData %>% pull(beta1)
optBeta1 = optimalData %>% pull(beta1)
optBeta2 = 1
optBeta3 = optimalData %>% pull(beta3)
optLambda1 = optimalData %>% pull(lambda1)
optLambda2 = optimalData %>% pull(lambda2)
optLambda3 = optimalData %>% pull(lambda3)

time_DBlist = seq(min(time_vec),max(time_vec),by=0.01)
DBbound1 = DecisionBoundary(time_DBlist,beta_vec = c(optimalData %>% pull(beta1),1,optimalData %>% pull(beta3)),lambda_vec = c(optimalData %>% pull(lambda1),optimalData %>% pull(lambda2),optimalData %>% pull(lambda3)),j=1)
DBbound3 = DecisionBoundary(time_DBlist,beta_vec = c(optimalData %>% pull(beta1),1,optimalData %>% pull(beta3)),lambda_vec = c(optimalData %>% pull(lambda1),optimalData %>% pull(lambda2),optimalData %>% pull(lambda3)),j=3)

changePoint1 = time_DBlist[which.min(DBbound1>1)]
changePoint3 = time_DBlist[which.min(DBbound3<1)]
DBdata = data.frame(time = time_DBlist, postProbRatio1 =DBbound1, postProbRatio3 =DBbound3 )
DB1 = DBdata %>% ggplot(aes(x=time,y=postProbRatio1))+geom_line()+  theme_minimal()+
  geom_hline(yintercept = 1,lty=2)+
  annotate("text", x = changePoint1+10, y = 1.3,label=paste0("Time :",changePoint1))+
  labs(
    title = "Posterior Ratio 1/2",
    x = "",
    y = ""
  )+ scale_y_continuous(
    breaks = c(1,seq(2, 7, by = 5)) # y축 눈금을 5 단위로 설정
  )+
  geom_point(data = data.frame(x = c(changePoint1), y = c(1)), aes(x = x, y = y), color = "blue", size = 3)
DB3 = DBdata %>% ggplot(aes(x=time,y=postProbRatio3))+geom_line()+ 
  labs(
    title = "Posterior Ratio 3/2",
    x = "",
    y = ""
  ) +
  annotate("text", x = changePoint3-10, y = 5,label=paste0("Time :",changePoint3))+
  geom_point(data = data.frame(x = c(changePoint3), y = c(1)), aes(x = x, y = y), color = "blue", size = 3)+
  theme_minimal()+geom_hline(yintercept = 1,lty=2)

posteriorPlot = DB1+DB3

hzData = data.frame(time = time_DBlist , hz1 = hazardrate(time_DBlist,optBeta1,optLambda1),
                    hz2 = hazardrate(time_DBlist,optBeta2,optLambda2),
                    hz3 = hazardrate(time_DBlist,optBeta3,optLambda3))

hzPlot = hzData %>% ggplot(aes(x=time , y=hz1))+geom_line(color="red")+
  geom_line(aes(x=time ,y=hz2),color="green3")+geom_line(aes(x=time ,y=hz3),color="blue")+
  # scale_y_break(c(0.3,1), space = 0.3 , scales="free") +
  # geom_vline(xintercept = 7.01)+
  # geom_vline(xintercept = 79.48)+
  labs(
    title = "Hazard Rate",
    x = "",
    y = ""
  )+
  # annotate("text", x = changePoint1+10, y = 1.3,label=paste0("Time :",changePoint1))+
  # annotate("text", x = changePoint3-10, y = 1.3,label=paste0("Time :",changePoint3))+
  theme_minimal()

DecisionPlot = posteriorPlot+hzPlot

DecisionPlot


abs_diff_sum <- sapply(2:length(result_latentZ_mat), function(i) {
  sum(abs(result_latentZ_mat[[i]] - result_latentZ_mat[[i-1]]))
})

data.frame(Annealing = annealingSchedule , latentDiff=c(abs_diff_sum[1],abs_diff_sum)) %>% 
  ggplot(aes(x=Annealing,y=latentDiff))+geom_line()


dim(theta_df_full)

theta_df_full %>% mutate( inital_Beta = paste0(init_beta1 ,"_",init_beta2,"_",init_beta3,"_",init_pi1,"_",init_pi2,"_",init_pi3)) %>% 
  group_by(inital_Beta) %>%   filter(Qlike == max(Qlike)) %>%     # Qlike 최대값 필터링
  slice(1) %>%     ungroup()   %>% as.data.frame()       
unique(theta_df_full$inital_Beta)

theta_df_full %>% mutate( inital_Beta = paste0(init_beta1 ,"_",init_beta2,"_",init_beta3,"_",init_pi1,"_",init_pi2,"_",init_pi3)) %>% pull(inital_Beta) %>% 
  unique()
  

write.table(theta_df_full,paste0(unlist(strsplit(file_name, "\\."))[1],"_InitResult.txt"),row.names = F)




