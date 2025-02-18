source("DAEM_BarrierMethod_function.R")

# Reading Data
# file_name = 'Aarest_data.txt'
# file_name = 'RearDump.txt'
# file_name = 'SerumReversal.txt'
# file_name = 'FRT_censord.txt'
file_name = 'LFP.txt'
fdata = read.table(file_name,header = T)

# Data preprocessing
N = nrow(fdata)
k=2
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
# bpBaseSchedule = exp(seq(log(1e+3), log(1e+5), length.out = maxIterAnnealing))
bpBaseSchedule = seq(1e+2,1e+5,length.out=maxIterAnnealing)
betaTotSchedule = seq(1e+1,1e-8,length.out=maxIterAnnealing)


## result 기록 
result_latentZ_mat = list()
theta_df_full = NULL
theta_df = NULL
## initial Parameter : beta , lambda , pi
# initial_beta = c(0.5,1,2)
initial_beta = c(0.5,1)
# initial_pi_set = c(1,1,1)
initial_pi_set = c(1,1)
initial_pi = initial_pi_set / sum(initial_pi_set)
# initial_lambda = c(1e-5,1e-10)
initial_lambda = c(1,1e-5)

## setting parameter
beta_vec = initial_beta
pi_vec = initial_pi
lambda_vec = initial_lambda
latentZ_mat = Estep_result2(beta_vec,lambda_vec,pi_vec,alpha=1)

print(paste0(c("Beta1 at 1 is minus :",diffB_onlyB(1,latentZ_mat,j=1)<0,diffB_onlyB(1,latentZ_mat,j=1)) ,collapse = " / "))

diffB_onlyB(0,latentZ_mat,j=1)

for( ITerAnneal in 1:maxIterAnnealing){
  
  annealingPara = annealingSchedule[ITerAnneal]
  bpBase  = bpBaseSchedule[ITerAnneal]
  betaTot = betaTotSchedule[ITerAnneal]
  
  bp1=bpBase
  # bp3=bpBase
  
  for( iter in 1:maxEMIter){ 
    
    #### M-Step ####
    #### Pi : M-Step :  #### 
    new_pi = colSums(latentZ_mat)/N
    #### Beta : M-Step : bijection method #### 
    candi_before_vec = beta_vec
    candi_lambda_vec = lambda_vec

    
    # newton_raphson(beta_init = candi_before_vec[1], lambda = candi_lambda_vec[1] ,latentZ_mat, j =1 )
    new_beta1 = barrier_beta1(candi_before_vec[1],latentZ_mat,bp=bpBase)
    # barrierFunc_1
    # diffB_onlyB
    #### Update Parameter ####
    new_beta = c(new_beta1,1)
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
      printResult2()
      print(c("betaTot : ",betaTot))
      print("!!!!!!!!!!!!!!!!!!!! parameter diff Break !!!!!!!!!!!!!!!")
      theta_df_full = rbind(theta_df_full,
                            data.frame(beta1=beta_vec[1],
                                       lambda1=lambda_vec[1],lambda2=lambda_vec[2],
                                       diffBeta1=diffB_onlyB(beta_vec[1],latentZ_mat,j=1),
                                       diffBeta2=diffB_onlyB(beta_vec[2],latentZ_mat,j=2),
                                       # diffBeta3=diffB_onlyB(beta_vec[3],latentZ_mat,j=3),
                                       diffLambda1 = diffL(beta_vec[1],lambda_vec[1],latentZ_mat,j=1),
                                       # diffLambda3 = diffL(beta_vec[3],lambda_vec[3],latentZ_mat,j=3),
                                       pi1=pi_vec[1],
                                       pi2=pi_vec[2],
                                       # pi3=pi_vec[3],
                                       init_beta1 = initial_beta[1],
                                       init_beta2 = initial_beta[2],
                                       # init_beta3 = initial_beta[3],
                                       init_pi1   = initial_pi[1],
                                       init_pi2   = initial_pi[2],
                                       # init_pi3   = initial_pi[3],
                                       tempering = alpha_temper,
                                       Qlike = Qfunc(beta=beta_vec[1],lambda_vec[1],latentZ_mat,j=1)+Qfunc(beta=beta_vec[2],lambda_vec[2],latentZ_mat,j=2),
                                       Q1=Qfunc_onlyB(beta=beta_vec[1],latentZ_mat,j=1),
                                       Q2=Qfunc_onlyB(beta=beta_vec[2],latentZ_mat,j=2),
                                       # Q3=Qfunc_onlyB(beta=beta_vec[3],latentZ_mat,j=3),
                                       bpBase = bpBase,
                                       bp1=bp1
                                       # bp3=bp3
                            ))
      
      result_latentZ_mat[[ITerAnneal]]=latentZ_mat
      break}
    
    #### E-Step ####
    latentZ_mat = Estep_result2(beta_vec,lambda_vec,pi_vec,alpha=alpha_temper)
  }
}


# theta_df_full$beta3
theta_df_full = theta_df_full %>% mutate(across(-1, as.numeric))

plot(theta_df_full$Qlike)
# plot(theta_df_full$Qlike[80:100])
# plot(theta_df_full$Q1+theta_df_full$Q3+theta_df_full$Q2)
# plot(theta_df_full$beta1)
# plot(theta_df_full$beta3)
# plot(theta_df_full$diffBeta1)
# plot(theta_df_full$diffBeta3)
# plot(theta_df_full$bpBase)

theta_df_full %>% tail

appAnnealLimit = theta_df_full$tempering[max(which(theta_df_full$Qlike < 3))]
optimalData = theta_df_full %>% filter(tempering==appAnnealLimit)
OptPi1 = optimalData$pi1
OptPi2 = optimalData$pi2
# OptPi3 = optimalData$pi3
OptBeta1 = optimalData$beta1
# OptBeta3 = optimalData$beta3
OptLike = optimalData$Qlike
OptDiffBeta1 =  signif(optimalData$diffBeta1, 3) 
# OptDiffBeta3 =  signif(optimalData$diffBeta3, 3) 
qa = theta_df_full %>% ggplot(aes(x=tempering,y=Qlike))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(A) Likelihood Q("* beta*")"),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_point(aes(x=appAnnealLimit,y=OptLike),color="red",size=3)+
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  # annotate("text", x = appAnnealLimit, y = -245, label = appAnnealLimit, color = "red", 
  #          hjust = -0.1,size=3)+
  theme_minimal()


beta1_gg = theta_df_full %>% ggplot(aes(x=tempering,y=beta1))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(B) Convergence of " * beta[1]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_point(aes(x=appAnnealLimit,y=OptBeta1),color="red",size=3)+
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()


# beta3_gg = theta_df_full %>% ggplot(aes(x=tempering,y=beta3))+geom_line(lty=2)+geom_point()+
#   labs(
#     title = expression("(C) Convergence of " * beta[3]),
#     x = "Annealing parameter",
#     y = ""
#   ) +
#   geom_point(aes(x=appAnnealLimit,y=OptBeta3),color="red",size=3)+
#   geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
#   theme_minimal()


# beta3_ggDiff = theta_df_full %>% ggplot(aes(x=tempering,y=diffBeta3))+geom_line(lty=2)+geom_point()+
#   labs(
#     title = expression("(D) Derivative of " * beta[3]),
#     x = "Annealing parameter",
#     y = ""
#   ) +
#   geom_vline(xintercept = appAnnealLimit,color="red",linetype=3,size=1)+
#   annotate("text", x = appAnnealLimit-0.05, y = OptDiffBeta3-1e-5, label =OptDiffBeta3, size = 5) +
#   theme_minimal()+
#   theme(axis.text.y=element_text(size=10))

beta1_ggDiff = theta_df_full %>% ggplot(aes(x=tempering,y=diffBeta1))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(E) Derivative of " * beta[1]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  annotate("text", x = appAnnealLimit+0.05, y = OptDiffBeta1-0.01, label =OptDiffBeta1, size = 5) +
  theme_minimal()

theta_df_full %>% ggplot(aes(x=tempering,y=diffBeta2))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(E) Derivative of " * beta[1]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  annotate("text", x = appAnnealLimit+0.05, y = OptDiffBeta1-0.01, label =OptDiffBeta1, size = 5) +
  theme_minimal()

pi_gg = theta_df_full %>%
  ggplot(aes(x = tempering)) +
  geom_line(aes(y = pi1, linetype = "pi1"), size = 1, color = "black") +
  geom_line(aes(y = pi2, linetype = "pi2"), size = 1, color = "black") +
  # geom_line(aes(y = pi3, linetype = "pi3"), size = 1, color = "black") +
  scale_linetype_manual(
    values = c("pi1" = 3, "pi2" = 4),
    labels = c(expression(pi[1]), expression(pi[2])),
    name = "" ) +
  labs(
    title = expression("(F) Convergence of " * pi["k"]),
    x = "Annealing parameter",
    y = ""
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_vline(xintercept = appAnnealLimit, color = "red", linetype = 3, size = 1) +
  theme_minimal()+
  theme(
    legend.text = element_text(size = 15),         # 🔥 범례 텍스트 크기 조정
    legend.title = element_text(size = 10),        # 🔥 범례 제목 크기 조정
    legend.key.size = unit(1, "cm")              # 🔥 범례 키(선 스타일 크기) 조정
  )




title_grob <- textGrob(
  "Limited ", 
  gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
)

parameterPlot = grid.arrange(
  qa,beta1_gg,beta1_ggDiff,pi_gg,
  ncol = 4,            # 열의 개수
  # heights = c(1, 1)     # 행 높이 비율
  top = title_grob
)


theta_df_full %>% filter(tempering>=appAnnealLimit) 
optimalData %>% pull(beta1)
optBeta1 = optimalData %>% pull(beta1)
optBeta2 = 1
# optBeta3 = optimalData %>% pull(beta3)
optLambda1 = optimalData %>% pull(lambda1)
optLambda2 = optimalData %>% pull(lambda2)
# optLambda3 = optimalData %>% pull(lambda3)

time_DBlist = seq(min(time_vec),max(time_vec),by=0.1)
DBbound1 = DecisionBoundary(time_DBlist,beta_vec = c(optimalData %>% pull(beta1),1,0),lambda_vec = c(optimalData %>% pull(lambda1),optimalData %>% pull(lambda2),0),j=1)
# DBbound3 = DecisionBoundary(time_DBlist,beta_vec = c(optimalData %>% pull(beta1),1,optimalData %>% pull(beta3)),lambda_vec = c(optimalData %>% pull(lambda1),optimalData %>% pull(lambda2),optimalData %>% pull(lambda3)),j=3)

# pi_vec = c(pi_vec[2],pi_vec[1])
changePoint1 = time_DBlist[which.min(DBbound1>1)]
# changePoint3 = time_DBlist[which.min(DBbound3<1)]
DBdata = data.frame(time = time_DBlist, postProbRatio1 =DBbound1 )
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
# DB3 = DBdata %>% ggplot(aes(x=time,y=postProbRatio3))+geom_line()+ 
#   labs(
#     title = "Posterior Ratio 3/2",
#     x = "",
#     y = ""
#   ) +
#   annotate("text", x = changePoint3-10, y = 5,label=paste0("Time :",changePoint3))+
#   geom_point(data = data.frame(x = c(changePoint3), y = c(1)), aes(x = x, y = y), color = "blue", size = 3)+
#   theme_minimal()+geom_hline(yintercept = 1,lty=2)

posteriorPlot = DB1

hzData = data.frame(time = time_DBlist , hz1 = hazardrate(time_DBlist,optBeta1,optLambda1),
                    hz2 = hazardrate(time_DBlist,optBeta2,optLambda2))

which((pi_vec[1]/pi_vec[2])*(hzData$hz1/hzData$hz2)*exp(-lambda_vec[1]*(time_DBlist^beta_vec[1]) + lambda_vec[2]*time_DBlist)<1)

hzPlot = hzData %>% ggplot(aes(x=time , y=hz1))+geom_line(color="red")+
  geom_line(aes(x=time ,y=hz2),color="green3")+
  # scale_y_break(c(0.3,2), space = 0.3 , scales="free") +
  # geom_vline(xintercept = 7.01,lty=2)+
  # geom_vline(xintercept = 79.48,lty=2)+
  labs(
    title = "Hazard Rate",
    x = "",
    y = ""
  )+
  # coord_cartesian(ylim = c(0,0.3))+
  # annotate("text", x = changePoint1, y = 0.3,hjust=-0.2,label=paste0("Time :",changePoint1))+
  # annotate("text", x = changePoint3, y = 0.3,hjust=1.1,label=paste0("Time :",changePoint3))+
  theme_minimal()

DecisionPlot = posteriorPlot+hzPlot

grid.arrange(
  DB1,DB3,hzPlot,
  ncol=3,
  top = textGrob(
    "Aarest Decision Boundary", 
    gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
  )
)


write.table(theta_df_full,paste0(paste0(unlist(strsplit(file_name, "\\."))[1]),"_DAEM.txt"),row.names = F)


