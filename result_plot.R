
# roadResult = read.table("Proposed_Result_RearDump.txt")
# roadResult = read.table("Proposed_Result_FRT_censord.txt")
roadResult = read.table("Proposed_Result_Aarest_data.txt")
# roadResult = read.table("Proposed_Result_SerumReversal.txt")

roadResult = read.table("Proposed_Result_NoDAEM_Aarest_data.txt")
# roadResult = read.table("Proposed_Result_NoDAEM_FRT_censord.txt")
# roadResult = read.table("Proposed_Result_NoDAEM_RearDump.txt")
# roadResult = read.table("Proposed_Result_NoDAEM_SerumReversal.txt")
optFilename = roadResult$data_Name %>% unique()



roadResult %>% ggplot(aes(x=Iter , y=Qlike))+geom_point()+geom_line()+labs(title="Expected log-likelihood",y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))

roadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))+geom_hline(yintercept = 1,lty=2,color="red")
roadResult %>% ggplot(aes(x=Iter , y=1))    +geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))+geom_hline(yintercept = 1,lty=2,color="red")
roadResult %>% ggplot(aes(x=Iter , y=lambda1  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=lambda2  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=lambda3  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))




A1=roadResult %>% ggplot(aes(x=Iter , y=Qlike))+geom_point()+geom_line()+labs(title="Expected log-likelihood",y="",x="Iteration")+theme(axis.text = element_text(size=12))
A2=roadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A3=roadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A4=roadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))+geom_hline(yintercept = 1,lty=2,color="red")
A5=roadResult %>% ggplot(aes(x=Iter , y=1))    +geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A6=roadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))+geom_hline(yintercept = 1,lty=2,color="red")
A7=roadResult %>% ggplot(aes(x=Iter , y=lambda1  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A8=roadResult %>% ggplot(aes(x=Iter , y=lambda2  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A9=roadResult %>% ggplot(aes(x=Iter , y=lambda3  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A10=roadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A11=roadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A12=roadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))




title_grob <- textGrob(
  optFilename, 
  gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
)

grid.arrange(
  A1,A2,A3,
  A4,A5,A6,
  A7,A8,A9,
  A10,A11,A12,
  ncol = 3,           # 열의 개수
  # heights = c(1, 1)     # 행 높이 비율
  top = title_grob
)


###############################################


optimalData = roadResult %>% filter(tempering==0.999999999) 
optFilename = roadResult$data_Name %>% unique()
optBeta1 = optimalData %>% pull(beta1)
optBeta2 = 1
optBeta3 = optimalData %>% pull(beta3)
optLambda1 = optimalData %>% pull(lambda1)
optLambda2 = optimalData %>% pull(lambda2)
optLambda3 = optimalData %>% pull(lambda3)
optPi1 = optimalData %>% pull(pi1)
optPi2 = optimalData %>% pull(pi2)
optPi3 = optimalData %>% pull(pi3)

time_DBlist = seq(min(time_vec),max(time_vec),by=0.1)

DBbound1 = DecisionBoundary(time_DBlist,beta_vec = c(optBeta1,1,optBeta3),lambda_vec = c(optLambda1,optLambda2,optLambda3),pi_vec=c(optPi1,optPi2,optPi3),j=1)
DBbound3 = DecisionBoundary(time_DBlist,beta_vec = c(optBeta1,1,optBeta3),lambda_vec = c(optLambda1,optLambda2,optLambda3),pi_vec=c(optPi1,optPi2,optPi3),j=3)


changePoint1 = time_DBlist[which.min(DBbound1>1)]
changePoint3 = time_DBlist[which.min(DBbound3<1)]
DBdata = data.frame(time = time_DBlist, postProbRatio1 =DBbound1, postProbRatio3 = DBbound3 )
DB1 = DBdata %>% ggplot(aes(x=time,y=postProbRatio1))+geom_line()+  theme_minimal()+
  geom_hline(yintercept = 1,lty=2)+
  annotate("text", x = changePoint1+10, y = 1.3,label=paste0("Time :",changePoint1))+
  labs(
    title = "Posterior Ratio (infant vs constant)",
    x = "",
    y = ""
  )+ scale_y_continuous(
    breaks = c(1,seq(2, 7, by = 5)) # y축 눈금을 5 단위로 설정
  )+
  geom_point(data = data.frame(x = c(changePoint1), y = c(1)), aes(x = x, y = y), color = "blue", size = 3)
DB3 = DBdata %>% ggplot(aes(x=time,y=postProbRatio3))+geom_line()+ 
  labs(
    title = "Posterior Ratio (wear-out vs constant)",
    x = "",
    y = ""
  ) +
  annotate("text", x = changePoint3-10, y = 3,label=paste0("Time :",changePoint3))+
  geom_point(data = data.frame(x = c(changePoint3), y = c(1)), aes(x = x, y = y), color = "blue", size = 3)+
  theme_minimal()+geom_hline(yintercept = 1,lty=2)

posteriorPlot = DB1+DB3

hzData = data.frame(time = time_DBlist , hz1 = hazardrate(time_DBlist,optBeta1,optLambda1),
                    hz2 = hazardrate(time_DBlist,optBeta2,optLambda2),
                    hz3 = hazardrate(time_DBlist,optBeta3,optLambda3))

hzPlot = hzData %>% ggplot(aes(x=time , y=hz1))+geom_line(color="red")+
  geom_line(aes(x=time ,y=hz2),color="green3")+geom_line(aes(x=time ,y=hz3),color="blue")+
  # scale_y_break(c(0.3,2), space = 0.3 , scales="free") +
  geom_vline(xintercept = changePoint1,lty=2)+
  geom_vline(xintercept = changePoint3,lty=2)+
  labs(
    title = "Hazard Rate",
    x = "",
    y = ""
  )+
  annotate("text", x = changePoint1, y = 0.01,hjust=-0.2,label=paste0("Time :",changePoint1))+
  annotate("text", x = changePoint3, y = 0.01,hjust=-0.3,label=paste0("Time :",changePoint3))+
  theme_minimal()

DecisionPlot = posteriorPlot+hzPlot

grid.arrange(
  DB1,DB3,
  ncol=2,
  top = textGrob(
    "Posterior Ratio in the Phase of Device G Failure", 
    gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
  )
)


########################################################

overtime_DBlist = seq(min(time_vec),max(time_vec)*3,by=0.1)
OverTime1 = DecisionBoundary(overtime_DBlist,beta_vec = c(optBeta1,1,optBeta3),lambda_vec = c(optLambda1,optLambda2,optLambda3),pi_vec=c(optPi1,optPi2,optPi3),j=1)
OverTime3 = DecisionBoundary(overtime_DBlist,beta_vec = c(optBeta1,1,optBeta3),lambda_vec = c(optLambda1,optLambda2,optLambda3),pi_vec=c(optPi1,optPi2,optPi3),j=3)

OverTime = data.frame( time = overtime_DBlist , posterior1=OverTime1,posterior3=OverTime3 )

OverTime %>% ggplot(aes(x=time, y=posterior1))+
  geom_line()+geom_hline(yintercept = 1,lty=2)


OverTime %>% ggplot(aes(x=time, y=posterior1)) +
  geom_rect(aes(xmin=0, xmax=300, ymin=-Inf, ymax=Inf), fill="lightgray", alpha=0.3) + 
  geom_line() + 
  geom_hline(yintercept = 1, lty=2) +
  annotate("text", x=150, y=max(OverTime$posterior1, na.rm=TRUE)*1.01, 
           label="Observed Data Range (0-300)", vjust=-0.5, size=4)+
  # annotate("text", x = changePoint1+10, y = 1.3,label=paste0("Time :",changePoint1))+
  labs(
    title = "",
    x = "",
    y = ""
  )+
  geom_vline(xintercept = 300,lty=2)+
  theme(plot.title = element_text(hjust = 0.5))

OverTime %>% ggplot(aes(x=time, y=posterior3)) +
  geom_rect(aes(xmin=0, xmax=300, ymin=-Inf, ymax=Inf), fill="lightgray", alpha=0.3) + 
  geom_line() + 
  geom_hline(yintercept = 1, lty=2) +
  annotate("text", x=150, y=max(OverTime$posterior3, na.rm=TRUE)*1.01, 
           label="Observed Data Range (0-300)", vjust=-0.5, size=4)+
  # annotate("text", x = changePoint1+10, y = 1.3,label=paste0("Time :",changePoint1))+
  labs(
    title = "",
    x = "",
    y = ""
  )+
  geom_vline(xintercept = 300,lty=2)+
  theme(plot.title = element_text(hjust = 0.5))

