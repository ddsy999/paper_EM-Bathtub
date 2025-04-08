

# roadResult = read.table("Proposed_Result_RearDump.txt")
# roadResult = read.table("Proposed_Result_Aarest_data.txt")
# roadResult = read.table("Proposed_Result_FRT_censord.txt")
roadResult = read.table("Proposed_Result_SerumReversal.txt")

# roadResult = read.table("Proposed_Result_NoDAEM_Aarest_data.txt")
# roadResult = read.table("Proposed_Result_NoDAEM_FRT_censord.txt")
# roadResult = read.table("Proposed_Result_NoDAEM_RearDump.txt")
# roadResult = read.table("Proposed_Result_NoDAEM_SerumReversal.txt")
optFilename = roadResult$data_Name %>% unique()
optimalEnd = which.min(abs(roadResult$diffBeta1)+abs(roadResult$diffBeta3))
# optimalEnd = 99
optimalData = roadResult[optimalEnd ,]


roadResult %>% ggplot(aes(x=Iter , y=Qlike))+geom_point()+geom_line()+labs(title="Expected log-likelihood",y="",x="Iteration")+theme(axis.text = element_text(size=8))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# roadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# 
# roadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=1))    +geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[2]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=lambda1  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[1]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=lambda2  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[2]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=lambda3  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[3]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="Iteration")+theme(axis.text = element_text(size=8))
# roadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="Iteration")+theme(axis.text = element_text(size=8))

# roadResult %>% tail
roadResult = roadResult[50:100,]

# A1=roadResult %>% ggplot(aes(x=Iter , y=Qlike))+geom_point()+geom_line()+labs(title="Expected log-likelihood",y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A2=roadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10) )+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A3=roadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A4=roadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A5=roadResult %>% ggplot(aes(x=Iter , y=1))    +geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[2]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A6=roadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A7=roadResult %>% ggplot(aes(x=Iter , y=lambda1  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[1]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A8=roadResult %>% ggplot(aes(x=Iter , y=lambda2  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[2]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A9=roadResult %>% ggplot(aes(x=Iter , y=lambda3  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[3]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)+scale_y_continuous(labels = label_scientific(digits = 1))
# A10=roadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A11=roadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
# A12=roadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="Iteration")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 10))+geom_vline(xintercept = optimalEnd,color="red",lty=2)

A1=roadResult %>% ggplot(aes(x=Iter , y=Qlike))+geom_point()+geom_line()+labs(title="Expected log-likelihood",y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A2=roadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12) )+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A3=roadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A4=roadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A5=roadResult %>% ggplot(aes(x=Iter , y=1))    +geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[2]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A6=roadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A7=roadResult %>% ggplot(aes(x=Iter , y=lambda1  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[1]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A8=roadResult %>% ggplot(aes(x=Iter , y=lambda2  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[2]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A9=roadResult %>% ggplot(aes(x=Iter , y=lambda3  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[3]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)+scale_y_continuous(labels = label_scientific(digits = 1))
A10=roadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A11=roadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)
A12=roadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="")+theme(axis.text = element_text(size=8), plot.title = element_text(size = 12))+geom_vline(xintercept = optimalEnd,color="red",lty=2)



library(tidyr)
# 필요한 열만 추출
pi_data <- roadResult %>%
  select(Iter, pi1, pi2, pi3) %>%
  pivot_longer(cols = starts_with("pi"), names_to = "Component", values_to = "Value")

# 수렴 그래프 그리기
p_pi <- ggplot(pi_data, aes(x = Iter, y = Value, color = Component)) +
  geom_line(size = 0.5) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = optimalEnd,color="red",lty=2)+
  labs(
    title = "Convergence of π",
    x = "Iteration",
    y = expression(pi)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top"
  )

# print(p_pi)

# 필요한 열만 추출 및 정리
grad_data <- roadResult %>%
  select(Iter, diffBeta1, diffBeta3) %>%
  pivot_longer(
    cols = -Iter,
    names_to = "Component",
    values_to = "Value"
  )
roadResult$diffLambda3

# 이름 보기 좋게 변경 (선택 사항)
grad_data$Component <- recode(grad_data$Component,
                              diffBeta1 = "∇β₁",
                              diffBeta3 = "∇β₃",
                              diffLambda1 = "∇λ₁",
                              diffLambda3 = "∇λ₃")


# 변환된 값 추가
grad_data_log <- grad_data %>%
  mutate(LogValue = log1p(abs(Value)))  # log1p(x) = log(1 + x)

# 시각화
p_grad_log <- ggplot(grad_data_log, aes(x = Iter, y = LogValue, color = Component)) +
  geom_line(size = 0.5) +
  geom_vline(xintercept = optimalEnd,color="red",lty=2)+
  geom_point(size = 1.5) +
  labs(
    title = "Gradient Convergence", #log(1 + |∇|) t
    x = "Iteration",
    y = expression(log(1 + abs(nabla)))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top"
  )



grid.arrange(
  A4, A6, p_grad_log+facet_wrap(.~Component, scales = "free_y"), #nullGrob() #+facet_wrap(.~Component, scales = "free_y")
  p_pi,A1,hz_resultPlot,
  ncol = 4,
  top = textGrob(
    "", 
    gp = gpar(fontsize = 15, fontface = "bold", col = "Black")
  ),
  layout_matrix = rbind(
    c(1, 2, 3,3),
    c(4, 5, 6,6)  # 구분선 행
  )
)


#####################################



optimalData
optimalData %>% tail
###############################################
# Over Time 
###############################################
dim(roadResult)
# optimalData = roadResult %>% filter(tempering==0.999999999)
# optimalData = roadResult[find_last_stable_index(roadResult$diffBeta1,0.1) ,]
# optimalData = roadResult[ 15,]
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






fdata = read.table(paste0(optFilename,".txt"),header = T)
dataName = tools::file_path_sans_ext(optFilename)

# Data preprocessing
N = nrow(fdata)
k=3 
event_vec = as.numeric(fdata[,2])
time_vec = as.numeric(fdata[,1])


surv_obj <- Surv(time = time_vec, event = event_vec)
fit <- survfit(surv_obj ~ 1)

# 시간과 누적 생존율
times <- fit$time
surv_probs <- fit$surv

# 누적 hazard (대략적 추정)
cumhaz <- -log(surv_probs)

# 시간 구간별 변화량
delta_time <- diff(c(0, times))
delta_hazard <- diff(c(0, cumhaz))
hazard_rate <- delta_hazard / delta_time

# 시간 중간값 (또는 사건 발생 시간)
plot_times <- times

hazard_df <- data.frame(
  time = plot_times,
  hazard = hazard_rate
)
hazard_df = hazard_df %>% filter(hazard<Inf)
hazard_df$Component <- "Empirical hazard rate"

hz_long <- hzData %>%
  pivot_longer(cols = starts_with("hz"), names_to = "Component", values_to = "Hazard")

hz_plot_df <- bind_rows(
  hz_long, 
  hazard_df %>% rename(Hazard = hazard)
)

hz_resultPlot = ggplot(hz_plot_df, aes(x = time, y = Hazard, color = Component)) +
  geom_line(data = filter(hz_plot_df, Component != "Empirical hazard rate"), size = 1) +
  geom_point(data = filter(hz_plot_df, Component == "Empirical hazard rate"), size = 2,alpha=0.2) +
  labs(
    title = "Hazard function",
    x = "",
    y = "",
    color = "Component"
  ) +
  # geom_vline(xintercept = changePoint1, lty = 2) +
  # geom_vline(xintercept = changePoint3, lty = 2) +
  # annotate("text", x = changePoint1, y = 0.03, hjust = -0.2, label = paste0("Time: ", changePoint1)) +
  # annotate("text", x = changePoint3, y = 0.03, hjust = 1, label = paste0("Time: ", changePoint3)) +
  scale_color_manual(
    values = c(
      "hz1" = "red",
      "hz2" = "green3",
      "hz3" = "blue",
      "Empirical hazard rate" = "black"
    ),
    labels = c(
      "hz1" = "Hazard 1",
      "hz2" = "Hazard 2",
      "hz3" = "Hazard 3",
      "Empirical hazard rate" = "Empirical"
    )
  )+
  theme_minimal()+
  # coord_cartesian(ylim = c(0,2))+
  guides(color = guide_legend(ncol = 2))+
  theme(legend.position = "top",legend.title = element_blank())  # ⬅️ 범례 제거!



hz_resultPlot


# 
# DecisionPlot = posteriorPlot+hzPlot
# 
# grid.arrange(
#   DB1,DB3,
#   ncol=2,
#   top = textGrob(
#     "Posterior Ratio in the Phase of Device G Failure", 
#     gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
#   )
# )




######################################################




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


###############################################
# No DAEM 
###############################################

roadResult = read.table("Proposed_Result_Aarest_data.txt")
NoDAEMroadResult = read.table("Proposed_Result_NoDAEM_Aarest_data.txt")

roadResult = read.table("Proposed_Result_FRT_censord.txt")
NoDAEMroadResult = read.table("Proposed_Result_NoDAEM_FRT_censord.txt")

roadResult = read.table("Proposed_Result_RearDump.txt")
NoDAEMroadResult = read.table("Proposed_Result_NoDAEM_RearDump.txt")

roadResult = read.table("Proposed_Result_SerumReversal.txt")
NoDAEMroadResult = read.table("Proposed_Result_NoDAEM_SerumReversal.txt")

B1=roadResult %>% ggplot(aes(x=Iter , y=Qlike))+geom_point()+geom_line()+labs(title="log-likelihood",y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
B2=roadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
B3=roadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
B4=roadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))+geom_hline(yintercept = 1,lty=2,color="red")
B5=roadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
B10=roadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
B11=roadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
B12=roadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))

C1=NoDAEMroadResult %>% ggplot(aes(x=Iter , y=Qlike))+geom_point()+geom_line()+labs(title="log-likelihood",y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
C2=NoDAEMroadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
C3=NoDAEMroadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
C4=NoDAEMroadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))+geom_hline(yintercept = 1,lty=2,color="red")
C5=NoDAEMroadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
C10=NoDAEMroadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
C11=NoDAEMroadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))
C12=NoDAEMroadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="Iteration")+theme(axis.text = element_text(size=10), plot.title = element_text(size = 12, face = "bold"))+scale_y_continuous(labels = number_format(accuracy = 0.001))

array(0,c(2,3))

grid.arrange(
  B1, B2, B3,
  hline(),  # 첫 번째 구분선
  B4, B5,nullGrob(),
  hline(),  # 첫 번째 구분선
  B10,B11,B12,
  ncol = 3,
  top = textGrob(
    "", 
    gp = gpar(fontsize = 5, fontface = "bold", col = "Black")
  ),
  heights = c(1, 0.05, 1,0.05,1),  # 구분선 높이를 작게 설정
  layout_matrix = rbind(
    c(1, 2, 3),
    c(4, 4, 4),  # 구분선 행
    c(5, 6, 7),
    c(8, 8, 8),  # 구분선 행
    c(9, 10, 11)
  )
)


grid.arrange(
  C1, C2, C3,
  hline(),  # 첫 번째 구분선
  C4, C5,nullGrob(),
  hline(),  # 첫 번째 구분선
  C10,C11,C12,
  ncol = 3,
  top = textGrob(
    "", 
    gp = gpar(fontsize = 5, fontface = "bold", col = "Black")
  ),
  heights = c(1, 0.05, 1,0.05,1),  # 구분선 높이를 작게 설정
  layout_matrix = rbind(
    c(1, 2, 3),
    c(4, 4, 4),  # 구분선 행
    c(5, 6, 7),
    c(8, 8, 8),  # 구분선 행
    c(9, 10, 11)
  )
)




###############################################
# Hazard rate 
###############################################
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(dplyr)

source("DAEM_BarrierMethod_function.R")

# Reading Data
# file_name = 'Aarest_data.txt'
# file_name = 'FRT_censord.txt'
# file_name = 'RearDump.txt'
# file_name = 'SerumReversal.txt'
file_name = 'LFP.txt'

fdata = read.table(file_name,header = T)
dataName = tools::file_path_sans_ext(file_name)

# Data preprocessing
N = nrow(fdata)
k=3 
event_vec = as.numeric(fdata[,2])
time_vec = as.numeric(fdata[,1])


surv_obj <- Surv(time = time_vec, event = event_vec)
fit <- survfit(surv_obj ~ 1)

# 시간과 누적 생존율
times <- fit$time
surv_probs <- fit$surv

# 누적 hazard (대략적 추정)
cumhaz <- -log(surv_probs)

# 시간 구간별 변화량
delta_time <- diff(c(0, times))
delta_hazard <- diff(c(0, cumhaz))
hazard_rate <- delta_hazard / delta_time

# 시간 중간값 (또는 사건 발생 시간)
plot_times <- times
# 
# 
# initial_beta = c(0.5,1,2)
# initial_lambda = initial_lambda_calc(time_vec,event_vec,initial_beta,ratio1=0.2,ratio3=0.9)
# initial_lambda[3]=0.3
hazard_df <- data.frame(
  time = plot_times,
  hazard = hazard_rate
  # hazard1 = initial_beta[1]*initial_lambda[1]*unique(time_vec)^(initial_beta[1]-1),
  # hazard1 = initial_beta[2]*initial_lambda[2]*unique(time_vec)^(initial_beta[2]-1),
  # hazard3 = initial_beta[3]*initial_lambda[3]*unique(time_vec)^(initial_beta[3]-1)
)
hazard_df = hazard_df[1:(length(times)-1),]
hazardPoint = ggplot(hazard_df, aes(x = time, y = hazard)) +
  geom_point(color = "black", size = 2) +
  geom_line()+
  # geom_line(aes(y = hazard1), color = "red") +
  labs(title = "LFP",
       x = "Time",
       y = "Hazard Rate") +
  theme_minimal()
hz5 = hazardPoint

hz1+hz2+hz4+hz5


grid.arrange(
  hz1,hz2,hz4,hz5,
  top = textGrob(
    "Empirical Hazard Rate Plot", 
    gp = gpar(fontsize = 15, fontface = "bold", col = "Black")
  ),
  ncol = 4)

# 
# hazard_long <- hazard_df %>%
#   select(time, hazard1, hazard2 = hazard1.1, hazard3) %>%
#   pivot_longer(cols = starts_with("hazard"),
#                names_to = "Model",
#                values_to = "Hazard")
# 
# ggplot() +
#   geom_point(data = hazard_df, aes(x = time, y = hazard), color = "blue") +
#   geom_line(data = hazard_long, aes(x = time, y = Hazard, color = Model), size = 1) +
#   labs(title = "Estimated Hazard Rate with Model-Based Lines",
#        x = "Time",
#        y = "Hazard Rate") +
#   theme_minimal()



########################################################################
## TTT-plot 
##########################################

source("DAEM_BarrierMethod_function.R")

# Reading Data
# file_name = 'Aarest_data.txt'
file_name = 'FRT_censord.txt'
# file_name = 'RearDump.txt' #제외
# file_name = 'SerumReversal.txt'
# file_name = 'LFP.txt'

fdata = read.table(file_name,header = T)
dataName = tools::file_path_sans_ext(file_name)

# 고장 데이터만 정렬
failures <- sort(fdata$time[fdata$event == 1])
n <- length(failures)
T_total <- sum(failures)

# TTT 좌표 계산
x_vals <- (1:n) / n
y_vals <- sapply(1:n, function(i) {
  (sum(failures[1:i]) + (n - i) * failures[i]) / T_total
})

ttt_df <- data.frame(x = x_vals, y = y_vals)

# TTT plot
TTT2=ggplot(ttt_df, aes(x = x, y = y)) +
  geom_step(direction = "hv") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Device G",
    x = "Proportion of Failures (i/n)",
    y = "TTT Transform"
  ) +
  theme_minimal()

TTT1
TTT2
TTT3 #제외
TTT4
TTT5

TTT1+TTT2+TTT4+TTT5



grid.arrange(
  TTT1,TTT2,TTT4,TTT5,
  top = textGrob(
    "Total Time on Test (TTT) Plot", 
    gp = gpar(fontsize = 15, fontface = "bold", col = "Black")
  ),
  ncol = 4)


#########################################################












