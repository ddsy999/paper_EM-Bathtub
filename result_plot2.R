
# roadResult = read.table("Proposed_Result_RearDump.txt")
# roadResult = read.table("Proposed_Result_FRT_censord.txt")
roadResult = read.table("Proposed_Result_LFP.txt")
# roadResult = read.table("Proposed_Result_SerumReversal.txt")

# roadResult = read.table("Proposed_Result_NoDAEM_Aarest_data.txt")

optFilename = roadResult$data_Name %>% unique()



roadResult %>% ggplot(aes(x=Iter , y=Q1+Q2))+geom_point()+geom_line()+labs(title="Expected log-likelihood",y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
# roadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))

roadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=1))    +geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
# roadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))+geom_hline(yintercept = 1,lty=2,color="red")
roadResult %>% ggplot(aes(x=Iter , y=lambda1  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=lambda2  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
# roadResult %>% ggplot(aes(x=Iter , y=lambda3  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
roadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
# roadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))




A1=roadResult %>% ggplot(aes(x=Iter , y=Q1+Q2))+geom_point()+geom_line()+labs(title="Expected log-likelihood",y="",x="Iteration")+theme(axis.text = element_text(size=12))
A2=roadResult %>% ggplot(aes(x=Iter , y=diffBeta1  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
# A3=roadResult %>% ggplot(aes(x=Iter , y=diffBeta3  ))+geom_point()+geom_line()+labs(title=expression("Gradient of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A4=roadResult %>% ggplot(aes(x=Iter , y=beta1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A5=roadResult %>% ggplot(aes(x=Iter , y=1))    +geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
# A6=roadResult %>% ggplot(aes(x=Iter , y=beta3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*beta[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A7=roadResult %>% ggplot(aes(x=Iter , y=lambda1  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A8=roadResult %>% ggplot(aes(x=Iter , y=lambda2  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
# A9=roadResult %>% ggplot(aes(x=Iter , y=lambda3  ))+geom_point()+geom_line()+labs(title=expression("Convergence of "*lambda[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A10=roadResult %>% ggplot(aes(x=Iter , y=pi1))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[1]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
A11=roadResult %>% ggplot(aes(x=Iter , y=pi2))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[2]),y="",x="Iteration")+theme(axis.text = element_text(size=12))
# A12=roadResult %>% ggplot(aes(x=Iter , y=pi3))+geom_point()+geom_line()+labs(title=expression("Convergence of "*pi[3]),y="",x="Iteration")+theme(axis.text = element_text(size=12))





# grid.arrange(
#   A1,A2,A3,
#   A4,A5,A6,
#   A7,A8,A9,
#   A10,A11,A12,
#   ncol = 3,           # 열의 개수
#   # heights = c(1, 1)     # 행 높이 비율
#   top =  textGrob(
#     optFilename, 
#     gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
#   )
# )

hline <- function() {
  grid::linesGrob(y = unit(c(0.5, 0.5), "npc"), # y축 위치 (중앙)
                  x = unit(c(0, 1), "npc"), # 왼쪽 끝에서 오른쪽 끝까지
                  gp = gpar(col = "black", lwd = 2)) # 검은색, 두께 2pt
}

grid.arrange(
  A1, A2,
  hline(),  # 첫 번째 구분선
  A4, A5, 
  hline(),  # 두 번째 구분선
  A7, A8, 
  hline(),  # 세 번째 구분선
  A10, A11, 
  ncol = 2,
  top = textGrob(
    "", 
    gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
  ),
  heights = c(1, 0.05, 1, 0.05, 1, 0.05, 1),  # 구분선 높이를 작게 설정
  layout_matrix = rbind(
    c(1, 2),
    c(3, 3),  # 구분선 행
    c(4, 5),
    c(6,6),  # 구분선 행
    c(7,8),
    c(9,9), # 구분선 행
    c(10,11)
  )
)

###############################################
# Over Time 
###############################################

optimalData = roadResult %>% filter(tempering==0.999999999) 
optFilename = roadResult$data_Name %>% unique()
optBeta1 = optimalData %>% pull(beta1)
optBeta2 = 1
# optBeta3 = optimalData %>% pull(beta3)
optLambda1 = optimalData %>% pull(lambda1)
optLambda2 = optimalData %>% pull(lambda2)
# optLambda3 = optimalData %>% pull(lambda3)
optPi1 = optimalData %>% pull(pi1)
optPi2 = optimalData %>% pull(pi2)
# optPi3 = optimalData %>% pull(pi3)

time_DBlist = seq(min(time_vec),max(time_vec),by=1)

DBbound1 = DecisionBoundary(time_DBlist,beta_vec = c(optBeta1,1,1),lambda_vec = c(optLambda1,optLambda2,1),pi_vec=c(optPi1,optPi2,1),j=1)
# DBbound3 = DecisionBoundary(time_DBlist,beta_vec = c(optBeta1,1,optBeta3),lambda_vec = c(optLambda1,optLambda2,optLambda3),pi_vec=c(optPi1,optPi2,optPi3),j=3)


changePoint1 = max(time_DBlist[(DBbound1>1)])
# changePoint3 = time_DBlist[which.min(DBbound3<1)]
DBdata = data.frame(time = time_DBlist, postProbRatio1 = DBbound1 )
DB1 = DBdata %>% ggplot(aes(x=time,y=postProbRatio1))+geom_line()+  theme_minimal()+
  geom_hline(yintercept = 1,lty=2)+
  # annotate("text", x = changePoint1+10, y = 1.3,label=paste0("Time :",changePoint1))+
  labs(
    title = "Posterior Ratio (infant vs constant)",
    x = "",
    y = ""
  )+
  geom_point(data = data.frame(x = c(changePoint1), y = c(1)), aes(x = x, y = y), color = "blue", size = 3)
  # coord_cartesian(xlim =c(changePoint1-0.1,changePoint1+0.1),ylim = c(0,2) )

posteriorPlot = DB1

hzData = data.frame(time = time_DBlist , hz1 = hazardrate(time_DBlist,optBeta1,optLambda1),
                    hz2 = hazardrate(time_DBlist,optBeta2,optLambda2))

hzPlot = hzData %>% ggplot(aes(x=time , y=hz1))+geom_line(color="red")+
  geom_line(aes(x=time ,y=hz2),color="green3")+
  # scale_y_break(c(0.3,2), space = 0.3 , scales="free") +
  geom_vline(xintercept = changePoint1,lty=2)+
  # geom_vline(xintercept = changePoint3,lty=2)+
  labs(
    title = "Hazard Rate",
    x = "",
    y = ""
  )+
  # annotate("text", x = changePoint1, y = 0.01,hjust=-0.2,label=paste0("Time :",changePoint1))+
  # annotate("text", x = changePoint3, y = 0.01,hjust=-0.3,label=paste0("Time :",changePoint3))+
  # coord_cartesian(xlim = c(3000,4000))+
  theme_minimal()

hzPlot

DecisionPlot = posteriorPlot+hzPlot

grid.arrange(
  DB1,
  ncol=1,
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

fdata = read.table(file_name,header = T)
dataName = tools::file_path_sans_ext(file_name)

# Data preprocessing
N = nrow(fdata)
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
hazardOny_df <- data.frame(
  time = plot_times,
  hazard = hazard_rate
)

hazardOnlyPoint = ggplot(hazardOny_df, aes(x = time, y = hazard)) +
  geom_point(color = "blue", size = 2) +
  labs(title = "Estimated Hazard Rate",
       x = "Time",
       y = "Hazard Rate") +
  theme_minimal()
hazardOnlyPoint

initial_beta = c(0.1,1)
initial_lambda = initial_lambda_func2(time_vec,event_vec,initial_beta,ratio1=0.01,ratio3=0.9)
# initial_lambda[3]=0.3
hazard_df <- data.frame(
  time = plot_times,
  hazard = hazard_rate,
  hazard1 = initial_beta[1]*initial_lambda[1]*unique(time_vec)^(initial_beta[1]-1),
  hazard2 = initial_beta[2]*initial_lambda[2]*unique(time_vec)^(initial_beta[2]-1)
)

hazardPoint = ggplot(hazard_df, aes(x = time, y = hazard)) +
  geom_point(color = "blue", size = 2) +
  geom_line(aes(y = hazard1), color = "red") +
  geom_line(aes(y = hazard2), color = "green") +
  labs(title = "Estimated Hazard Rate",
       x = "Time",
       y = "Hazard Rate") +
  theme_minimal()

hazardPoint

hazard_long <- hazard_df %>%
  select(time, hazard1, hazard2 = hazard1.1, hazard3) %>%
  pivot_longer(cols = starts_with("hazard"),
               names_to = "Model",
               values_to = "Hazard")

ggplot() +
  geom_point(data = hazard_df, aes(x = time, y = hazard), color = "blue") +
  geom_line(data = hazard_long, aes(x = time, y = Hazard, color = Model), size = 1) +
  labs(title = "Estimated Hazard Rate with Model-Based Lines",
       x = "Time",
       y = "Hazard Rate") +
  theme_minimal()
