library(reshape2)
library(dplyr)
library(ggplot2)
library(patchwork)
library(survival)
library(ggbreak)
library(cowplot)
library(gridExtra)
library(grid)
library(scales)
library(foreach)
library(doParallel)

# DFdata = read.table("Aarest_data_DAEM.txt",header = T)
# DFdata = read.table("FRT_censord_DAEM.txt",header = T)
DFdata = read.table("RearDump_DAEM.txt",header = T)
DFdata %>% ggplot(aes(x=tempering,y=Qlike))+geom_line(lty=2)+geom_point()

appAnnealLimit = DFdata$tempering[max(which(DFdata$Qlike < -150))]
optimalData = DFdata %>% filter(tempering==appAnnealLimit)
OptPi1 = optimalData$pi1
OptPi2 = optimalData$pi2
OptPi3 = optimalData$pi3
OptBeta1 = optimalData$beta1
OptBeta3 = optimalData$beta3
OptLike = optimalData$Qlike
OptDiffBeta1 =  signif(optimalData$diffBeta1, 3) 
OptDiffBeta3 =  signif(optimalData$diffBeta3, 3) 



P1 = DFdata %>% ggplot(aes(x=tempering,y=Qlike))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(A) Likelihood Q("* beta*")"),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_point(aes(x=appAnnealLimit,y=OptLike),color="red",size=3)+
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  annotate("text", x = appAnnealLimit, y = OptLike-2, label = appAnnealLimit, color = "red",
           hjust = 1.2,size=3)+
  theme_minimal()
P1 
P2 = DFdata %>% ggplot(aes(x=tempering,y=beta1))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(B) Convergence of " * beta[1]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_point(aes(x=appAnnealLimit,y=OptBeta1),color="red",size=3)+
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()
P2
P3 = DFdata %>% ggplot(aes(x=tempering,y=beta3))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(C) Convergence of " * beta[3]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_point(aes(x=appAnnealLimit,y=OptBeta3),color="red",size=3)+
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  theme_minimal()
P3 

P4 = DFdata %>% ggplot(aes(x=tempering,y=diffBeta1))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(D) Derivative of " * beta[1]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",lty=3,size=1)+
  geom_point(aes(x=appAnnealLimit,y=OptDiffBeta3),color="red",size=3)+
  annotate("text", x = appAnnealLimit-0.1, y = OptDiffBeta1+0.01, label =OptDiffBeta1, size = 4) +
  theme_minimal()
P4

P5 = DFdata %>% ggplot(aes(x=tempering,y=diffBeta3))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(E) Derivative of " * beta[3]),
    x = "Annealing parameter",
    y = ""
  ) +
  geom_vline(xintercept = appAnnealLimit,color="red",linetype=3,size=1)+
  geom_point(aes(x=appAnnealLimit,y=OptDiffBeta3),color="red",size=3)+
  annotate("text", x = appAnnealLimit-0.2, y = OptDiffBeta3-1e-6, label =OptDiffBeta3, size = 4) +
  theme_minimal()+
  theme(axis.text.y=element_text(size=10))
P5

title_grob <- textGrob(
  "Rear Dump - DAEM with Barrier Method", 
  gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
)

grid.arrange(
  P1,P2,P3,P4,P5,
  ncol = 5 ,        # 열의 개수
  # heights = c(1, 1)     # 행 높이 비율
  top = title_grob
)

#########################################

# DFdataOEM = read.table("Aarest_OEM.txt",header = T)
# DFdataOEM = read.table("FRT_OEM.txt",header = T)
DFdataOEM = read.table("RearDump_OEM.txt",header = T)
DFdataOEM %>% tail

C1 = DFdataOEM %>% ggplot(aes(x=iter,y=sumQfunc  ))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(A) Likelihood Q("* beta*")"),
    x = "Iteration",
    y = ""
  ) +
  theme_minimal()

C2 = DFdataOEM %>% ggplot(aes(x=iter,y=beta1 ))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(B) Convergence of " * beta[1]),
    x = "Iteration",
    y = ""
  ) + theme_minimal()

C3 = DFdataOEM %>% ggplot(aes(x=iter,y=beta3 ))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(C) Convergence of " * beta[3]),
    x = "Iteration",
    y = ""
  ) + theme_minimal()


C4 = DFdataOEM %>% ggplot(aes(x=iter,y=diffB_beta1 ))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(D) Derivative of " * beta[1]),
    x = "Iteration",
    y = ""
  ) + theme_minimal()


C5 = DFdataOEM %>% ggplot(aes(x=iter,y=diffB_beta3 ))+geom_line(lty=2)+geom_point()+
  labs(
    title = expression("(E) Derivative of " * beta[3]),
    x = "Iteration",
    y = ""
  ) + theme_minimal()


title_grob <- textGrob(
  "Rear Dump - Original EM", 
  gp = gpar(fontsize = 16, fontface = "bold", col = "Black")
)

grid.arrange(
  C1,C2,C3,C4,C5,
  ncol = 5 ,        # 열의 개수
  # heights = c(1, 1)     # 행 높이 비율
  top = title_grob
)

