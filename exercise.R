storeRatio=numeric()
for(i in 1:length(result_latentZ_mat)){
    ddd = result_latentZ_mat[[i]]
    storeRatio[i] = ddd[1,1]/ddd[1,3]
}

plot(log(storeRatio[1:15],base=10))
ddff = data.frame(x=1:20,y=storeRatio[1:20])
p1 = ddff%>%ggplot(aes(x=x,y=log(y,10)))+geom_line()+
    geom_point()+
  labs(
    y = expression(log[10](Z[11] / Z[13])),
     x = expression(Iteration)
  )+theme_minimal() +
  geom_vline(xintercept = 15 ,lty=2,size=1,color="red")+
  theme(
    # panel.grid.major.x = element_line(color = "gray80", size = 0.5), # x축 격자만 표시
    panel.grid.major.y = element_line(color = "gray80", size = 0.5), # y축 격자 숨기기
    axis.text.y = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title = element_text(size=13)
    )+ggtitle("Figure 1")
p1
p2 = theta_df[1:20,]%>% mutate(x=1:20) %>% ggplot(aes(x=as.numeric(x),y=as.numeric(beta1)))+
geom_line()+geom_point()+
  labs(
    y = expression(beta[1]),
    x = expression(Iteration)
  )+geom_hline(yintercept = 1,color="black")+
  geom_vline(xintercept = 15 ,lty=2,size=1,color="red")+
  theme_minimal() +theme(
    # panel.grid.major.x = element_line(color = "gray80", size = 0.5), # x축 격자만 표시
    panel.grid.major.y = element_line(color = "gray80", size = 0.5), # y축 격자 숨기기
    axis.text.y = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title = element_text(size=15)
  )+ggtitle("Figure 2")

p3=theta_df %>% ggplot(aes(x=1:nrow(theta_df),y=as.numeric(beta1)))+geom_line()+
  geom_line(aes(y=as.numeric(beta3)))+
  theme_minimal()+
  labs(
    y ="",
    x = expression(Iteration)
  )+  theme_minimal() +theme(
    # panel.grid.major.x = element_line(color = "gray80", size = 0.5), # x축 격자만 표시
    panel.grid.major.y = element_line(color = "gray80", size = 0.5), # y축 격자 숨기기
    axis.text.y = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title = element_text(size=15)
  )+annotate(
    "text", 
    x = 70, y = 9.1, # 특정 좌표 지정
    label = expression(beta[1]), # 텍스트 내용
    size = 8, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+annotate(
    "text", 
    x = 70, y = 4.5, # 특정 좌표 지정
    label = expression(beta[3]), # 텍스트 내용
    size = 8, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+annotate(
    "text", 
    x = 70, y = 1.5, # 특정 좌표 지정
    label = expression(beta[2]==1), # 텍스트 내용
    size = 8, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+geom_hline(yintercept = 1)+geom_vline(xintercept = 15,lty=2,color="red")+
  annotate(
    "text", 
    x = 7.5, y = 1.5, # 특정 좌표 지정
    label = expression("Iter\n15"), # 텍스트 내용
    size = 4, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+ggtitle("Figure 3")

(p1/p2)|p3

