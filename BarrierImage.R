


Gf = function(x,r=1){
  (x*exp(-x))^r
}

A=Gf(seq(1,5,length.out=1000),0.8)
B=Gf(seq(1,5,length.out=1000),0.1)
data.frame(x=seq(1,5,length.out=1000),A,B)


x=seq(0,1,length.out=1000)
fx1 = function(x){
  (-10*x^3-20*x^2+20)/200
}
bf1 = function(x){
  1/x-1/(1-x)
}
ddf = data.frame(x=x,y1=fx1(x),b1=bf1(x)) 
p1=ddf %>% ggplot(aes(x=x,y=y1))+geom_line()+geom_hline(yintercept = 0,size=1)+
  annotate(
    "text", 
    x = 0.8, y = -0.02, # 특정 좌표 지정
    label =expression(beta[1]^"*"), # 텍스트 내용
    size = 6, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+
  annotate(
    "text", 
    x = 0.5, y = -0.02, # 특정 좌표 지정
    label =expression("0.5"), # 텍스트 내용
    size = 4, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+coord_cartesian(ylim = c(-0.1,0.1))+theme_minimal()+ggtitle(expression("Derivative function "~ Q[1]))+  
  labs(
    y = "",
    x = ""
  )+  theme(
    text = element_text(family = "serif") # 폰트를 'serif'로 설정, 크기는 14
  )

p2=ddf %>% ggplot(aes(x=x,y=b1))+geom_line()+geom_hline(yintercept = 0,size=1)+coord_cartesian(ylim = c(-10,10))+
  annotate(
    "text", 
    x = 0.5, y = -1, # 특정 좌표 지정
    label =expression("0.5"), # 텍스트 내용
    size = 4, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+theme_minimal()+
  labs(
    y = "",
    x = ""
  )+ggtitle(expression("Derivative Barrier function of " ~ Q[1]))+  theme(
    text = element_text(family = "serif") # 폰트를 'serif'로 설정, 크기는 14
  )


x=seq(1,3,length.out=1000)
fx2 = function(x){
  (-1*x^3-2*x^2+20)
}
bf2 = function(x){
  1/(x-1)
}
ddf1 = data.frame(x=x,y1=fx2(x),b1=bf2(x)) 
p11=ddf1 %>% ggplot(aes(x=x,y=y1))+geom_line()+geom_hline(yintercept = 0,size=1)+
  coord_cartesian(ylim = c(-15,15))+
  theme_minimal()+
  labs(    y = "",    x = ""  )+
  annotate(
    "text", 
    x = 2.1, y = -3, # 특정 좌표 지정
    label =expression(beta[3]^"*"), # 텍스트 내용
    size = 6, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+ggtitle(expression("Derivative function "~Q[3]))+  theme(
    text = element_text(family = "serif") # 폰트를 'serif'로 설정, 크기는 14
  )
p22=ddf1 %>% ggplot(aes(x=x,y=b1))+geom_line()+geom_hline(yintercept = 0,size=1)+
  coord_cartesian(ylim = c(-30,30))+
  theme_minimal()+
  labs(    y = "",    x = ""  )+ggtitle(expression("Derivative Barrier function of " ~ Q[3]))+  theme(
    text = element_text(family = "serif") # 폰트를 'serif'로 설정, 크기는 14
  )




ddf3 = data.frame(x=x,y1=fx2(x),y2=(fx2(x)+bf2(x)*10)) 
p31=ddf3%>% ggplot(aes(x=x,y=y2))+geom_line()+geom_hline(yintercept = 0,size=1)+
  geom_line(aes(x=x,y=y1),lty=2)+
  annotate(
    "text", 
    x = 2.1, y = -5, # 특정 좌표 지정
    label =expression(beta[3]^"*"), # 텍스트 내용
    size = 5, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+
  annotate(
    "text", 
    x = 1.2, y = 25, # 특정 좌표 지정
    label =expression(~ frac(d, d~beta) ~ Q[3]), # 텍스트 내용
    size = 3, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+
  annotate(
    "text", 
    x = 2, y = 25, # 특정 좌표 지정
    label =expression(~ frac(d, d~beta) ~ BQ[3]), # 텍스트 내용
    size = 3, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+
  annotate(
    "text", 
    x = 2.5, y = 5, # 특정 좌표 지정
    label =  expression(hat(beta)[3]), # 텍스트 내용
    size = 5, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+coord_cartesian(ylim = c(-30,30))+theme_minimal()+ggtitle(expression("Derivative function "~ BQ[3]))+  
  labs(
    y = "",
    x = ""
  )+  theme(
    text = element_text(family = "serif") # 폰트를 'serif'로 설정, 크기는 14
  )


x=seq(0,1,length.out=1000)
fx1 = function(x){
  (-10*x^3-20*x^2+20)/200
}
bf1 = function(x){
  1/x-1/(1-x)
}
ddf4 = data.frame(x=x,y1=fx1(x),y2=(fx1(x)+bf1(x)/20)) 
p32=ddf4 %>% ggplot(aes(x=x,y=y2))+geom_line()+geom_hline(yintercept = 0,size=1)+
  geom_line(aes(x=x,y=y1),lty=2)+coord_cartesian(ylim = c(-0.2,0.2))+
  annotate(
    "text", 
    x = 0.5, y = -0.02, # 특정 좌표 지정
    label =expression("0.5"), # 텍스트 내용
    size = 3, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+
  annotate(
    "text", 
    x = 0.85, y = -0.05, # 특정 좌표 지정
    label =expression(beta[1]^"*"), # 텍스트 내용
    size = 5, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+
  annotate(
    "text", 
    x = 0.6, y = -0.05, # 특정 좌표 지정
    label =  expression(hat(beta)[1]), # 텍스트 내용
    size = 5, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+theme_minimal()+ggtitle(expression("Derivative function "~ BQ[1]))+  
  labs( y = "",x = "")+  theme(
    text = element_text(family = "serif") # 폰트를 'serif'로 설정, 크기는 14
  )+
  annotate(
    "text", 
    x = 0.42, y = 0.19, # 특정 좌표 지정
    label =expression(~ frac(d, d~beta) ~ BQ[1]), # 텍스트 내용
    size = 3, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )+
  annotate(
    "text", 
    x = 0.05, y = 0.15, # 특정 좌표 지정
    label =expression(~ frac(d, d~beta) ~ Q[1]), # 텍스트 내용
    size = 3, # 텍스트 크기
    color = "black", # 텍스트 색상
    fontface = "bold" # 텍스트 스타일
  )


(p1+p2+p32)/(p11+p22+p31)








