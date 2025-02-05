DD = read.table("Aarest_data_InitResult.txt",header = T)
DD = read.table("FRT_censord_InitResult.txt",header = T)
DD = read.table("RearDump_InitResult.txt",header = T)
DD = read.table("SerumReversal_InitResult.txt",header = T)




CC = read.table("Aarest_data_EMInitResult.txt",header = T)
CC = read.table("FRT_censord_EMInitResult.txt",header = T)
CC = read.table("RearDump_EMInitResult.txt",header = T)
CC = read.table("SerumReversal_EMInitResult.txt",header = T)

dim(CC)


# ggplot2 패키지를 사용하여 Qlike 값을 그룹별로 플로팅
library(ggplot2)

# 
# 
# 
# ggplot(DD, aes(x = tempering, y = beta1)) +
#   geom_point(size = 3, color = "blue",alpha=0.1) + # 포인트 추가
#   geom_line(aes(group = group,color =  group), linetype = "dashed") + # 그룹별 연결선
#   # facet_wrap(~ group, ncol = 3, scales = "free_y") + 
#   labs(
#     title = "Qlike by Groups",
#     x = "",
#     y = "Qlike"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) # X축 레이블 회전
# 
# DD$group_bin <- cut(match(DD$g, unique_values), breaks = 5, labels = c("G1", "G2", "G3", "G4", "G5"))
# 
# ggplot(DD, aes(x = tempering, y = beta1, color = as.factor(group))) +
#   geom_line() +
#   facet_wrap(~ group_bin) +  # 5개씩 묶어서 출력
#   theme_minimal()

DD$g <- paste(DD$init_beta1, DD$init_beta2, DD$init_beta3, DD$init_pi1, DD$init_pi2, DD$init_pi3, sep = "_")
unique_values <- unique(DD$g)  # 유니크한 값 추출
DD$group <- as.character( match(DD$g, unique_values) )  # 순서대로 번호 매김

DD %>% tail(10)
DDA = DD %>% filter(tempering == max(DD[which(DD$Qlike< -565),"tempering" ]))
dim( DD)
OptAnneal =DDA$tempering %>% unique
OptQlike = DDA$Qlike[1]



ggplot(DD, aes(x = tempering, y = Qlike)) +
  geom_point(size = 3, color = "blue") + # 포인트 추가
  geom_line(aes(group = group), color = "red", linetype = "dashed") + # 그룹별 연결선
  labs(
    title = "Qlike by Groups",
    x = "",
    y = "Qlike"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # X축 레이블 회전

ggplot(DD, aes(x = tempering, y = beta1)) +
  geom_point(size = 3, color = "blue") + # 포인트 추가
  geom_line(aes(group = group), color = "red", linetype = "dashed") + # 그룹별 연결선
  labs(
    title = "Qlike by Groups",
    x = "",
    y = "Qlike"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # X축 레이블 회전

Plot_DAEM = ggplot(DD, aes(x = tempering , y = Qlike)) +
  geom_point(size = 3, aes(color = group))+
  geom_point(x=OptAnneal,y=OptQlike , color="red",size=3)+
  geom_line(aes(color=group))+ labs(
    title = paste0(unique(DD$filename) ),
    x = "Annealing",
    y = "Qlike"
  ) +
  ggplot(DDA, aes(x = group , y = beta1)) +
  geom_point(size = 3, aes(color=group)) + # 포인트 추가
  geom_line(aes(group = group), color = "red", linetype = "dashed") + # 그룹별 연결선
  # facet_wrap(~ group, ncol = 3, scales = "free_y") + 
  labs(
    title = " beta1",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none") +

  ggplot(DDA, aes(x = group , y = beta3)) +
  geom_point(size = 3, aes(color=group),size=3) + # 포인트 추가
  geom_line(aes(group = group), color = "red", linetype = "dashed") + # 그룹별 연결선
  # facet_wrap(~ group, ncol = 3, scales = "free_y") + 
  labs(
    title = " beta3",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")





CC$group <- paste(CC$init_beta1, CC$init_beta2, CC$init_beta3, CC$init_pi1, CC$init_pi2, CC$init_pi3, sep = "_")
unique_values <- unique(CC$g)  # 유니크한 값 추출
CC$group <- as.character( match(CC$g, unique_values) )  # 순서대로 번호 매김
Plot_OEM =ggplot(CC, aes(x = group, y = Qlike)) +
  geom_point(size = 3, aes(color=group)) + # 포인트 추가
  labs(
    title = paste0(unique(DD$filename)," likelihood by Initial OEM"),
    x = "",
    y = "Qlike"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ # X축 레이블 회전
ggplot(CC, aes(x = group, y = beta1)) +
  geom_point(size = 3,  aes(color=group)) + # 포인트 추가
  # facet_wrap(~ group, ncol = 3, scales = "free_y") + 
  labs(
    title = " beta1",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  # geom_hline(yintercept = 1,lty=2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")+
ggplot(CC, aes(x = group, y = beta3)) +
  geom_point(size = 3,  aes(color=group)) + # 포인트 추가
  # facet_wrap(~ group, ncol = 3, scales = "free_y") + 
  labs(
    title = "beta3",
    x = "",
    y = ""
  ) +
  # geom_hline(yintercept = 1,lty=2)+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")

Plot_OEM
Plot_DAEM

library(ggpp)
DFF = DD %>% select(group,init_beta1 ,init_beta2 ,init_beta3 , init_pi1  ,init_pi2  ,init_pi3) %>% distinct %>% 
  mutate(across(starts_with("init_pi"), ~ round(.x, 2)))


ggplot() +
  geom_table(aes(x = 1, y = 1, label = list(DFF))) +
  theme_void()
