library(ggplot2)
library(dplyr)

# 시간(Time) 범위 설정
time <- seq(0, 100, length.out = 50000)

# Hazard Rate 함수 정의 (욕조 곡선 형태)
hazard_rate <- function(t) {
  ifelse(t < 20, 
         0.09 + 0.9 * exp(-0.2 * t),               # 초기 고장 (burn-in)
         ifelse(t < 80, 
                rep(0.1, length(t)),              # 정상 수명
                0.1 + 0.002 * (t - 80)^2))        # 마모 고장 (wear-out)
}

# 데이터 프레임 생성
df <- data.frame(
  time = time,
  hazard = hazard_rate(time)
)

# ggplot으로 그리기
ggplot(df, aes(x = time, y = hazard)) +
  geom_line(size = 1) +
  geom_vline(xintercept = c(20, 80), linetype = "dashed") +
  annotate("text", x = 10, y = 0.8, label = "Early\n(burn-in)", size = 6) +
  annotate("text", x = 50, y = 0.18, label = "Constant\n(useful, random)", size = 6) +
  annotate("text", x = 88, y = 0.8, label = "Wear-out", size = 6) +
  annotate("segment", x = 20, xend = 80, y = 0.9, yend = 0.9,
           arrow = arrow(ends = "both", length = unit(0.2, "cm")),
           size = 1) +
  annotate("text", x = 50, y = 0.95, label = "Useful life", size = 6) +
  annotate("segment", x = 0, xend = 20, y = 0.9, yend = 0.9,
           arrow = arrow(ends = "both", length = unit(0.2, "cm")),
           size = 1) +
  annotate("text", x = 10, y = 0.95, label = "Burn-in Process", size = 5) +
  labs(
    title = "                                                     Bathtub Curve",
    x = "Time/Cycles",
    y = "Hazard Rate"
  ) +
  theme_minimal(base_size = 14)+
  theme(
    axis.text = element_blank(),        # 눈금 텍스트 제거
    axis.ticks = element_blank(),       # 눈금선 제거
    panel.grid = element_blank(),       # 배경 격자 제거
    panel.border = element_blank()
  ) +
  geom_hline(yintercept = 0, size = 0.5) +
  geom_vline(xintercept = 0, size = 0.5)+
  geom_point(x = 20, y = 0.1, size = 5,alpha=0.1, color = "red")+
  geom_point(x = 80, y = 0.1, size = 5,alpha=0.1, color = "blue")
