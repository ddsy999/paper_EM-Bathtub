
file_name = 'dieselGeneratorFans.txt'
file_name = 'Aarest_data.txt'
file_name = 'RearDump.txt'
file_name = 'SerumReversal.txt'
# file_name = 'FRT_censord.txt'
# file_name = 'throttle.txt'
fdata = read.table(file_name,header = T)
event_vec = fdata[,2] %>% as.numeric()
time_vec = fdata[,1]%>% as.numeric()

surv_obj <- Surv(time_vec, event_vec)
# Kaplan-Meier 모델 적합
km_fit <- survfit(surv_obj ~ 1)

# Nelson-Aalen 추정 (누적 위험 함수 추정)
na_fit <- basehaz(coxph(surv_obj ~ 1), centered = FALSE)

# 위험률(hazard rate) 계산: 누적 위험 변화율
hazard_rate <- diff(na_fit$hazard) / diff(na_fit$time)
hazard_time <- na_fit$time[-1]

# 데이터프레임 생성
hazard_data <- data.frame(
  time = hazard_time,
  hazard = hazard_rate
)

# 플롯 그리기
ggplot(hazard_data, aes(x = time, y = hazard)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(
    title = "Hazard Rate Point Plot (Derived from KM/Nelson-Aalen)",
    x = "Time",
    y = "Hazard Rate"
  ) +
  theme_minimal()












# TTT 계산 함수
calculate_TTT <- function(failure_times, total_items) {
  # 정렬된 고장 시간
  sorted_times <- sort(failure_times)
  
  # TTT 계산을 위한 초기화
  TTT_values <- numeric(length(sorted_times))
  
  # TTT 계산
  for (i in seq_along(sorted_times)) {
    if (i == 1) {
      TTT_values[i] <- sorted_times[i] * total_items
    } else {
      TTT_values[i] <- TTT_values[i - 1] + 
        (sorted_times[i] - sorted_times[i - 1]) * (total_items - (i - 1))
    }
  }
  
  # 데이터 프레임 반환
  data.frame(
    Failure_Time = sorted_times,
    TTT = TTT_values,
    Proportion = seq_along(sorted_times) / total_items
  )
}




# 주어진 데이터
total_items <- length(time_vec) # 총 부품 수

# TTT 계산
TTT_data <- calculate_TTT(time_vec, total_items)

TTT_data %>% ggplot(aes(x=Proportion,y=TTT/max(TTT_data$TTT)))+geom_line()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")








