library(ggplot2)
library(dplyr)

# --- 필요한 벡터들: 사용자 제공 ---
# 예: 10개의 데이터 포인트
DHIter = 20

latentZ_mat <- result_latentZ_mat[[DHIter]]
pi_vec <- theta_df_full[DHIter,] %>% select(pi1,pi2,pi3)
bpBase = theta_df_full[DHIter,] %>% pull(bpBase)
bpBase = 1
theta_df_full[,"lambda1"]
tail(theta_df_full)

# start_beta1 =  theta_df_full[DHIter-1,] %>% pull(beta1)
# start_beta3 =  theta_df_full[DHIter-1,] %>% pull(beta3)

start_beta1 =  0.5
start_beta3 =  2

# --- 사용자 정의 Q 함수 ---
Qfunc <- function(beta = 1, lambda, latentZ_mat, j = 1) {
  sum(latentZ_mat[, j] * event_vec * log(lambda)) +
    sum(latentZ_mat[, j] * event_vec * log(beta)) +
    sum(latentZ_mat[, j] * event_vec * (beta - 1) * log(time_vec)) -
    sum(latentZ_mat[, j] * lambda * time_vec^beta) +
    sum(latentZ_mat[, j] * log(pi_vec[j]))
}

Qfunc_contour <- function(beta1, beta3, latentZ_mat) {
  new_beta <- c(beta1, NA, beta3)
  new_lambda <- sapply(c(1, 3), function(i) {
    sum(latentZ_mat[, i] * event_vec) / sum(latentZ_mat[, i] * (time_vec^new_beta[i]))
  })
  Qfunc(beta1, new_lambda[1], latentZ_mat, j = 1) +
    Qfunc(beta3, new_lambda[2], latentZ_mat, j = 3)
}


Qfunc_contour_with_Barrier <- function(beta1, beta3, latentZ_mat,bpBase) {
  new_beta <- c(beta1, NA, beta3)
  new_lambda <- sapply(c(1, 3), function(i) {
    sum(latentZ_mat[, i] * event_vec) / sum(latentZ_mat[, i] * (time_vec^new_beta[i]))
  })
  Qfunc(beta1, new_lambda[1], latentZ_mat, j = 1) +(1/bpBase)*log(beta1*(1-beta1))+
    Qfunc(beta3, new_lambda[2], latentZ_mat, j = 3)+(1/bpBase)*log((beta3-1))
}



# beta1, beta3 그리드 생성
beta1_seq <- seq(0.01, 1.5, length.out = 200)
beta3_seq <- seq(0.01, 6, length.out = 100)

grid_df <- expand.grid(beta1 = beta1_seq, beta3 = beta3_seq)

# Q값 계산
grid_df$Qvalue <- mapply(function(b1, b3) Qfunc_contour(b1, b3, latentZ_mat),
                         grid_df$beta1, grid_df$beta3)

grid_df$QvalueBarrier <- mapply(
  function(b1, b3) Qfunc_contour_with_Barrier(b1, b3, latentZ_mat, bpBase),
  grid_df$beta1, grid_df$beta3
)

restricted_df <- data.frame(
  xmin = c(1, min(grid_df$beta1)),
  xmax = c(max(grid_df$beta1), 1),
  ymin = c(min(grid_df$beta3), min(grid_df$beta3)),
  ymax = c(max(grid_df$beta3), 1)
)

# 최대점
max_point <- grid_df %>% filter(Qvalue == max(Qvalue))
max_pointBp <- grid_df %>% filter(QvalueBarrier == max(QvalueBarrier, na.rm = TRUE))
start_point_df = data.frame(beta1=start_beta1,beta3=start_beta3,Qvalue=1,QvalueBarrier=1)
# 등고선 그리기
percent_val = DHIter

q20=ggplot(grid_df, aes(x = beta1, y = beta3, z = Qvalue)) +
  geom_contour(binwidth = 5, alpha = 0.5) +
  geom_rect(data = restricted_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey60", alpha = 0.3) +
  geom_point(data = max_point, aes(x = beta1, y = beta3),
             color = "red", size = 3) +
  geom_text(data = max_point, aes(label = "Max", x = beta1, y = beta3),
            hjust = -0.2, vjust = -0.5, color = "red", fontface = "bold", size = 4) +
  # geom_point(data = start_point_df, aes(x = beta1, y = beta3),
  #            color = "darkgreen", size = 3, alpha = 0.8) +
  # geom_text(data = max_point, aes(label = "Initial", x = start_beta1, y = start_beta3),
  #           hjust = -0.2, vjust = -0.5, color = "darkgreen", fontface = "bold", size = 4) +
  labs(
    title = bquote("(A) "*"Q"*"(β"*")  Contour ("* .(percent_val) *"% D-H Iteration)"),
    x = expression(beta[1]),
    y = expression(beta[3])
  ) +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold")
  )



b20 = ggplot(grid_df, aes(x = beta1, y = beta3, z = QvalueBarrier)) +
  geom_contour(binwidth = 5, alpha = 0.5) +
  geom_rect(data = restricted_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey60", alpha = 0.3) +
  geom_point(data = max_pointBp, aes(x = beta1, y = beta3),
             color = "red", size = 3) +
  geom_text(data = max_pointBp, aes(label = "Max", x = beta1, y = beta3),
            hjust = -0.2, vjust = -0.5, color = "red", fontface = "bold", size = 4) +
  # geom_point(data = start_point_df, aes(x = beta1, y = beta3),
  #            color = "darkgreen", size = 3, alpha = 0.8) +
  # geom_text(data = max_point, aes(label = "Initial", x = start_beta1, y = start_beta3),
  #           hjust = -0.2, vjust = -0.5, color = "darkgreen", fontface = "bold", size = 4) +
  labs(
    title = bquote("(C) "*"Q"*"(β"*")  Contour With Barrier ("* .(percent_val) *"% D-H Iteration)"),
    x = expression(beta[1]),
    y = expression(beta[3])
  ) +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold")
  )



(q20+q90)/(b20+b90)
